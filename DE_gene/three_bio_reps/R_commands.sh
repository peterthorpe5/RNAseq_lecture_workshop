# PLEASE note this is not a shell script but the colours when the file ends in 
# .sh and notepad++ agree with me :)


# load the package
# if not installed
install.packages("BiocManager")
BiocManager::install("edgeR")

# load the package
library("edgeR")

# set the working directory to where you want:
# in Rstudio you can do this easily with File, new project
setwd("/Users/pjt6/Documents/RNAseq_workshop/DE_gene/three_bio_reps")

# check it
getwd()

# see what is in the directory 
dir()

# load in the data
data <- read.delim("Mcerasi.genes.counts.matrix", header=T, row.names=1)

# group the replicas
group <- factor(c(1,1,1,2,2,2))

# store the data in a list-based object.
rnaseqMatrix <- DGEList(counts=data,group=group)

# have a little look at the data 
head(rnaseqMatrix)

# filter very low expression genes as these do not contribute and negatively afffect the stats
keep <- filterByExpr(rnaseqMatrix)

rnaseqMatrix <- rnaseqMatrix[keep,,keep.lib.sizes=FALSE]

# to account for sequencing depth, calcNormFactors finds a set of scaling facotrs for lib sizes.
# this minimised the log fold change between samples for most genes
# Note this is not FRPM or TPM normalisation, raw values need to be given to EdgeR, as these are
# needed to estimate the mean-variance relationship between the samples
rnaseqMatrix <- calcNormFactors(rnaseqMatrix)

# write a table of the lib size and normalisation factors. Look at how these are different. 
rnaseqMatrix$samples$eff.lib.size = rnaseqMatrix$samples$lib.size * rnaseqMatrix$samples$norm.factors
write.table(rnaseqMatrix$samples, file="example.matrix.TMM_info.txt", quote=F, sep="\t", row.names=F)

# have a look
rnaseqMatrix$samples

# group are your samples
design <- model.matrix(~group)

rnaseqMatrix <- estimateDisp(rnaseqMatrix,design)

# To perform quasi-likelihood F-tests: (better for low numbers of reps)
fit <- glmQLFit(rnaseqMatrix,design)

# have a little look
fit

qlf <- glmQLFTest(fit,coef=2)

topTags(qlf)


tTags = topTags(qlf,n=NULL)

result_table = tTags$table

result_table = data.frame(sampleA="cherry", sampleB="gallium", result_table)

result_table$logFC = -1 * result_table$logFC

write.table(result_table, file='M.cerasi_cherry_vs_gallium.GLM.edgeR.DE_results', sep='	', quote=F, row.names=T)

write.table(rnaseqMatrix, file='M.cerasi_cherry_vs_gallium.GLM.edgeR.count_matrix', sep='	', quote=F, row.names=T)


#To perform likelihood ratio tests: (better with GLM)
fit <- glmFit(rnaseqMatrix,design)

lrt <- glmLRT(fit,coef=2)

topTags(lrt)

#### classicL EDGR

conditions = factor(c(rep("cherry", 3), rep("gallium", 3)))

exp_study = DGEList(counts=rnaseqMatrix, group=conditions)

exp_study = calcNormFactors(exp_study)

exp_study = estimateCommonDisp(exp_study)

exp_study = estimateTagwiseDisp(exp_study)

# this is the stats test used to compare these using 1 factor. Anymore you need to use a GLM
et = exactTest(exp_study, pair=c("cherry", "gallium"))

tTags = topTags(et,n=NULL)

result_table = tTags$table

result_table = data.frame(sampleA="cherry", sampleB="gallium", result_table)

result_table$logFC = -1 * result_table$logFC

write.table(result_table, file='M.cerasi_cherry_vs_gallium.edgeR.DE_results', sep='	', quote=F, row.names=T)

write.table(rnaseqMatrix, file='M.cerasi_cherry_vs_gallium.edgeR.count_matrix', sep='	', quote=F, row.names=T)
# some functions:

plot_MA = function(logCounts, logFoldChange, FDR, xlab="logCounts", ylab="logFC", title="MA plot", pch=20) {

    plot(logCounts, logFoldChange, col=ifelse(FDR<=0.05, "red", "black"), xlab=xlab, ylab=ylab, main=title, pch=pch);;

}


plot_Volcano = function(logFoldChange, FDR, xlab="logFC", ylab="-1*log10(FDR)", title="Volcano plot", pch=20) {

   plot(logFoldChange, -1*log10(FDR), col=ifelse(FDR<=0.05, "red", "black"), xlab=xlab, ylab=ylab, main=title, pch=pch);

}


plot_MA_and_Volcano = function(logCounts, logFoldChange, FDR, xlab="logCounts", ylab="logFC", title="MA plot") {

    def.par = par(no.readonly = TRUE) # save default, for resetting...

    gridlayout = matrix(c(1:4),nrow=2,ncol=2, byrow=TRUE);
    layout(gridlayout, widths=c(1,1,1,1), heights=c(1,1,1,1)) 

    #plot_MA(logCounts, logFoldChange, FDR);
    #plot_Volcano(logFoldChange, FDR);

    # draw again, but use a smaller dot for data points
    plot_MA(logCounts, logFoldChange, FDR, pch='.');
    plot_Volcano(logFoldChange, FDR, pch='.');
    

    par(def.par)   
        
    
}

####
pdf("Volcano.pdf")

plot_MA_and_Volcano(result_table$logCPM, result_table$logFC, result_table$FDR)

dev.off()
# load the sample types

samples_data = read.table("samples_described.txt", header=F, check.names=F, fill=T)
samples_data = samples_data[samples_data[,2] != '',]
sample_types = as.character(unique(samples_data[,1]))
rep_names = as.character(samples_data[,2])
data = data[, colnames(data) %in% samples_data[,2], drop=F ]
nsamples = length(sample_types)
sample_colors = rainbow(nsamples)
names(sample_colors) = sample_types
sample_type_list = list()
for (i in 1:nsamples) {
    samples_want = samples_data[samples_data[,1]==sample_types[i], 2]
    sample_type_list[[sample_types[i]]] = as.vector(samples_want)
}
sample_factoring = colnames(data)
for (i in 1:nsamples) {
    sample_type = sample_types[i]
    replicates_want = sample_type_list[[sample_type]]
    sample_factoring[ colnames(data) %in% replicates_want ] = sample_type}

# load in the heat map function


primary_data = read.table("Mcerasi.genes.counts.matrix", header=T, com='', sep="\t", row.names=1, check.names=F)
primary_data = as.matrix(primary_data)
initial_matrix = as.matrix(primary_data)


# to do a PCA
pdf("minRow10.CPM.log2.principal_components.pdf")
data = as.matrix(primary_data)
# Z-scale and center the genes across all the samples for PCA
prin_comp_data = initial_matrix
prin_comp_data = log2(prin_comp_data+1)
prin_comp_data = scale(prin_comp_data)
prin_comp_data = t(scale(t(prin_comp_data), center=TRUE, scale=F)) # just center trans expr level, retain original effect size.
pca = prcomp(t(prin_comp_data), center = FALSE, scale. = FALSE)
pc_pct_variance = (pca$sdev^2)/sum(pca$sdev^2)
def.par <- par(no.readonly = TRUE) # save default, for resetting...
gridlayout = matrix(c(1:4),nrow=2,ncol=2, byrow=TRUE);
layout(gridlayout, widths=c(1,1));
write.table(pca$rotation, file="minRow10.CPM.log2.PCA.loadings", quote=F, sep="	")
write.table(pca$x, file=".minRow10.CPM.log2.PCA.scores", quote=F, sep="	")
for (i in 1:(max(3,2)-1)) {
     xrange = range(pca$x[,i])
     yrange = range(pca$x[,i+1])
     samples_want = rownames(pca$x) %in% sample_type_list[[sample_types[1]]]
     pc_i_pct_var = sprintf("(%.2f%%)", pc_pct_variance[i]*100)
     pc_i_1_pct_var = sprintf("(%.2f%%)", pc_pct_variance[i+1]*100)
     plot(pca$x[samples_want,i], pca$x[samples_want,i+1], xlab=paste('PC',i, pc_i_pct_var), ylab=paste('PC',i+1, pc_i_1_pct_var), xlim=xrange, ylim=yrange, col=sample_colors[1])
     for (j in 2:nsamples) {
         samples_want = rownames(pca$x) %in% sample_type_list[[sample_types[j]]]
         points(pca$x[samples_want,i], pca$x[samples_want,i+1], col=sample_colors[j], pch=j)
     }
     plot.new()
     legend('topleft', as.vector(sample_types), col=sample_colors, pch=1:nsamples, ncol=2)
 }
 
 par(def.par)
 pcscore_mat_vals = pca$rotation[,1:3]
 pcscore_mat = matrix_to_color_assignments(pcscore_mat_vals, col=colorpanel(256,'purple','black','yellow'), by='row')
 colnames(pcscore_mat) = paste('PC', 1:ncol(pcscore_mat))
 dev.off()
 
