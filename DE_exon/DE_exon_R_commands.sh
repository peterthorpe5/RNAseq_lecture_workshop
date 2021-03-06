#R
# further reading
# https://bioconductor.org/packages/release/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.html
#to install the DEXSEQ package:

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DEXSeq")

# set the working directory to where you want:
# in Rstudio you can do this easily with File, new project
# note: windows users still need the path with / characters
setwd("/Users/pjt6/Documents/RNAseq_workshop/DE_exon")

# Make sure all the files are decompressed. (gunzip file.gz)

# check the current working directory
getwd()

# load the R package that does the DE exon analysis
library("DEXSeq")

inDir = '/Users/pjt6/Documents/RNAseq_workshop/DE_exon'

# this will load all files in the specified directory that end in 
# .exon.counts. Make sure they are decompressed. 
countFiles = list.files(inDir, pattern=".exon.counts$", full.names=TRUE)

# show is the files that it found
basename(countFiles)

# this is the gff file, (genome feature files) which contains the exon - to gene info
# already prepared as the program likes it using their scripts. 
flattenedFile = list.files(inDir, pattern="gff$", full.names=TRUE)

basename(flattenedFile)


# create a sample table of our conditons and the replicas
# pairend is the library prep type.			
sampleTable = data.frame(row.names = c("Mc_PR_Cherry_1", "Mc_PR_Cherry_2", "Mc_PR_Cherry_3", "Mc_PR_Galium_1", "Mc_PR_Galium_2", "Mc_PR_Galium_3", "Mc_PR_cress_1", "Mc_PR_cress_2", "Mc_PR_cress_3"), condition = c("Mc_PR_Cherry", "Mc_PR_Cherry", "Mc_PR_Cherry", "Mc_PR_Galium", "Mc_PR_Galium", "Mc_PR_Galium", "Mc_PR_cress", "Mc_PR_cress", "Mc_PR_cress"),libType = c("paired-end", "paired-end", "paired-end", "paired-end", "paired-end", "paired-end", "paired-end", "paired-end", "paired-end"))
			
# see the sample table
sampleTable

######################################################################################################

# this is a more complex experimental design where the exon and sample (plant host type) is compared
# (condition:exon is an interation term under stood in R statistics). 
# DEXSeqDataSetFromHTSeq is an object which will be set up to hold all the data for susequent steps
dxd = DEXSeqDataSetFromHTSeq(countFiles,sampleData=sampleTable,design= ~ sample + exon + condition:exon,flattenedfile=flattenedFile)

# create a file with gene name, grep genome.gff | cut -f gene_name > geneIDsinsubset.txt
# if you want to run this on a smaller number of genes then you can control the list here
# note the gene names were change din the gff, which is why this doesnt work. 
# can get them from the gff if you wish?
genesForSubset = read.table(file.path(inDir, "geneIDsinsubset.txt"),stringsAsFactors=FALSE)[[1]]

# This will fail, dont worry! just continue
dxd = dxd[geneIDs( dxd ) %in% genesForSubset,]

colData(dxd)


# lots at the first 5 entries
# counts per condition per exon, per gene
head( counts(dxd), 5)
# samples explicitly stated
head( featureCounts(dxd), 5) 
# extra info
head ( rowRanges(dxd), 3 )

sampleAnnotation( dxd )

####################################################################################################################################################################################
4) normalisation

# can use multiple core if using  package BiocParallel,  and  implemented  the BPPARAM

# This may not work on Windows
BiocManager::install("BiocParallel")
library(BiocParallel)
BPPARAM= MulticoreParam(workers=4)

# on windows - this doesnt work!!
#SnowParam(workers=4, type="SOCK")

# uses DEseq2 method (control for the different lib sizes - sequencing depths)
# This is essential in all DE analysis. If you have sample A with HUGE cov
# then this would alter the relative, but not the actual results, thus the data needs to be normalised according
# to this relative depths of each sequencing library
dxd = estimateSizeFactors( dxd )


# dispersion estimation -  estimate the technical variation - THIS IS SUPER SLOW
# This is essential in all DE analysis as this models the variation, of what it
# think is technical variation and thus needs to be accounted for (versus real
# biological data.
# This is SLOOOOOW 
# in the backgroun it samples data points (normally genes) which do not differ
# by large amounts. Thus it avoid genes that differ by large amounts. Over a number
# of sampling sets the variation between replicas and the whole dataset can be estimated.
# Thus the more reps the greater the ability to estimate the techincal variation and thus
# account for this in the real model
HOW="The fitting proceeds as follows: for each gene, an estimate of the dispersion is 
found which maximizes the Cox Reid-adjusted profile likelihood 
(the methods of Cox Reid-adjusted profile likelihood maximization for estimation 
of dispersion in RNA-Seq data were developed by McCarthy, et al. (2012); 
a trend line capturing the dispersion-mean relationship 
 is fit to the maximum likelihood estimates; a normal prior is determined for the log 
 dispersion estimates centered on the predicted value from the trended fit with variance 
 equal to the difference between the observed variance of the log dispersion estimates 
 and the expected sampling variance; finally maximum a posteriori dispersion estimates 
 are returned. This final dispersion parameter is used in subsequent tests. The final 
 dispersion estimates can be accessed from an object using dispersions. The 
 fitted dispersion-mean relationship is also used in varianceStabilizingTransformation. 
 All of the intermediate values (gene-wise dispersion estimates, fitted dispersion 
 estimates from the trended fit, etc.) are stored in mcols(dds), with information 
 about these columns in mcols(mcols(dds)).
The log normal prior on the dispersion parameter has been proposed by Wu, et al. (2012) 
and is also implemented in the DSS package."


dxd = estimateDispersions( dxd)

# multi core -  you can run this in multi threaded mode. - Never worked for me!!
dxd = estimateDispersions( dxd, BPPARAM=BPPARAM)

# this graph represents per-exon dispertion estimates vs the mean normalised count.
# Then the fitted line is added which has been used to normalise the data set for 
# technical variation. 
# In more detail: plots the per-exon dispersion estimates versus the mean normalised count, 
# the resulting fitted values and the a posteriori (shrinked) dispersion estimates
plotDispEsts( dxd )

# Now we have worked out the dipersion of the data (estimates of ...), technical variation,
# we have worked ou the size factors which is a result of different sequecing depths of 
# the raw data. A generalised linear model is used to fit the data. 
# test for DE (exon) expression for each exon in each gene with the 
# formula ~sample + exon + condition:exon and compare it to the 
# smaller model (the null model) ~ sample + exon
dxd = testForDEU( dxd )

		# multi core
		dxd = testForDEU( dxd, BPPARAM=BPPARAM)



# estimate fold change per exon based on the condition. This is based on the formula and the GLM. 
dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition")

# get the results and store in a new variable called dxr1 -  we dont want to write over the old 
# variable name
dxr1 = DEXSeqResults( dxd )

# see the datadxr1

dxr1

# get a description of each coloumn:
mcols(dxr1)$description

# write ou the full results table
write.table(dxr1, file = "Dexseq_results.out", append = FALSE, quote = TRUE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")


#From this object, we can ask how many genes are significant with a false discovery rate of 5%:
table ( dxr1$padj < 0.05 )

#From this object, we can ask how many genes are significant with a false discovery rate of 1%:
table ( dxr1$padj < 0.01 )

# ask how many genes are affected (FDR 0.001)
table ( tapply( dxr1$padj < 0.001, dxr1$groupID, any ) )

# MA lot: plot the log of fold change vs the average normalised count per exon. Red are exons
# considered significant. 
plotMA(dxr1, cex=0.8)



####################################################################################################################################################################################
4) visualization
# put in the appropriate gene name of interest!!!!!!!

# to plot the results for one specific gene of interest ....
plotDEXSeq( dxr1, "1056", displayTranscripts=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )


# WRITE TO BROWSEABLE HTML USING DEXSeqHTML - THIS IS THE COOL BIT. right at the end :)
# this will plot out results for all genes. Delete it when you are finished with it (end of module)
# as you will have many thousands of files here which will slow your comp down. 
DEXSeqHTML( dxr1, FDR=0.001)

# go into the "DEXSeqReport" folder. Open the html with firefox - GO WILD!