
# can use this to run de with multi conditions. 
# RNentropy: an entropy-based tool for the detection of significant variation of gene expression across multiple RNA-Seq experiments	
https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gky055/4829696
http://159.149.160.56/RNentropy/	



with replicates:

# TR_COL	1
# GENE_COL	1
# EXP	cherry	2,3,4
# EXP	cress	5,6,7
# EXP	galium	8,9,10
# END

Mca00001	129.66	136.89	128.87	191.56	173.34	175.69	162.94	177.84	172.69
Mca00002	18.98	23.80	23.48	37.03	24.87	30.70	28.00	28.64	32.55
Mca00003	1.63	2.65	2.16	2.55	2.17	2.04	1.43	1.76	2.01
Mca00004	66.15	63.04	67.37	53.61	51.01	53.84	52.50	50.44	48.58
Mca00005	0.00	0.00	0.00	2.71	0.00	0.00	0.16	0.27	0.00
Mca00006	1.95	0.97	1.15	3.72	1.91	1.05	1.55	1.53	0.94



# TO RUN THE dE ANALYSIS: NEED TO FORMAT THE INFILE FIRST by putting the # INFO at the top of the matrix file. 
# infile has to be pre normalised. 
/storage/home/users/pjt6/shelf_apps/apps/RNentropy1.1.1/RNentropy -f Mc_cherry_cress_galium.matrix.TMM_normalized.FPKM



# parse results: 
# ./select_results Mc_cherry_cress_galium123.genes.counts.matrix.TMM_normalized.FPKM.summary.res GPV_threshold LPV_threshold sample_num rep_num

# BELOW WE HAVE 3 REPS FROM 3 SAMPLES (3, 3) ...... 0.001 is the post HOC multi gene correction threshold
/storage/home/users/pjt6/shelf_apps/apps/RNentropy1.1.1/select_results Mc_cherry_cress_galium.matrix.TMM_normalized.FPKM.summary.res 0.001 0.001 3 3
