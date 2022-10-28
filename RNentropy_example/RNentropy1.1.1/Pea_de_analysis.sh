cd /home/pt40963/scratch/Downloads/RNentropy1.1.1

# TO RUN THE dE ANALYSIS: NEED TO FORMAT THE INFILE FIRST
./RNentropy -f Published_pea_N116_12345_PS01_12345_avr_12345_vir12345.published_models.genes.counts.matrix.TMM_normalized.FPKM

 ./RNentropy -f Api_N116_12345_PS01_12345_avr_12345_vir12345.genes.counts.matrix.TMM_normalized.FPKM

 
./RNentropy -f Mc_cherry_cress_galium123.genes.counts.matrix.TMM_normalized.FPKM

 
# parse results: 
# ./select_results Mc_cherry_cress_galium123.genes.counts.matrix.TMM_normalized.FPKM.summary.res GPV_threshold LPV_threshold sample_num rep_num

# BELOW WE HAVE 3 REPS FROM 3 SAMPLES ...... 
./select_results Published_pea_N116_12345_PS01_12345_avr_12345_vir12345.published_models.genes.counts.matrix.TMM_normalized.FPKM.summary.res 0.001 0.001 4 5


./select_results  Api_N116_12345_PS01_12345_avr_12345_vir12345.genes.counts.matrix.TMM_normalized.FPKM.summary.res 0.001 0.001 4 5




./select_results Mc_cherry_cress_galium123.genes.counts.matrix.TMM_normalized.FPKM.summary.res 0.001 0.001 3 3


./select_results Published_pea_N116_12345_PS01_12345_avr_12345_vir12345.published_models.genes.counts.matrix.TMM_normalized.FPKM.summary.res 0.01 0.01 4 5


./select_results  Api_N116_12345_PS01_12345_avr_12345_vir12345.genes.counts.matrix.TMM_normalized.FPKM.summary.res 0.01 0.01 4 5
