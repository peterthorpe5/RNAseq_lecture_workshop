
conda activate python27

htseq-count  -s no  -r pos -f bam --type gene --idattr ID ppenAligned.sortedByCoord.out.bam genome.gff3 > condition_rep_RNAseq.counts



