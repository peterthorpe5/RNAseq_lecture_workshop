cd $HOME
STAR --genomeDir star_indicies/ \
--limitBAMsortRAM 455554136874 --runThreadN 16 \
--readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 7  \
--outFilterMultimapNmax 5 --outFileNamePrefix DARKTELLEFT_1 \
--readFilesIn WTCHG_328669_209_paired_R1.fq.gz WTCHG_328669_209_paired_R2.fq.gz 
