#!/bin/bash
#$ -V ## pass all environment variables to the job, VERY IMPORTANT
#$ -N Trimmo ## job name
#$ -S /bin/bash ## shell where it will run this job
#$ -j y ## join error output to normal output
#$ -cwd ## Execute the job from the current working directory
#$ -pe multi 2 # I am a broken line!!!


cd $HOME/genome_assembly_workshop 

###########################################################################
# this is how we ran fastqc (run this again on the trimmed reads)
#fastqc ecoli_R1.fastq.gz


###########################################################################
# this is how we will trim the reads - trimmomatic
#java -jar trimmomatic-0.35.jar PE -phred33 input_forward.fq.gz input_reverse.fq.gz output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
#
#This will perform the following:
#
#    Remove adapters (ILLUMINACLIP:TruSeq3-PE.fa:2:30:10)
#    Remove leading low quality or N bases (below quality 3) (LEADING:3)
#    Remove trailing low quality or N bases (below quality 3) (TRAILING:3)
#    Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15 (SLIDINGWINDOW:4:15)
#    Drop reads below the 36 bases long (MINLEN:36)
# \ chacter means you can look at one long command.  over multiple lines. 

java -jar /shelf/training/Trimmomatic-0.38/trimmomatic-0.38.jar PE -summary trim_summary.txt \
-threads 2 -phred33 ./reads/subsampled_R1.fastq.gz ./reads/subsampled_R2.fastq.gz subsampled_R1_paired.fastq.gz \
subsampled_R1_unpaired.fastq.gz subsampled_R2_paired.fastq.gz subsampled_R2_unpaired.fastq.gz \
ILLUMINACLIP:/shelf/training/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:45 
