#!/bin/bash
#$ -V ## pass all environment variables to the job, VERY IMPORTANT
#$ -N FastQCTraining ## job name
#$ -S /bin/bash ## shell where it will run this job
#$ -j y ## join error output to normal output
#$ -cwd ## Execute the job from the current working directory

# this line takes us into the correct directory where the data is
cd $HOME/genome_assembly_workshop/

# this line loads fastqc
echo "loading Fastqc module"
module load FASTQC


# this line run fastqc on a data file.
# * is a wild card to match anything, for us R1 and R2. 
fastqc $HOME/genome_assembly_workshop/reads/subsampled_R*.fastq.gz
echo "finished"



