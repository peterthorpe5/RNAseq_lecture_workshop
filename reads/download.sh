#!/bin/bash
#$ -V ## pass all environment variables to the job, VERY IMPORTANT
#$ -N FastQCTraining ## job name
#$ -S /bin/bash ## shell where it will run this job
#$ -j y ## join error output to normal output
#$ -cwd ## Execute the job from the current working directory
#$ -l hostname=marvin  # can only download from this node

# this line takes us into the correct directory where the data is
cd $HOME/genome_assembly_workshop/

# dowload all from here ...https://www.ebi.ac.uk/ena/browser/view/PRJEB24338

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/003/ERR2239823/ERR2239823_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/003/ERR2239823/ERR2239823_2.fastq.gz

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/004/ERR2239824/ERR2239824_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/004/ERR2239824/ERR2239824_2.fastq.gz