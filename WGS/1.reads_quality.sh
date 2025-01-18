#!/bin/bash -l
#SBATCH -D /path/to/raw/fastq/files  
#SBATCH -o /path/to/raw/fastq/files/logs/qualityLog-%j.txt
#SBATCH -e /path/to/raw/fast1/files/logs/qualityLog-%j.err
#SBATCH -t 24:00:00
#SBATCH -J reads_quality
#SBATCH --nodes=1
#SBATCH --ntasks 8
#SBATCH --mem 24g
#SBATCH --mail-type=ALL # if you want emails, otherwise remove


module load fastqc

mkdir metrics

fastqc *.fastq.gz -o metrics

###after running go to metrics directory and apply multiqc . to summarize results
###it is neccesary to install multiqc
###for further information https://multiqc.info/docs/
