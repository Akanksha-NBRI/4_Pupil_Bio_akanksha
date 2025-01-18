#!/bin/bash -l
#SBATCH -D /scratch/WGS/
#SBATCH -o /scratch/WGS/SNPCall-Log-%j.txt
#SBATCH -e /scratch/WGS/SNPCall-Log-%j.err
#SBATCH -t 12:00:00
#SBATCH -J Variant-Call
#SBATCH --partition=smp
#SBATCH --mem 8g
#SBATCH --mail-type=ALL # if you want emails, otherwise remove

module load samtools/1.13
# index file 
samtools faidx Reference/Ref_genes.fasta

java -Xmx4g -jar /home/asingh3/TOOLS/picard.jar CreateSequenceDictionary REFERENCE=Reference/Ref_genes.fasta OUTPUT=Reference/Ref_genes.dict 

#Variant calling
module load miniconda
conda activate GATK

gatk Mutect2 -R Reference/Ref_genes.fasta -I Bam_files_genes/PA220KH-lib09-P19-Tumor_S2_L001.final.bam -I Bam_files_genes/PA221MH-lib09-P19-Norm_S1_L001.final.bam -normal PA221MH-lib09-P19-Norm_S1_L001 -O Pupil_mutect2_genes.vcf.gz

gatk FilterMutectCalls -R Reference/Ref_genes.fasta -V Pupil_mutect2_genes.vcf.gz -O Pupil_mutect2_genes.filtered.vcf.gz

#Estimate statsistic
bcftools stats -s - Pupil_mutect2_genes.vcf.gz -O Pupil_mutect2_genes.filtered.vcf.gz > Bcftools.stats

#mutation rate estimation
python MutationRate.py Pupil_mutect2_genes.filtered.vcf.gz > Mutation.rate
