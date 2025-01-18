#!/bin/bash -l
#SBATCH -D /scratch/asingh3/Indian_Amaranth/WGS/
#SBATCH -o /scratch/asingh3/Indian_Amaranth/WGS/HaplotypeCaller-Log-%j.txt
#SBATCH -e /scratch/asingh3/Indian_Amaranth/WGS/HaplotypeCaller-Log-%j.err
#SBATCH -t 10:00:00
#SBATCH -J Indian_amaranth_Mapping_V3
#SBATCH --nodes=1
#SBATCH --ntasks 8
#SBATCH --mem 12g
#SBATCH --array=0-2
#SBATCH --mail-type=ALL # if you want emails, otherwise remove




module use /opt/rrzk/modules/experimental
module load bwamem2/2.2.1
module load samtools/1.13

#gunzip Reference/GCA_000001405.15_GRCh38_genomic.fna.gz
REFERENCE=Reference/Ref_genes.fasta
bwa-mem2 index $REFERENCE

PROVIDER=ILLUMINA

OUTPUTPATH=Bam_files_genes 
FASTQPATH=Input
mkdir -p $OUTPUTPATH
mkdir -p $OUTPUTPATH/metrics/


#for SAMPLES in /projects/ag-stetter/sequencing_data/indian_amaranth/set1/raw_data/*
IFS=$'\n' FILES=($(cat Samples.txt))

SAMPLES=${FILES[$SLURM_ARRAY_TASK_ID]}

#for ((i=0;i<${#SAMPLES[@]};++i));

#do 
	echo Maping reads of ${SAMPLES[i]}
SORTED_NAME=${OUTPUTPATH}/${SAMPLES[i]}.bam
echo $SORTED_NAME

bwa-mem2 mem -t 8 -R '@RG\tID:'${SAMPLES[i]}'\tSM:'${SAMPLES[i]}'\tCN:'${PROVIDER}'\tPL:illumina' $REFERENCE ${FASTQPATH}/${SAMPLES[i]}_R1_001.fastq.gz ${FASTQPATH}/${SAMPLES[i]}_R2_001.fastq.gz | samtools sort -O bam -o ${SORTED_NAME}


#echo mark duplicates
DEDUP_NAME=${OUTPUTPATH}/${SAMPLES[i]}.final.bam
METRICS_FILE=${OUTPUTPATH}/metrics/${SAMPLES[i]}.txt
java -Xmx4g -jar /home/asingh3/TOOLS/picard.jar MarkDuplicates INPUT=$SORTED_NAME OUTPUT=$DEDUP_NAME METRICS_FILE=$METRICS_FILE
samtools index $DEDUP_NAME

echo calculate samtools flagstat
samtools flagstat ${DEDUP_NAME} > ${OUTPUTPATH}/metrics/${SAMPLES[i]}.flagstat

echo removing sorted bam
rm $SORTED_NAME


