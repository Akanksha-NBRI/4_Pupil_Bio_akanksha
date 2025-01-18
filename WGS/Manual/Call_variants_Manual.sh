#!/bin/bash

# Input files and reference genome
NORMAL_BAM="WGS/Bam_files_genes/PA221MH-lib09-P19-Norm_S1_L001.final.bam"
TUMOR_BAM="WGS/Bam_files_genes/PA220KH-lib09-P19-Tumor_S2_L001.final.bam"
REFERENCE_GENOME="WGS/WithNEW-REF/WGS/Reference/Ref_genes.fasta"

# Output files
NORMAL_VCF="normal_variants.vcf"
TUMOR_VCF="tumor_variants.vcf"
DIFFERENTIAL_VCF="differential_variants.vcf"
METRICS_R_SCRIPT="metrics_calculation.R"


#Call variants for normal and tumor samples
echo "Calling variants for normal sample..."
bcftools mpileup -Ou -f $REFERENCE_GENOME $NORMAL_BAM | bcftools call -mv -Ov -o $NORMAL_VCF
bgzip $NORMAL_VCF
tabix -p vcf $NORMAL_VCF.gz
echo "Calling variants for tumor sample..."
bcftools mpileup -Ou -f $REFERENCE_GENOME $TUMOR_BAM | bcftools call -mv -Ov -o $TUMOR_VCF
bgzip $TUMOR_VCF
tabix -p vcf $TUMOR_VCF.gz
#Identify differential variants between tumor and normal
echo "Identifying differential variants..."
bcftools isec -p output_diff -Ov $NORMAL_VCF.gz $TUMOR_VCF.gz
mv output_diff/0003.vcf $DIFFERENTIAL_VCF
#Filter the file for quality
bcftools filter -e "QUAL<30 || DP<10" -Ov -o filtered_$DIFFERENTIAL_VCF $DIFFERENTIAL_VCF
bcftools filter -e "QUAL<30 || DP<10" -Ov -o filtered_$NORMAL_VCF $NORMAL_VCF.gz
bcftools filter -e "QUAL<30 || DP<10" -Ov -o filtered_$TUMOR_VCF $TUMOR_VCF.gz

## EStimate statistics
bcftools stats -s - filtered_$DIFFERENTIAL_VCF $DIFFERENTIAL_VCF > Differential.stats
bcftools stats -s -filtered_$NORMAL_VCF $NORMAL_VCF.gz > Normal.stats
bcftools stats -s -filtered_$TUMOR_VCF $TUMOR_VCF.gz > Tumor.stats



