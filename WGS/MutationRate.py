import pysam

# Input files
filtered_vcf = "/home/alle/Desktop/Personal/Pupib_Bio/WGS/WithNEW-REF/WGS/Pupil_mutect2_genes.filtered.vcf.gz"
effective_genome_size_mb = 141  # Example: 30 Mb for exonic regions in WES

# Count somatic variants
vcf = pysam.VariantFile(filtered_vcf)
somatic_variant_count = sum(1 for rec in vcf)

# Calculate mutation rate
mutation_rate = somatic_variant_count / effective_genome_size_mb
print(f"Total Somatic Variants: {somatic_variant_count}")
print(f"Mutation Rate: {mutation_rate:.2f} mutations/Mb")

