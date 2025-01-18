library(VariantAnnotation)

# Load the VCF file
vcf <- readVcf("Normal.Variants.recode.vcf", genome="Reference/Ref_genes.fasta")

info_data <- info(vcf)  # Extract INFO fields
geno_data <- geno(vcf)  # Extract FORMAT fields
# Extract strand bias information
strand_bias <- info_data$SB
strand_bias[!is.na(strand_bias)]

# Extract allele depth (AD) and calculate the ratio of alternate allele to total reads
allele_depth <- geno_data$AD
alt_to_total_ratio <- sapply(allele_depth, function(x) x[2] / sum(x))
alt_to_total_ratio <- na.omit(alt_to_total_ratio)  # Remove NA values

# Summarize the ratio
summary(alt_to_total_ratio)
# Extract Read Position Rank Sum
read_position_bias <- info_data$ReadPosRankSum
summary(read_position_bias[!is.na(read_position_bias)])

# Extract Quality by Depth (QD)
quality_by_depth <- info_data$QD
summary(quality_by_depth[!is.na(quality_by_depth)])
# Extract depth of coverage
depth_of_coverage <- geno_data$DP
summary(depth_of_coverage[!is.na(depth_of_coverage)])
#Plot seq errors
hist(alt_to_total_ratio, breaks = 50, main = "Histogram of Alternate Allele Ratio",
     xlab = "Alternate Allele Depth / Total Depth", col = "blue")
dev.copy(png,"AA_Depth.png")
dev.off()
median(alt_to_total_ratio)
background <- median(alt_to_total_ratio)


# Confidence threshold: Reads required per million to detect a mutation
calculate_rpm <- function(background, confidence=0.99) {
  # Poisson approximation: RPM required to exceed background
  # log(1-confidence) ~ Poisson cumulative prob at given rate
  # lambda = RPM * background
  reads_required <- qpois(confidence, lambda=1 / background)
  rpm_required <- reads_required / 1e6
  return(rpm_required)
}


# Calculate RPM for various confidence levels
confidence_levels <- c(0.95, 0.99, 0.999)
rpm_values <- sapply(confidence_levels, calculate_rpm, background=median_background)

# Output results
cat("Reads per million (RPM) required for mutation detection at various confidence levels:\n")
for (i in seq_along(confidence_levels)) {
  cat(sprintf("Confidence %.1f%%: %.2f RPM\n", confidence_levels[i] * 100, rpm_values[i]))
}

# Save results to a file
writeLines(c(
  paste("Median Background Mutation Level (VAF):", median_background),
  paste("Confidence Levels:", paste(confidence_levels, collapse=", ")),
  paste("RPM Values:", paste(round(rpm_values, 2), collapse=", "))
), con="background_mutation_results.txt")

