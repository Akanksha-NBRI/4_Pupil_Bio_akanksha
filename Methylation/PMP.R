# libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# load dataset
data <- read.table("PupilBioTest_PMP_revA.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

# total PMP occurrences across tissues
data <- data %>%
  mutate(Total_Count = rowSums(select(., `X.000`, `X.001`, `X.010`, `X.011`, `X.100`, `X.101`, `X.110`, `X.111`)))

# Identify PMPs for each tissue
# Aggregate PMP counts per tissue
pmp_tissue_summary <- data %>%
  group_by(CpG_Coordinates, Tissue) %>%
  summarise(
    Total_Count = sum(Total_Count, na.rm = TRUE),
    .groups = "drop"
  )

# Normalize tissues
pmp_specificity <- pmp_tissue_summary %>%
  group_by(CpG_Coordinates) %>%
  mutate(
    Total_Count_Others = sum(Total_Count) - Total_Count,  
    Specificity_Score = Total_Count / (Total_Count + Total_Count_Others),  # Specificity formula
    Tissue_Rank = if_else(Tissue == "cfDNA", 1, 0)  
  )

write.csv(pmp_specificity, file="PMP_Specificity_RAW.csv", row.names=FALSE)

# high-specificity PMPs
high_specificity_pmps <- pmp_specificity %>%
  filter(Tissue_Rank == 1) %>%
  arrange(desc(Specificity_Score)) %>%
  filter(Specificity_Score > 0.9)  # Threshold for high specificity

write.csv(high_specificity_pmps, "high_specificity_pmps.csv", row.names = FALSE)


#Plotting
histogram <- ggplot(high_specificity_pmps, aes(x = Specificity_Score)) +
  geom_histogram(binwidth = 0.01, fill = "blue", color = "black", alpha = 0.7) +
  labs(
    title = "Distribution of PMP Specificity Scores",
    x = "Specificity Score",
    y = "Frequency"
  ) +
  theme_minimal()

top_pmps <- high_specificity_pmps %>%
  filter(Tissue == "cfDNA") %>%
  arrange(desc(Specificity_Score)) %>%
  head(10)  # Top 10 PMPs

barplot <- ggplot(top_pmps, aes(x = reorder(CpG_Coordinates, Specificity_Score), y = Specificity_Score, fill = Tissue)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    title = "Top 10 PMPs with Highest Specificity for cfDNA",
    x = "CpG Coordinates",
    y = "Specificity Score"
  ) +
  theme_minimal()

threshold_data <- high_specificity_pmps %>%
  arrange(desc(Specificity_Score)) %>%
  mutate(Cumulative_Count = row_number())  # Count PMPs as threshold decreases

cumulative_plot <- ggplot(threshold_data, aes(x = Specificity_Score, y = Cumulative_Count)) +
  geom_line(color = "darkred", size = 1) +
  labs(
   title = "Threshold Optimization: PMP Count vs Specificity",
    x = "Specificity Threshold",
    y = "Number of PMPs Selected"
  ) +
  theme_minimal()

#Save Plots
ggsave("specificity_histogram.png", plot = histogram, width = 8, height = 6)
ggsave("top_pmps_barplot.png", plot = barplot, width = 8, height = 6)
ggsave("threshold_optimization_plot.png", plot = cumulative_plot, width = 8, height = 6)


#P-value calculation
specificity_data <- read.csv("high_specificity_pmps.csv")

# contingency table and compute p-values for each PMP
pmp_p_values <- specificity_data %>%
  rowwise() %>%
  mutate(
    Counts_Tissue1 = Total_Count,  
    Counts_Other = Total_Count_Others,  
 #    Fisher's Exact Test
    P_Value = fisher.test(matrix(c(Counts_Tissue1, Counts_Other, sum(Counts_Tissue1), sum(Counts_Other)),
                                 nrow = 2))$p.value
  )

# Step 4: Adjust p-values using Benjamini-Hochberg (FDR)
pmp_p_values <- pmp_p_values %>%
  mutate(Adjusted_P_Value = p.adjust(P_Value, method = "BH"))

# Step 5: Save the results
write.csv(pmp_p_values, "pmp_with_pvalues.csv", row.names = FALSE)

# Display results
print(pmp_p_values)


###################################################################

#pvalue using permutation test
set.seed(123)  # Ensure reproducibility

# Define permutation function
permutation_test <- function(obs_count, other_counts, n_permutations = 100) {
  total = obs_count + other_counts
  null_dist <- replicate(n_permutations, sum(sample(total, obs_count)))
  mean(null_dist >= obs_count) 
}

# Apply permutation test to each PMP
pmp_permutation <- specificity_data %>%
  rowwise() %>%
  mutate(
    P_Value = permutation_test(Total_Count, Total_Count_Others, n_permutations = 100)
  )

# Adjust p-values (Benjamini-Hochberg)
pmp_permutation <- pmp_permutation %>%
  mutate(Adjusted_P_Value = p.adjust(P_Value, method = "BH"))

# Save results
write.csv(pmp_permutation, "pmp_with_permutation_pvalues.csv", row.names = FALSE)

significant_pmps <- pmp_permutation %>%
  filter(Adjusted_P_Value < 0.05)
write.csv(significant_pmps, "significant_pmps.csv", row.names = FALSE)

#visualize volcano plot

volcano_plot <- ggplot(pmp_permutation, aes(x = Specificity_Score, y = -log10(Adjusted_P_Value))) +
  geom_point(aes(color = Adjusted_P_Value < alpha), alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  labs(
    title = "Volcano Plot of PMP Specificity and Significance",
    x = "Specificity Score",
    y = "-log10(Adjusted P-Value)",
    color = "Significant"
  ) +
  theme_minimal()

# Save and display the plot
ggsave("volcano_plot_pmps.png", plot = volcano_plot, width = 8, height = 6)
print(volcano_plot)


