#Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)

#Load the dataset
data <- read.table("PupilBioTest_PMP_revA.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

# Calculate total reads per row
data <- data %>%
  mutate(Total_Count = rowSums(select(., `X.000`, `X.001`, `X.010`, `X.011`, `X.100`, `X.101`, `X.110`, `X.111`)))

# Split CpG_Coordinates into individual positions
expanded_data <- data %>%
  separate_rows(CpG_Coordinates, sep = ":") %>%
  mutate(Coverage_Per_CpG = Total_Count / 3)  # Divide total count across 3 CpG positions

# Aggregate coverage for each CpG position and tissue
cpg_coverage_tissue <- expanded_data %>%
  group_by(Tissue, CpG_Coordinates) %>%
  summarise(Total_Coverage = sum(Coverage_Per_CpG, na.rm = TRUE), .groups = 'drop') %>%
  arrange(Tissue, as.numeric(CpG_Coordinates))

write.csv(cpg_coverage_tissue, "cpg_coverage_by_tissue.csv", row.names = FALSE)

#Calculate summary statistics
summary_stats <- cpg_coverage_tissue %>%
  group_by(Tissue) %>%
  summarise(
    Median_Coverage = median(Total_Coverage, na.rm = TRUE),
    Mean_Coverage = mean(Total_Coverage, na.rm = TRUE),
    SD_Coverage = sd(Total_Coverage, na.rm = TRUE),
    CV_Coverage = SD_Coverage / Mean_Coverage * 100
  )

write.csv(summary_stats, "median_and_cv_by_tissue.csv", row.names = FALSE)

# Visualizations
boxplot <- ggplot(cpg_coverage_tissue, aes(x = Tissue, y = Total_Coverage, fill = Tissue)) +
  geom_boxplot() +
  labs(title = "Coverage Distribution by Tissue", x = "Tissue", y = "CpG Coverage") +
  theme_minimal()

barplot <- ggplot(summary_stats, aes(x = Tissue, y = Median_Coverage, fill = Tissue)) +
  geom_bar(stat = "identity") +
  labs(title = "Median Coverage by Tissue", x = "Tissue", y = "Median Coverage") +
  theme_minimal()

scatterplot <- ggplot(summary_stats, aes(x = Tissue, y = Mean_Coverage)) +
  geom_point(size = 4, color = "blue") +
  geom_errorbar(aes(ymin = Mean_Coverage - SD_Coverage, ymax = Mean_Coverage + SD_Coverage), width = 0.2) +
  geom_text(aes(label = paste0("CV: ", round(CV_Coverage, 1), "%")), vjust = -1.5, color = "darkred") +
  labs(title = "Mean Coverage and CV by Tissue", x = "Tissue", y = "Mean Coverage") +
  theme_minimal()

ggsave("boxplot_coverage_distribution.png", plot = boxplot, width = 8, height = 6)
ggsave("barplot_median_coverage.png", plot = barplot, width = 8, height = 6)
ggsave("scatterplot_mean_cv_coverage.png", plot = scatterplot, width = 8, height = 6)

