# Step 1: Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# Step 2: Load the dataset
data <- read.table("PupilBioTest_PMP_revA.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

# Step 3: Calculate total PMP occurrences across tissues
data <- data %>%
  mutate(Total_Count = rowSums(select(., `X.000`, `X.001`, `X.010`, `X.011`, `X.100`, `X.101`, `X.110`, `X.111`)))

pmp_data <- data

# Compute Total Reads and VRF for each PMP
pmp_data <- pmp_data %>%
  rowwise() %>%
  mutate(
    Total_Reads = sum(c_across(`X.000`:`X.111`)),  
    Variant_Reads = Total_Reads - `X.000`,       
    VRF = Variant_Reads / Total_Reads )      


# Group by Tissue and calculate mean VRF for each PMP
mean_vrf_by_tissue <- pmp_data %>%
  group_by(CpG_Coordinates, Tissue) %>%
  summarize(
    Mean_VRF = mean(VRF, na.rm = TRUE),      
    Median_VRF = median(VRF, na.rm = TRUE),  
    SD_VRF = sd(VRF, na.rm = TRUE))

write.csv(mean_vrf_by_tissue, "mean_vrf_by_tissue.csv", row.names = FALSE)


# Create a boxplot of VRF by Tissue
vrf_plot <- ggplot(mean_vrf_by_tissue, aes(x = Tissue, y = Mean_VRF, fill = Tissue)) +
  geom_boxplot(alpha = 0.7) +
  labs(
    title = "Variant Read Fraction (VRF) by Tissue",
    x = "Tissue",
    y = "Mean VRF") +
  theme_minimal()

ggsave("vrf_boxplot.png", plot = vrf_plot, width = 8, height = 6)

###############Confidence with  simulation
set.seed(123)
depth_data <- data.frame(
  Sequencing_Depth = rep(seq(10, 100, by=10), each=10),
  Specificity_Confidence = rep(seq(0.1, 1, by=0.1), each=10) + rnorm(100, 0, 0.05)
)

# Plot
ggplot(depth_data, aes(x=Sequencing_Depth, y=Specificity_Confidence)) +
  geom_point() +
  geom_smooth(method="loess", color="blue") +
  labs(title="Impact of Sequencing Depth on Specificity Confidence",
    x="Sequencing Depth (Total Reads)",
    y="Specificity Confidence") +
  theme_minimal()
ggsave("plot.png")

#### 1 million read#####################################

# Filter for Tissue #2 and select top 10 PMPs with highest specificity
top_10_pmps <- pmp_data %>%
 filter(Tissue == "Islet") %>%
 arrange(desc(VRF)) %>%
  slice(1:10)  # Select top 10

# Define total sequencing depth
total_reads <- 1e6  # 1 million reads

#Calculate required reads for confident Tissue #2 identification
top_10_pmps <- top_10_pmps %>%
  mutate(Required_Reads = total_reads * VRF,  
   CI_Lower = Required_Reads - 1.96 * sqrt(Required_Reads), 
   CI_Upper = Required_Reads + 1.96 * sqrt(Required_Reads))

# Save the read thresholds 
write.csv(top_10_pmps, "top_10_pmps_read_thresholds.csv", row.names = FALSE)
# plotting
ggplot(top_10_pmps, aes(x = reorder(CpG_Coordinates, -Required_Reads), y = Required_Reads)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.4, color = "black") +
  labs(title = "Required Reads for Confidently Calling Tissue #2",
    x = "PMP (CpG Coordinates)",
    y = "Required Reads") +
  theme_minimal() +
  coord_flip()



############Validation
# Split PMP coordinates into individual CpG sites
pmp_data <- pmp_data %>%
  mutate(CpG_Sites = strsplit(as.character(CpG_Coordinates), ":")) %>%
  unnest_longer(CpG_Sites)

# Calculate specificity (VRF) for individual CpG sites
individual_cpgs <- pmp_data %>%
  group_by(CpG_Sites, Tissue) %>%
  summarize(
    VRF_CpG = mean(VRF, na.rm = TRUE),  # Mean VRF for individual CpG
    .groups = "drop"
  )
write.csv(individual_cpgs, file="Individual_CpGs_out.csv")

# Summarize PMP specificity
pmp_specificity <- pmp_data %>%
  group_by(CpG_Coordinates, Tissue) %>%
  summarize(
    VRF_PMP = mean(VRF, na.rm = TRUE),  # Mean VRF for PMP
    .groups = "drop"
  )
write.csv(pmp_specificity, file="VRF_PMP_Spec-out-validatayion.csv")

# Merge PMP and individual CpG data
comparison_data <- inner_join(
  individual_cpgs,
  pmp_specificity,
  by = c("CpG_Sites" = "CpG_Coordinates", "Tissue")
)

write.csv(comparison_data, file="Comparison_Data_out.csv")

# Plotggplot(data, aes(x = VRF_CpG, y = Specificity_Score, color = Tissue)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(title = "Specificity Score vs VRF for CpG Sites",
    x = "VRF (Variant Read Fraction)",
    y = "Specificity Score") +
  theme_minimal() +
  scale_color_manual(values = c("Islet" = "blue"))  
ggsave("Validation.png")
