# Install packages
install.packages("readxl")

# Load required libraries
library(readxl)
library(ggplot2)
library(tidyr)
library(dplyr)

#Read in the Excel file
data <- read_excel("GTscore_Summary_comparison.xlsx", sheet = "locus_summary")
to_remove <- read_excel("GTscore_Summary_comparison.xlsx", sheet = "to_remove")

# Optional: check your column names
print(colnames(data))

#Reshape to long format
data_long <- data %>%
  pivot_longer(cols = c(primerprobeprop_EH, primerprobeprop_OG),
               names_to = "GenotypeType",
               values_to = "Proportion") %>%
  mutate(Proportion = as.numeric(Proportion))  # <-- Force numeric

# Filter data_long for the loci to be removed
highlight_loci <- data_long %>%
  filter(locus %in% to_remove$locus)
highlight_loci <- inner_join(data, to_remove, by = "locus")
highlight_long <- inner_join(highlight_loci, data_long, by = "locus")
#################################################################################
#Step 1: Let's create a histogram overlaying the proportion of primer probe reads captured by each locus for both methods (EH and OG)
histogram_loci <- ggplot(data_long, aes(x = Proportion, fill = GenotypeType)) +
  geom_histogram(position = "identity", bins = 30, alpha = 0.6, color = "black") +
  
  # Black base for outline
  geom_rug(data = highlight_long, aes(x = Proportion),
           sides = "b", size = 1.5, color = "black", inherit.aes = FALSE) +
  
  # Goldenrod overlay for visibility
  geom_rug(data = highlight_long, aes(x = Proportion),
           sides = "b", size = 0.8, color = "goldenrod", inherit.aes = FALSE) +
  
  scale_fill_manual(values = c("primerprobeprop_EH" = "red", "primerprobeprop_OG" = "blue"),
                    labels = c("EH", "OG")) +
  labs(
    title = "Histogram of Target Capture Rates (Proportions) across Loci",
    x = "Primer Probe Proportion",
    y = "Count (# of loci)",
    fill = "Method"
  ) +
  theme_minimal()#This figure shows us that EH methods performed poorly compared to OG methods

histogram_loci

ggsave(filename = "histogram_loci.jpeg", plot = histogram_loci, dpi = 300, width = 8, height = 6, units = "in")

#################################################################################
#Step 2: Let's plot a scatter plot to compare primer probe reads across loci, for both EH and OG methods
# Define the columns of interest
num_cols <- c("primerprobeprop_OG", "primerprobeprop_EH")

# Ensure numeric conversion of columns in 'data'
data <- data %>%
  mutate(across(all_of(num_cols), ~ suppressWarnings(as.numeric(trimws(.)))))

#Ensure numeric conversion of columns in 'highlight_loci'
highlight_loci <- highlight_loci %>%
  mutate(across(all_of(num_cols), ~ suppressWarnings(as.numeric(trimws(.)))))

#Filter out rows in 'data' with NA in any of the numeric columns
plot_df <- data %>%
  filter(if_all(all_of(num_cols), ~ !is.na(.)))

# Optional: check for bad rows (you can comment this out if not needed)
# bad_rows <- data %>%
#   filter(if_any(all_of(num_cols), is.na))
# View(bad_rows)

# Create scatter plot
scatter_loci <- ggplot(plot_df, aes(x = primerprobeprop_OG, y = primerprobeprop_EH)) +
  geom_point(color = "darkgreen", alpha = 0.7) +
  geom_point(
    data = highlight_loci,
    aes(x = primerprobeprop_OG, y = primerprobeprop_EH),
    color = "black", shape = 21, fill = "yellow", size = 3, stroke = 1.2
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  scale_x_continuous(
    name = "Original GTseq method (OG)",
    breaks = seq(0, 1, by = 0.2),
    limits = c(0, 1)
  ) +
  scale_y_continuous(
    name = "Emma's GTseq method (EH)",
    breaks = seq(0, 1, by = 0.2),
    limits = c(0, 1)
  ) +
  labs(title = "Primer Probe Proportion Across Loci: EH vs OG") +
  theme_minimal()

# Print and save the plot
scatter_loci
ggsave(filename = "scatter_loci.jpeg", plot = scatter_loci, dpi = 300, width = 8, height = 6, units = "in")

################################################################################
#Step 3: Let's create a histogram overlaying the proportion of primer probe reads captured by each individual for both methods (EH and OG) 
data2 <- read_excel("GTscore_Summary_comparison.xlsx", sheet = "individual_summary")

histogram_individuals <- ggplot(data2, aes(x = Primer.Probe.Proportion, fill = treatment)) +
  geom_histogram(position = "identity", bins = 30, alpha = 0.5, color = "black") +
  scale_fill_manual(values = c("OG" = "blue", "EH" = "red")) +
  labs(
    title = "Primer Probe Proportions for Individuals Under Two Treatments",
    x = "Primer Probe Proportion",
    y = "Number of Individuals"
  ) +
  theme_minimal() #This figure shows us that EH and OG methods performed somewhat 
#similarly at the low end, with OG methods producing a few high primer-probe results
ggsave(filename = "histogram_individuals.jpeg", plot = histogram_individuals, dpi = 300, width = 8, height = 6, units = "in")


#scatter plot for individual summary
# Step 1: Pivot wider to get one row per sample with separate EH and OG columns
data_wide <- data2 %>%
  pivot_wider(names_from = treatment, values_from = Primer.Probe.Proportion)

# Step 2: Ensure numeric (just in case)
data_wide$OG <- as.numeric(data_wide$OG)
data_wide$EH <- as.numeric(data_wide$EH)

###############################################################################
#Compare genotyping rate for each individual across both methods

indiv_genotype <- read_excel("concatenate_check.xlsx", sheet = "comparison")


# Step 2: Convert relevant columns to numeric
numeric_indiv_genotype <- indiv_genotype %>%
  mutate(
    genotypingrate_OG = as.numeric(genotypingrate_OG),
    genotypingrate_EH = as.numeric(genotypingrate_EH)
  )

# Step 3: Reshape data to long format for bar plot
data_long_numeric_indiv_genotype <- numeric_indiv_genotype %>%
  pivot_longer(
    cols = c(genotypingrate_OG, genotypingrate_EH),
    names_to = "Method",
    values_to = "GenotypingRate"
  ) %>%
  mutate(Method = recode(Method,
                         genotypingrate_OG = "OG",
                         genotypingrate_EH = "EH"))

# Step 4: Calculate mean genotyping rate for each method
mean_rates <- data_long_numeric_indiv_genotype %>%
  group_by(Method) %>%
  summarise(mean_rate = mean(GenotypingRate, na.rm = TRUE))

# Step 5: Create bar plot with horizontal dashed lines for means
indiv_genotype_plot <- ggplot(data_long_numeric_indiv_genotype,
                              aes(x = lab_sampleID, y = GenotypingRate, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(data = mean_rates,
             aes(yintercept = mean_rate, color = Method),
             linetype = "dashed", linewidth = 1) +
  scale_fill_manual(values = c("OG" = "blue", "EH" = "red")) +
  scale_color_manual(values = c("OG" = "blue", "EH" = "red"), guide = "none") +
  labs(
    title = "Genotyping Rate per Individual by Method",
    x = "Individual (lab_sampleID)",
    y = "Genotyping Rate",
    fill = "Method"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

# Step 6: Display the plot
indiv_genotype_plot
# Save the plot as a high-resolution JPEG
ggsave(
  filename = "genotyping_rate_per_individual.jpg",
  plot = indiv_genotype_plot,
  width = 12,         # Width in inches (adjust if needed)
  height = 6,         # Height in inches (adjust if needed)
  dpi = 300,
  units = "in",
  device = "jpeg"
)
