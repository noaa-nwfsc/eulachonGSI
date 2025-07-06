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

#Let's create a histogram overlaying the proportion of primer probe reads captured by each locus for both methods (EH and OG)
histogram_loci <-ggplot(data_long, aes(x = Proportion, fill = GenotypeType)) +
  geom_histogram(position = "identity", bins = 30, alpha = 0.6, color = "black") +
  geom_rug(data = highlight_long, aes(x = Proportion, color = GenotypeType),
           sides = "b", size = 0.7, inherit.aes = FALSE) +
  scale_fill_manual(values = c("primerprobeprop_EH" = "red", "primerprobeprop_OG" = "blue"),
                    labels = c("EH", "OG")) +
  scale_color_manual(values = c("primerprobeprop_EH" = "red", "primerprobeprop_OG" = "blue")) +
  labs(
    title = "Histogram of Target Capture Rates (Proportions) across Loci",
    x = "Primer Probe Proportion",
    y = "Count (# of loci)",
    fill = "Method",
    color = "Method"
  ) +
  theme_minimal()#This figure shows us that EH methods performed poorly compared to OG methods
ggsave(filename = "histogram_loci.jpeg", plot = histogram_loci, dpi = 300, width = 8, height = 6, units = "in")


#Let's plot a scatter plot to compare primer probe reads across loci, for both EH and OG methods
data$primerprobeprop_OG <- as.numeric(data$primerprobeprop_OG)
data$primerprobeprop_EH <- as.numeric(data$primerprobeprop_EH)

stopifnot(all(c("primerprobeprop_OG", "primerprobeprop_EH") %in% colnames(highlight_loci)))

scatter_loci <- ggplot(data, aes(x = primerprobeprop_OG, y = primerprobeprop_EH)) +
  geom_point(color = "darkgreen", alpha = 0.7) +
  # Overlay to_remove loci with different style
  geom_point(data = highlight_loci, aes(x = primerprobeprop_OG, y = primerprobeprop_EH),
             color = "black", shape = 21, fill = "yellow", size = 3, stroke = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  scale_x_continuous(
    name = "Original GTseq method (OG)",
    breaks = seq(0, 1, by = 0.2)
  ) +
  scale_y_continuous(
    name = "Emma's GTseq method (EH)",
    breaks = seq(0, 1, by = 0.2)
  ) +
  labs(title = "Primer Probe Proportion Across Loci: EH vs OG") +
  theme_minimal() #The plot shows most points falling below the slope line, which 
#indicates that the OG methods performed better for the majority of loci
ggsave(filename = "scatter_loci.jpeg", plot = scatter_loci, dpi = 300, width = 8, height = 6, units = "in")

#Now let's repeat these but for the individual summary instead
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
# Step 3: Create scatter plot
scatter_individuals <- ggplot(data2, aes(x = sample, y = treatme)) +
  geom_point(color = "darkgreen", alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  scale_x_continuous(
    name = "Original GTseq method (OG)",
    breaks = seq(0, 1, by = 0.2)
  ) +
  scale_y_continuous(
    name = "Emma's GTseq method (EH)",
    breaks = seq(0, 1, by = 0.2)
  ) +
  labs(title = "Primer Probe Proportion for Individuals: EH vs OG") +
  theme_minimal() #This plot shows, unlike the locus summary plot, a large 
#cluster of individuals that EH methods worked better for, and then a number of 
#samples spread over a large range that performed better with OG methods
ggsave(filename = "scatter_individuals", plot = scatter_individuals, dpi = 300, width = 8, height = 6, units = "in")


#histogram for individual summary of off-target reads
ggplot(data_long2, aes(x = Proportion, fill = GenotypeType)) +
  geom_histogram(position = "identity", bins = 30, alpha = 0.6, color = "black") +
  scale_fill_manual(values = c("Primer.Probe.ProportionEH" = "red", "Primer.Probe.ProportionOG" = "blue"),
                    labels = c("EH", "OG")) +
  labs(
    title = "Histogram of Target Capture Rates (Proportions) for Individuals",
    x = "Primer Probe Proportion",
    y = "Count (# of loci)",
    fill = "Method"
  ) +
  theme_minimal()

