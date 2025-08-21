library(readr)
library(dplyr)

# Read file (tab-delimited assumed)
df <- read_delim("singleSNP_summary.txt", delim = "\t", show_col_types = FALSE)

# --- Make numeric: handles "1,234", "87%", etc. ---
df <- df %>%
  mutate(
    AvgReadDepth_raw   = as.character(AvgReadDepth),
    GenotypeRate_raw   = as.character(GenotypeRate),
    AvgReadDepth       = readr::parse_number(AvgReadDepth_raw),
    GenotypeRate       = readr::parse_number(GenotypeRate_raw)
  )

# If GenotypeRate looks like percentages (0–100), convert to proportion (0–1)
geno_max <- suppressWarnings(max(df$GenotypeRate, na.rm = TRUE))
if (is.finite(geno_max) && geno_max > 1 && geno_max <= 100) {
  df <- df %>% mutate(GenotypeRate = GenotypeRate / 100)
}

# --- Quartiles ---
depth_q <- quantile(df$AvgReadDepth, probs = c(0.25, 0.75), na.rm = TRUE)
geno_q  <- quantile(df$GenotypeRate, probs = c(0.25, 0.75), na.rm = TRUE)

# Top/bottom lists (by metric)
depth_bottom <- df %>% filter(AvgReadDepth <= depth_q[1]) %>% pull(Locus_ID)
depth_top    <- df %>% filter(AvgReadDepth >= depth_q[2]) %>% pull(Locus_ID)

geno_bottom  <- df %>% filter(GenotypeRate <= geno_q[1])  %>% pull(Locus_ID)
geno_top     <- df %>% filter(GenotypeRate >= geno_q[2])  %>% pull(Locus_ID)

# Overlaps
overlap_top    <- intersect(depth_top, geno_top)
overlap_bottom <- intersect(depth_bottom, geno_bottom)

# Save lists
write_lines(depth_top,        "depth_top.txt")
write_lines(depth_bottom,     "depth_bottom.txt")
write_lines(geno_top,         "geno_top.txt")
write_lines(geno_bottom,      "geno_bottom.txt")
write_lines(overlap_top,      "overlap_top.txt")
write_lines(overlap_bottom,   "overlap_bottom.txt")

# Annotated master table
combined <- df %>%
  mutate(
    Depth_Group = case_when(
      Locus_ID %in% depth_top    ~ "top_depth",
      Locus_ID %in% depth_bottom ~ "bottom_depth",
      TRUE ~ "other"
    ),
    Genotype_Group = case_when(
      Locus_ID %in% geno_top     ~ "top_geno",
      Locus_ID %in% geno_bottom  ~ "bottom_geno",
      TRUE ~ "other"
    )
  )
write_delim(combined, "singleSNP_with_quartiles.txt", delim = "\t")

# Counts
summary_counts <- tibble(
  group = c("depth_top","depth_bottom","geno_top","geno_bottom","overlap_top","overlap_bottom"),
  count = c(length(depth_top), length(depth_bottom), length(geno_top), length(geno_bottom),
            length(overlap_top), length(overlap_bottom))
)
print(summary_counts)
write_delim(summary_counts, "quartile_counts.txt", delim = "\t")

# ------------------------------
# Min and Max for each quartile
# ------------------------------

quartile_ranges <- tibble(
  Metric = c("AvgReadDepth","AvgReadDepth","GenotypeRate","GenotypeRate"),
  Quartile = c("top","bottom","top","bottom"),
  Min = c(
    min(df$AvgReadDepth[df$AvgReadDepth >= depth_q[2]], na.rm = TRUE),
    min(df$AvgReadDepth[df$AvgReadDepth <= depth_q[1]], na.rm = TRUE),
    min(df$GenotypeRate[df$GenotypeRate >= geno_q[2]], na.rm = TRUE),
    min(df$GenotypeRate[df$GenotypeRate <= geno_q[1]], na.rm = TRUE)
  ),
  Max = c(
    max(df$AvgReadDepth[df$AvgReadDepth >= depth_q[2]], na.rm = TRUE),
    max(df$AvgReadDepth[df$AvgReadDepth <= depth_q[1]], na.rm = TRUE),
    max(df$GenotypeRate[df$GenotypeRate >= geno_q[2]], na.rm = TRUE),
    max(df$GenotypeRate[df$GenotypeRate <= geno_q[1]], na.rm = TRUE)
  )
)

print(quartile_ranges)
write_delim(quartile_ranges, "quartile_min_max.txt", delim = "\t")

library(ggplot2)
# Quartile values
depth_q <- quantile(combined$AvgReadDepth, probs = c(0.25, 0.75), na.rm = TRUE)
bottom_q <- depth_q[1]
top_q    <- depth_q[2]

# Determine x-axis limits (focus on main range)
x_min <- min(combined$AvgReadDepth, na.rm = TRUE)
x_max <- quantile(combined$AvgReadDepth, 0.99, na.rm = TRUE)  # exclude extreme outliers

# Histogram with shaded quartiles and fewer x-axis labels
depth_plot <- ggplot(combined, aes(x = AvgReadDepth)) +
  # Histogram in light gray
  geom_histogram(binwidth = 20, fill = "lightgray", color = "black") +
  
  # Shaded bottom quartile
  annotate("rect", xmin = -Inf, xmax = bottom_q, ymin = 0, ymax = Inf, 
           alpha = 0.2, fill = "orange") +
  
  # Shaded top quartile
  annotate("rect", xmin = top_q, xmax = Inf, ymin = 0, ymax = Inf, 
           alpha = 0.2, fill = "skyblue") +
  
  # Vertical dashed lines at quartile cutoffs
  geom_vline(xintercept = bottom_q, linetype = "dashed", color = "orange", size = 1) +
  geom_vline(xintercept = top_q, linetype = "dashed", color = "skyblue", size = 1) +
  
  # Focus x-axis on main range and reduce number of labels
  scale_x_continuous(limits = c(x_min, x_max), 
                     breaks = seq(0, x_max, by = 200)) +  # fewer labels
  
  theme_minimal() +
  labs(
    title = "Histogram of AvgReadDepth with Top/Bottom Quartiles",
    x = "AvgReadDepth",
    y = "Count"
  )
ggsave("AvgReadDepth_hist.jpeg", plot = depth_plot, width = 8, height = 5, dpi = 300)

# Histogram: GenotypeRate
geno_plot <- ggplot(combined, aes(x = GenotypeRate, fill = Genotype_Group)) +
  geom_histogram(binwidth = 0.05, alpha = 0.6, position = "identity", color = "black") +
  scale_fill_manual(values = c("top_geno" = "skyblue", 
                               "bottom_geno" = "orange", 
                               "other" = "grey")) +
  theme_minimal() +
  labs(title = "Histogram of GenotypeRate by Quartile",
       x = "GenotypeRate",
       y = "Count",
       fill = "Quartile")
ggsave("GenotypeRate_hist.jpeg", plot = geno_plot, width = 8, height = 5, dpi = 300)

###### now,  identify the loci as high or low quartile loci ######
library(dplyr)
library(readr)
library(stringr)

# ------------------------------
# Read data
# ------------------------------
primers <- read_csv("primers_cleaned.csv")
snp_summary <- read_delim("singleSNP_summary.txt", delim = "\t", show_col_types = FALSE)

# Make numeric
snp_summary <- snp_summary %>%
  mutate(
    AvgReadDepth = readr::parse_number(AvgReadDepth),
    GenotypeRate = readr::parse_number(GenotypeRate)
  )

# If GenotypeRate > 1, convert to proportion
if(max(snp_summary$GenotypeRate, na.rm = TRUE) > 1){
  snp_summary <- snp_summary %>%
    mutate(GenotypeRate = GenotypeRate / 100)
}

# ------------------------------
# Calculate quartiles
# ------------------------------
depth_q <- quantile(snp_summary$AvgReadDepth, probs = c(0.25, 0.75), na.rm = TRUE)
geno_q  <- quantile(snp_summary$GenotypeRate, probs = c(0.25, 0.75), na.rm = TRUE)

# ------------------------------
# Assign quartile membership per metric
# ------------------------------
snp_summary <- snp_summary %>%
  mutate(
    Depth_Quartile = case_when(
      AvgReadDepth <= depth_q[1] ~ "low",
      AvgReadDepth >= depth_q[2] ~ "high",
      TRUE                        ~ "middle"
    ),
    Genotype_Quartile = case_when(
      GenotypeRate <= geno_q[1] ~ "low",
      GenotypeRate >= geno_q[2] ~ "high",
      TRUE                        ~ "middle"
    )
  )

# ------------------------------
# Normalize names for joining
# ------------------------------
# Primers: remove -F/-R
primers <- primers %>%
  mutate(name_simple = str_replace(name_clean, "-[FR]$", ""))

# SNP summary: remove trailing _1 or similar
snp_summary <- snp_summary %>%
  mutate(Locus_simple = str_replace(Locus_ID, "_\\d+$", ""))

# ------------------------------
# Merge and combine quartile info
# ------------------------------
primers_quartiles <- primers %>%
  left_join(
    snp_summary %>% select(Locus_simple, Depth_Quartile, Genotype_Quartile),
    by = c("name_simple" = "Locus_simple")
  ) %>%
  # Create a single combined quartile column
  mutate(
    quartile = case_when(
      Depth_Quartile == "high" | Genotype_Quartile == "high" ~ "high",
      Depth_Quartile == "low"  | Genotype_Quartile == "low"  ~ "low",
      TRUE ~ "middle"
    )
  ) %>%
  select(-name_simple, -Depth_Quartile, -Genotype_Quartile)

# ------------------------------
# Save final CSV
# ------------------------------
write_csv(primers_quartiles, "primers_combined_quartiles.csv")


