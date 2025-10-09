# R script: Primer Reads by Locus colored by primer-pool quartile
# Requires: tidyverse (or at least readr, dplyr, ggplot2)
# Files (adjust paths if different):
# - /mnt/data/GTscore_locusSummary.txt
# - /mnt/data/primers_combined_quartiles.csv

library(tidyverse)

# 1) Read GTscore locus summary (tab-delimited)
gt <- read.delim("/mnt/data/GTscore_locusSummary.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Quick check of column names
# print(colnames(gt))
# Expected: "Locus", "Primer Reads", "Primer Probe Reads", "Primer Probe Proportion"

# 2) Read primer quartiles CSV
primers <- read.csv("/mnt/data/primers_combined_quartiles.csv", stringsAsFactors = FALSE)

# Quick check:
# print(head(primers))
# Expected columns include: name_clean, quartile

# 3) Normalize primer names in primers: remove trailing -R or -F
# Create a Locus column matching GTscore style (e.g., "Tpa_6453")
primers <- primers %>%
  mutate(name_clean = as.character(name_clean),
         Locus = sub("(-[RF])$", "", name_clean, perl = TRUE) )  # strips -R or -F at end

# 4) Prepare gt table: make sure Primer Reads is numeric and column name is a safe name
# If column name contains spaces like "Primer Reads", refer with backticks.
gt <- gt %>%
  mutate(`Primer Reads` = as.numeric(`Primer Reads`))

# 5) Join quartile information onto GTscore table by Locus
gt_joined <- gt %>%
  left_join(select(primers, Locus, quartile), by = "Locus") %>%
  mutate(quartile = ifelse(is.na(quartile), "unknown", quartile),
         quartile = factor(quartile, levels = c("high", "middle", "low", "unknown")))

# 6) Color mapping: high -> green, middle -> orange, low -> red, unknown -> grey
fill_colors <- c("high"   = "green",
                 "middle" = "orange",
                 "low"    = "red",
                 "unknown"= "grey70")

# 7) Plot: bar chart ordered by Primer Reads descending
library(ggplot2)
p <- ggplot(gt_joined, aes(x = reorder(Locus, -`Primer Reads`), y = `Primer Reads`, fill = quartile)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = fill_colors, na.value = "grey70", name = "Quartile") +
  labs(title = "Primer Reads by Locus (colored by primer-pool quartile)",
       subtitle = "Loci colored by primer concentration quartile used in primer pool",
       x = "Locus (ordered by Primer Reads)",
       y = "Primer Reads",
       caption = "GTscore locus summary and primers_combined_quartiles.csv used") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "right")

# Optional: use log10 y-axis to handle large dynamic range (uncomment if wanted)
p <- p + scale_y_log10()

# Print plot
print(p)

# Optionally save:
ggsave("/mnt/data/PrimerReads_by_Locus_quartile.png", p, width = 16, height = 6, dpi = 300)