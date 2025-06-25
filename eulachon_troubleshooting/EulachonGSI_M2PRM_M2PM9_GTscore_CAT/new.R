---
  title: "EulachonGSI_M2PRM_M2PM9_GTscore_CAT"
author: "MB"
date: "2024-06-18"
output: html_document
---
  
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r}
# Create output directories
if (!dir.exists("QCplots")) dir.create("QCplots")
if (!dir.exists("scatterPlots_loci")) dir.create("scatterPlots_loci")
if (!dir.exists("scatterPlots_samples")) dir.create("scatterPlots_samples")

# Load GTscore functions
stopifnot(file.exists("GTscore.R"))
source("GTscore.R")

# Load libraries
library(tidyverse)
library(msa)
library(pbapply)
library(BiocGenerics)
library(Biostrings)
library(BiocManager)
```

```{r}
# Read AmpliconReadCounter outputs
stopifnot(file.exists("LocusTable_singleSNPs.txt"))
stopifnot(file.exists("AlleleReads_singleSNPs.txt"))

singleSNP_locusTable <- read.delim("LocusTable_singleSNPs.txt", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
singleSNP_alleleReads <- read.delim("AlleleReads_singleSNPs.txt", header=TRUE, row.names=1, stringsAsFactors=FALSE, check.names=FALSE)

singleSNP_alleleReads_adjusted <- correctReads(singleSNP_locusTable, singleSNP_alleleReads)
write.table(singleSNP_alleleReads_adjusted, "singleSNP_alleleReads_adjusted.txt", quote=FALSE, sep="\t")
```

```{r}
# Genotyping
polyGenResults_singleSNP <- polyGen(singleSNP_locusTable, singleSNP_alleleReads_adjusted)
write.table(polyGenResults_singleSNP, "polyGenResults_singleSNP.txt", quote=FALSE, sep="\t")

# Clean column names
clean_names <- function(names_vec) basename(gsub("\\\\", "/", names_vec))
colnames(singleSNP_alleleReads_adjusted) <- clean_names(colnames(singleSNP_alleleReads_adjusted))
colnames(polyGenResults_singleSNP) <- clean_names(colnames(polyGenResults_singleSNP))

# Scatter plots
plotGenotypes_sample(singleSNP_locusTable, singleSNP_alleleReads_adjusted, polyGenResults_singleSNP, type='scatter', savePlot="Y", saveDir="scatterPlots_samples")
plotGenotypes(singleSNP_locusTable, singleSNP_alleleReads_adjusted, polyGenResults_singleSNP, type='scatter', savePlot="Y", saveDir="scatterPlots_loci")

# Summarize
singleSNP_summary <- summarizeGTscore(singleSNP_alleleReads, singleSNP_locusTable, polyGenResults_singleSNP)
write.table(singleSNP_summary, "singleSNP_summary.txt", quote=FALSE, sep="\t")
```

```{r}
# Sample summaries
GTscore_individualSummary <- read.delim("GTscore_individualSummary.txt", header=TRUE, stringsAsFactors=FALSE) %>% 
  mutate_if(is.numeric, round, digits=3)

singleSNP_sampleSummary <- summarizeSamples(polyGenResults_singleSNP, singleSNP_alleleReads_adjusted) %>% 
  mutate_if(is.numeric, round, digits=3)

write.table(singleSNP_sampleSummary, "singleSNP_sampleSummary.txt", quote=FALSE, sep="\t")
write.table(summary(singleSNP_sampleSummary), "singleSNP_sampleSummaryTable.txt", quote=FALSE, sep="\t")

GTscore_individualSummary <- merge(GTscore_individualSummary, singleSNP_sampleSummary, by.x="Sample", by.y="sample")
write.table(GTscore_individualSummary, "GTscore_individualSummary.txt", quote=FALSE, sep="\t")
write.table(summary(GTscore_individualSummary), "GTscore_individualSummaryTable.txt", quote=FALSE, sep="\t")
```

```{r}
# QA/QC thresholds
minGenotypeRate <- 0.73
minHeterozygosity <- 0.1
maxHeterozygosity <- 0.35
maxconScore <- 0.3
minconScore <- 0.01
matchThresh <- 0.8
commonThresh <- 0.6
```

```{r}
# Sample QC plots
library(ggplot2)

# Genotype rate
ggplot(GTscore_individualSummary) +
  geom_histogram(aes(x=GenotypeRate), binwidth=0.03) +
  geom_vline(xintercept=minGenotypeRate, lty="dashed") +
  labs(title="Sample Genotyping Rate", x="Genotype Rate", y="Count") +
  theme_bw()
ggsave("hist-SampleGTrateSNP.pdf", path="./QCplots", width=30, height=20, units="cm")

# Failures
sampleFails <- GTscore_individualSummary %>% filter(GenotypeRate <= minGenotypeRate)
write.table(format(sampleFails, digits=4, scientific=FALSE), "sampleFails.txt", quote=FALSE, sep="\t", row.names=FALSE)

# Other QC plots and filters omitted for brevity...
```

```{r}
# Problem sample identification
sampleWrongSpp <- GTscore_individualSummary %>% filter(Heterozygosity < minHeterozygosity)
sampleContaminants <- GTscore_individualSummary %>% filter(Heterozygosity >= maxHeterozygosity | conScore >= maxconScore)

sampleDups <- IDduplicateSamples(polyGenResults_singleSNP) %>% 
  filter(proportionMatch >= matchThresh & proportionCommon >= commonThresh) %>% 
  mutate_if(is.numeric, round, digits=3)

write.table(sampleWrongSpp, "sampleWrongSpp.txt", quote=FALSE, sep="\t", row.names=FALSE)
write.table(sampleContaminants, "sampleContaminants.txt", quote=FALSE, sep="\t", row.names=FALSE)
write.table(sampleDups, "sampleDuplicates.txt", quote=FALSE, sep="\t", row.names=FALSE)

# Blacklist
sampleBlacklist <- sampleFails
write.table(format(sampleBlacklist, digits=4, scientific=FALSE), "sampleBlacklist.txt", quote=FALSE, sep="\t", row.names=FALSE)
```

```{r}
# Clean data
GTscore_individualSummary_cleaned <- anti_join(GTscore_individualSummary, sampleBlacklist)
write.table(GTscore_individualSummary_cleaned, "GTscore_individualSummary_cleaned.txt", quote=FALSE, sep="\t")
write.table(summary(GTscore_individualSummary_cleaned), "GTscore_individualSummaryTable_cleaned.txt", quote=FALSE, sep="\t")
```

<!-- You can continue cleaning the rest of the read count plots section if needed. -->
  