#the purpose of this script is to clean the primer csv file associated with the
#eulachon gsi project in order to run PrimerPooler to determine if we need to split 
#our primers into two pools for PCR1. PrimerPooler determined that two pools were likely necessary
#and produced a list of 36 loci that were not aligning to the genome, thus I finish this script 
#by aligning the list of 36 primers to a genome fasta file to determine if they really would not
#align. This showed they did align, so I will double back and check my work with PrimerPooler. 

# load packages
library(readr)
library(dplyr)
library(stringr)

### Read data in from raw csvs; one contains primer info, the other contains a list 
#of loci no longer used that we want to remove
primerlist <- read_csv("primerlist.csv")
toremove   <- read_csv("toremove.csv", col_names = FALSE)

# Extract vector of IDs to remove
remove_ids <- toremove$X1

# Add a column extracting the prefix ("Tpa_####") from column G (7th column)
primerlist <- primerlist %>%
  mutate(prefix = str_extract(.[[7]], "^Tpa_\\d+"))

# Filter out rows whose prefix is in the removal list
filtered <- primerlist %>%
  filter(!prefix %in% remove_ids) %>%
  select(-prefix)

# Save result; this csv now represents our current pool of loci 
write_csv(filtered, "filtered_primers.csv")

### Now we want to clean the data to develop a fasta file for use with PrimerPooler
# Read filtered primer list
primers <- read_csv("filtered_primers.csv")

#Remove spaces from sequences in column H (8th column)
primers <- primers %>%
  mutate(seq_clean = gsub(" ", "", .[[8]]))

#Clean primer names in column G (7th column)
# Convert "Tpa_1234_456F" → "Tpa_1234-F"
primers <- primers %>%
  mutate(name_clean = str_replace(.[[7]], "^(Tpa_\\d+)_\\d+(F|R)$", "\\1-\\2"))

# Save cleaned table; this now has two columns "seq_clean" and "name_clean" that we will use to 
#build our fasta file
write_csv(primers, "primers_cleaned.csv")

#Now lets create our fasta file for use with PrimerPooler
# Read cleaned primer list
cleanedprimers <- read_csv("primers_cleaned.csv")

# Adjust: T = 20th col, U = 21st col
name_col <- 21  # U
seq_col  <- 20  # T

# Build FASTA content
fasta_lines <- c()
for (i in 1:nrow(cleanedprimers)) {
  header <- paste0(">", cleanedprimers[[i, name_col]])
  sequence <- cleanedprimers[[i, seq_col]]
  fasta_lines <- c(fasta_lines, header, sequence)
}

# Write all primers into one fasta file
writeLines(fasta_lines, "all_primers.fasta") #this fasta file contains our needed data

#BUT - we have not removed the illumina adapters from our sequences. This needs to be done for 
#PrimerPooler to work appropriately 

# Define adapters to remove; these are illumina adapters
adapters <- c(
  "CGACAGGTTCAGAGTTCTACAGTCCGACGATC",
  "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
)

# Function to remove adapters from each sequence in the "primers_cleaned.csv"
remove_adapters <- function(sequence, adapters) {
  for (ad in adapters) {
    sequence <- gsub(ad, "", sequence, fixed = TRUE)
  }
  return(sequence)
}

# Updated "cleanedprimers" and apply adapter removal to sequences in column T (20th column)
cleanedprimers <- cleanedprimers %>%
  mutate(seq_clean = sapply(.[[20]], remove_adapters, adapters = adapters))

# Optional: inspect cleaned sequences
head(cleanedprimers$seq_clean)

# Export all primers to a single FASTA file
name_col <- 21  # Column U, primer names
seq_col  <- "seq_clean"  # the cleaned sequences

fasta_lines <- c()
for (i in 1:nrow(cleanedprimers)) {
  header <- paste0(">", cleanedprimers[[i, name_col]])
  sequence <- cleanedprimers[[i, seq_col]]
  fasta_lines <- c(fasta_lines, header, sequence)
}

writeLines(fasta_lines, "all_primers_cleaned.fasta")
### we now have a fasta file that contains only the primers used in this current pool, with cleaned
#names and cleaned matching sequences that have the adapters removed. This file is ready to be
#plugged into PrimerPooler

#PrimerPooler returned 36 primers that would not align with the genome - let's double check that. 

###Are the 36 primers that were not ID'd by primerpooler not aligning with the genome?###
#identify 36 primers that did not match
# Create a vector of all primer names to search for
primer_names <- c(
  "Tpa_5636-F","Tpa_5636-R",
  "Tpa_6422-F","Tpa_6422-R",
  "Tpa_6016-F","Tpa_6016-R",
  "Tpa_6451-F","Tpa_6451-R",
  "Tpa_6220-F","Tpa_6220-R",
  "Tpa_5441-F","Tpa_5441-R",
  "Tpa_4506-F","Tpa_4506-R",
  "Tpa_5196-F","Tpa_5196-R",
  "Tpa_5360-F","Tpa_5360-R",
  "Tpa_4549-F","Tpa_4549-R",
  "Tpa_5540-F","Tpa_5540-R",
  "Tpa_4731-F","Tpa_4731-R",
  "Tpa_7848-F","Tpa_7848-R",
  "Tpa_2854-F","Tpa_2854-R",
  "Tpa_3679-F","Tpa_3679-R",
  "Tpa_4291-F","Tpa_4291-R",
  "Tpa_4230-F","Tpa_4230-R",
  "Tpa_3248-F","Tpa_3248-R",
  "Tpa_3879-F","Tpa_3879-R",
  "Tpa_4055-F","Tpa_4055-R",
  "Tpa_1703-F","Tpa_1703-R",
  "Tpa_1891-F","Tpa_1891-R",
  "Tpa_2006-F","Tpa_2006-R",
  "Tpa_2412-F","Tpa_2412-R",
  "Tpa_7652-F","Tpa_7652-R",
  "Tpa_7597-F","Tpa_7597-R",
  "Tpa_7757-F","Tpa_7757-R",
  "Tpa_7240-F","Tpa_7240-R",
  "Tpa_6937-F","Tpa_6937-R",
  "Tpa_7367-F","Tpa_7367-R",
  "Tpa_7604-F","Tpa_7604-R",
  "Tpa_7690-F","Tpa_7690-R",
  "Tpa_1042-F","Tpa_1042-R",
  "Tpa_0650-F","Tpa_0650-R",
  "Tpa_0569-F","Tpa_0569-R",
  "Tpa_1098-F","Tpa_1098-R"
)

# Filter the primers_cleaned data for these names
primers_notfound <- cleanedprimers %>%
  filter(name_clean %in% primer_names)

# Optional: write to a CSV
write_csv(primers_notfound, "selected_primers.csv")

### now compare to genome ###
library(Biostrings)

# Load the genome fasta
genome <- readDNAStringSet("Tpac_genome.fna")
primer_notfound_seqs <- DNAStringSet(primers_notfound$seq_clean)
primer_notfound_names <- primers_notfound$name_clean

# Check for exact matches in genome
matches <- lapply(primer_notfound_seqs, function(primer) {
  vmatchPattern(primer, genome)
})

# Identify primers with no matches
no_match <- sapply(matches, function(x) length(x) == 0)

# Output names of primers that cannot be aligned
primers_no_alignment <- primer_notfound_names[no_match]

# Print result
print(primers_no_alignment)

# Optional: save to CSV
write_csv(data.frame(primer = primers_no_alignment),
          "primers_no_alignment.csv")

## Now i want to create a new column in primers_cleaned.csv that identifies if a
#primer is in pool1 or pool2 according to PrimerPooler

# Function to read primer names from a pool file
read_pool <- function(file) {
  lines <- readLines(file)
  primers <- gsub("^>", "", lines[grepl("^>", lines)]) # remove ">" from headers
  return(primers)
}

pool1 <- read_pool("pool1.txt")
pool2 <- read_pool("pool2.txt")

check_pool <- function(primer) {
  in_pool1 <- primer %in% pool1
  in_pool2 <- primer %in% pool2
  
  if (in_pool1 && in_pool2) {
    return("both")
  } else if (in_pool1) {
    return("pool1")
  } else if (in_pool2) {
    return("pool2")
  } else {
    return("none")
  }
}

# Create a new dataframe with pool assignments
primers_with_pools <- cleanedprimers %>%
  mutate(Pool_Assignment = sapply(.[[21]], check_pool))

# Save as a new CSV (original primers_cleaned.csv stays the same)
write_csv(primers_with_pools, "primers_with_pools.csv")
