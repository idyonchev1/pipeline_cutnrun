library(tidyr)
library(dplyr)

# Read the count data
counts <- read.csv(snakemake@input[["counts"]], header=TRUE, sep="\t")
colnames(counts) <- gsub("^X\\.", "", colnames(counts))
colnames(counts) <- gsub("\\.filtered\\.$", "", colnames(counts))
print("Counts data:")
print(head(counts))
print(dim(counts))

# Process counts
counts <- counts %>%
  unite("combined", .chr., start., end., sep = "_")  # Updated column names
row.names(counts) <- counts$combined
counts <- counts %>% select(-combined)
print("Processed counts:")
print(head(counts))
print(dim(counts))

# Calculate total counts and normalization factors
sf <- colSums(counts)
sf_df <- data.frame(Sample = names(sf), sf = sf, normalizer = 1/sf)

# Clean up sample names in sf_df
sf_df$Sample <- gsub("^X\\.|\\.$", "", sf_df$Sample)

print("Size factors:")
print(head(sf_df))

# Read sample information
samples <- read.csv(snakemake@input[["samples"]])
print("Sample information:")
print(head(samples))
print(dim(samples))

# Merge size factors with sample information
sf_merged <- merge(sf_df, samples, by = "Sample")
print("Merged size factors and sample information:")
print(head(sf_merged))
print(dim(sf_merged))

# Calculate IgG normalizers
igg_normalizers <- sf_merged %>%
  filter(condition == "IgG") %>%
  select(replicate, line, igg_normalizer = normalizer)
print("IgG normalizers:")
print(igg_normalizers)

# Normalize to IgG
sf_normalized <- sf_merged %>%
  left_join(igg_normalizers, by = c("replicate", "line")) %>%
  mutate(normalized_to_igg = normalizer / igg_normalizer) %>%
  select(-igg_normalizer)
print("Normalized to IgG:")
print(head(sf_normalized))
print(dim(sf_normalized))

# Create export name
sf_normalized$export_name <- paste0(sf_normalized$line, "_", sf_normalized$condition, "_rep", sf_normalized$replicate)

# Write the output
write.table(sf_normalized, file = snakemake@output[["normalized"]], sep = "\t", quote = FALSE, row.names = FALSE)
