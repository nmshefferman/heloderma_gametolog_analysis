# ------
# Goal
#   Build summary tables identifying Z-W PAR best hits.
#   Then, identify Z-W PAR pairs based on whether the Z-W and W-Z best hits are reciprocal.
# Contents
#   1. Clean up data
#   2. W-Z best hits
#   3. Z-W best hits
#   4. Z-W PAR pairs
# ------

# Load in packages
library(tidyverse)

# Set working directory
setwd("/vf/users/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/PAR")

# === 1) Clean up data ===
# Read in BLAST tables
z_to_w <- read_tsv("blast_results/PAR_zGenes_wGenes.tsv", comment = "#", col_names = FALSE)
w_to_z <- read_tsv("blast_results/PAR_wGenes_zGenes.tsv", comment = "#", col_names = FALSE)

# Name columns
colnames(z_to_w)[1:12] <- c("Z", "W", "Percent identity", "Alignment length", "Number mismatch", "Number gap openings", "Start alignment query", "End alignment query", "Start alignment subject", "End alignment subject", "Evalue", "Bit score")
colnames(w_to_z)[1:12] <- c("W", "Z", "Percent identity", "Alignment length", "Number mismatch", "Number gap openings", "Start alignment query", "End alignment query", "Start alignment subject", "End alignment subject", "Evalue", "Bit score")

# Convert transcript IDs to gene IDs
z_to_w <- z_to_w %>%
  mutate(Z = str_remove(Z, "-mRNA-.*"),
         W = str_remove(W, "-mRNA-.*"),
         # Replace LOC_ with LOCW_ for W genes
         W = str_replace(W, "^LOC_", "LOCW_"))
w_to_z <- w_to_z %>%
  mutate(W = str_remove(W, "-mRNA-.*"),
         Z = str_remove(Z, "-mRNA-.*"),
         # Replace LOC_ with LOCW_ for W genes
         W = str_replace(W, "^LOC_", "LOCW_"))

# Read in gene coordinate files
z_coords <- read_tsv("gene_coords/Z_PAR_coords.gff3", comment = "#", col_names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"))
w_coords <- read_tsv("gene_coords/W_PAR_coords.gff3", comment = "#", col_names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"))

# Extract gene IDs from attributes column
z_coords <- z_coords %>% 
  mutate(Z_gene = str_extract(attributes, "LOC_[0-9]+")) %>% 
  select(Z_gene, Z_start = start, Z_end = end)
w_coords <- w_coords %>% 
  mutate(W_gene = str_extract(attributes, "LOCW_[0-9]+")) %>% 
  select(W_gene, W_start = start, W_end = end)

# === 2) W-Z best hits ===
z_best <- z_to_w %>% 
  group_by(Z) %>% 
  slice_min(order_by = Evalue, n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  select(Z, W_best = W, z_w_evalue = Evalue)

# Top two Z hits for each W gene
w_z_hits <- w_to_z %>% 
  group_by(W) %>% 
  arrange(Evalue, desc(`Alignment length`), .by_group = TRUE) %>% 
  #distinct(Z, .keep_all = TRUE) %>% 
  mutate(rank = row_number()) %>% 
  filter(rank <= 2) %>% 
  ungroup() %>% 
  left_join(z_best, by = "Z") %>% 
  mutate(reciprocal = (W == W_best))

# Summarize top two Z hits per W gene
w_z_summary <- w_z_hits %>% 
  group_by(W) %>% 
  summarize(
    Z_hit_1 = first(Z),
    Z_hit_1_evalue = first(Evalue),
    Z_hit_1_length = first(`Alignment length`),
    Z_hit_1_reciprocal = first(reciprocal),
    Z_hit_2 = if (n() >= 2) nth(Z, 2) else NA_character_,
    Z_hit_2_evalue = if (n() >= 2) nth(Evalue, 2) else NA_real_,
    Z_hit_2_length = if (n() >= 2) nth(`Alignment length`, 2) else NA_real_,
    Z_hit_2_reciprocal = if (n() >= 2) nth(reciprocal, 2) else NA, 
    .groups = "drop") 

# Final table
w_z_summary <- w_z_summary %>%
  left_join(w_coords, by = c("W" = "W_gene")) %>% 
  left_join(z_coords, by = c("Z_hit_1" = "Z_gene")) %>% 
  left_join(z_coords, by = c("Z_hit_2" = "Z_gene")) %>% 
  rename(
    Z_start_1 = Z_start.x,
    Z_end_1 = Z_end.x,
    Z_start_2 = Z_start.y,
    Z_end_2 = Z_end.y)

# Re-arrange columns
w_z_summary <- w_z_summary[, c(1, 10, 11, 2, 12, 13, 3, 4, 5, 6, 14, 15, 7, 8, 9)]

# === 3) Z-W best hits ===
# Best Z hit for each W gene
w_best <- w_to_z %>% 
  group_by(W) %>% 
  slice_min(order_by = Evalue, n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  select(W, Z_best = Z, w_z_evalue = Evalue)

# Top two W hits for each Z gene
z_w_hits <- z_to_w %>% 
  group_by(Z) %>% 
  arrange(Evalue, desc(`Alignment length`), .by_group = TRUE) %>% 
  #distinct(W, .keep_all = TRUE) %>% 
  mutate(rank = row_number()) %>% 
  filter(rank <= 2) %>% 
  ungroup() %>% 
  left_join(w_best, by = "W") %>% 
  mutate(reciprocal = (Z == Z_best))

# Summarize top two W hits per Z gene
z_w_summary <- z_w_hits %>% 
  group_by(Z) %>% 
  summarize(
    W_hit_1 = first(W),
    W_hit_1_evalue = first(Evalue),
    W_hit_1_length = first(`Alignment length`),
    W_hit_1_reciprocal = first(reciprocal),
    W_hit_2 = if (n() >= 2) nth(W, 2) else NA_character_,
    W_hit_2_evalue = if (n() >= 2) nth(Evalue, 2) else NA_real_,
    W_hit_2_length = if (n() >= 2) nth(`Alignment length`, 2) else NA_real_,
    W_hit_2_reciprocal = if (n() >= 2) nth(reciprocal, 2) else NA, 
    .groups = "drop")

# Final table
z_w_summary <- z_w_summary %>%
  left_join(z_coords, by = c("Z" = "Z_gene")) %>% 
  left_join(w_coords, by = c("W_hit_1" = "W_gene")) %>% 
  left_join(w_coords, by = c("W_hit_2" = "W_gene")) %>% 
  rename(
    W_start_1 = W_start.x,
    W_end_1 = W_end.x,
    W_start_2 = W_start.y,
    W_end_2 = W_end.y)

# Re-arrange columns
z_w_summary <- z_w_summary[, c(1, 10, 11, 2, 12, 13, 3, 4, 5, 6, 14, 15, 7, 8, 9)]

# === 4) Z-W PAR pairs ===
# Using each table of best hits, this table takes only reciprocal best hits. These are the most high-confidence gametologs.
w_z_summary_v2 <- w_z_summary %>%
  mutate(best_Z_gene = Z_hit_1) %>%
  # Join the Z-W table to compare the best W hit for each Z gene that is kept
  left_join(z_w_summary %>%
              select(Z, best_W_from_Z = W_hit_1),
            by = c("best_Z_gene" = "Z")) %>%
  # Identify reciprocal best-hit pairs (TRUE = W's best Z hit also identifies this W as its best W hit)
  mutate(reciprocal = !is.na(best_Z_gene) & (W == best_W_from_Z))

# Read in ref gff file to add gene names to final gametolog table
ref <- read_tsv("../../references/Heloderma_suspectum.final.reference.gff3.gz", comment = "#", col_names = c("seqid","source","type","start","end","score","strand","phase","attributes"), show_col_types = FALSE)

# Pull out gene names 
names_z <- ref %>%
  filter(type == "gene", seqid == "chrZ") %>%
  transmute(
    Z = str_extract(attributes, "LOC_[0-9]+"),
    Z_name = str_extract(attributes, "(?<=Note=Similar to ).*?(?= \\()")) %>%
  distinct()
names_w <- ref %>%
  filter(type == "gene", seqid == "chrW") %>%
  transmute(
    W = str_extract(attributes, "LOCW_[0-9]+"),
    W_name = str_extract(attributes, "(?<=Note=Similar to ).*?(?= \\()")) %>%
  distinct()

# Keep only true reciprocal best hits
Z_W_PARs <- w_z_summary_v2 %>%
  filter(reciprocal == TRUE) %>%
  select(W = best_W_from_Z, W_start, W_end, Z = best_Z_gene, Z_start = Z_start_1, Z_end = Z_end_1) %>% 
  # Join name columns to gametolog table
  left_join(names_w, by = "W") %>% 
  left_join(names_z, by = "Z") %>% 
  # Filter out any genes NOT in the PARs of either chromosome
  filter(Z_start <= 2209992, W_start <= 2205462) %>% 
  # Reorder columns
  select(c(1, 7, 2, 3, 4, 8, 5, 6))

# Save table
write_tsv(Z_W_PARs, "Z_W_PARs.tsv")
