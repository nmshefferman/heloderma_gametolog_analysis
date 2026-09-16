# ------
# Goal
#   Build a table that lists the names of the genes that are are orthologs to the same Gila Z genes.
#   This table will be used to make the MSA files needed to run IQ-Tree.
# Contents
#   1) Compare all orthologs to the same Gila Z gene across all species
#   2) Build table of orthologs to the shared Heloderma/Varanus Z-W genes
#   3) Separate combined table into analyses subsets
#     a) Make a table that has genes that are present in all species
#     b) Make a table that has genes that are in all species except Shinisaurus
# ------

# Load in packages
library(tidyverse)

# Set working directory
setwd("/vf/users/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/")

# === 1) Compare all orthologs to the same Gila Z gene across all species ===
# Read in Heloderma gametolog files
gila_z <- read_tsv("heloderma_only/heloderma_Z_W_gametologs_final.tsv", show_col_types = FALSE) %>% 
  rename(Heloderma = Z) %>% 
  left_join(read_tsv("outgroups/transcript_lengths/heloderma_transcript_lengths.tsv", col_names = c("Heloderma_transcript", "Heloderma", "Heloderma_len"), show_col_types = FALSE),
            by = "Heloderma") %>% 
  group_by(Heloderma) %>%
  slice_max(order_by = Heloderma_len, n = 1, with_ties = FALSE) %>%
  select(Heloderma, Heloderma_transcript)

# Join gila W genes
gila_w <- read_tsv("heloderma_only/heloderma_Z_W_gametologs_final.tsv", show_col_types = FALSE) %>% 
  rename(Heloderma = Z, Heloderma_W = W) %>% 
  left_join(read_tsv("outgroups/transcript_lengths/heloderma_w_transcript_lengths.tsv", col_names = c("Heloderma_W_transcript", "Heloderma_W", "Heloderma_W_len"), show_col_types = FALSE),
            by = "Heloderma_W") %>% 
  group_by(Heloderma) %>% 
  slice_max(order_by = Heloderma_W_len, n = 1, with_ties = FALSE) %>%
  select(Heloderma, Heloderma_W_transcript)

# Read in outgroup files
aniella <- read_tsv("outgroups/orthologs/anniella_ortho.tsv", show_col_types = FALSE) %>%
  # Keep only unique combinations of genes
  distinct(Heloderma, Anniella) %>% 
  left_join(read_tsv("outgroups/transcript_lengths/anniella_transcript_lengths.tsv", col_names = c("Anniella_transcript", "Anniella", "Anniella_len"), show_col_types = FALSE),
            by = "Anniella") %>%
  group_by(Heloderma) %>%
  slice_max(order_by = Anniella_len, n = 1, with_ties = FALSE) %>%
  select(Heloderma, Anniella_transcript)

anolis <- read_tsv("outgroups/orthologs/anolis_ortho.tsv", show_col_types = FALSE) %>%
  separate_rows(Anolis, sep = ",\\s*") %>% 
  distinct(Heloderma, Anolis) %>% 
  left_join(read_tsv("outgroups/transcript_lengths/anolis_transcript_lengths.tsv", col_names = c("Anolis_transcript", "Anolis", "Anolis_len"), show_col_types = FALSE),
            by = "Anolis") %>%
  group_by(Heloderma) %>%
  slice_max(order_by = Anolis_len, n = 1, with_ties = FALSE) %>%
  select(Heloderma, Anolis_transcript)

candoia <- read_tsv("outgroups/orthologs/candoia_ortho.tsv", show_col_types = FALSE) %>%
  separate_rows(Candoia, sep = ",\\s*") %>% 
  distinct(Heloderma, Candoia) %>% 
  left_join(read_tsv("outgroups/transcript_lengths/candoia_transcript_lengths.tsv", col_names = c("Candoia_transcript", "Candoia", "Candoia_len"), show_col_types = FALSE),
            by = "Candoia") %>%
  group_by(Heloderma) %>%
  slice_max(order_by = Candoia_len, n = 1, with_ties = FALSE) %>%
  select(Heloderma, Candoia_transcript)

elgaria <- read_tsv("outgroups/orthologs/elgaria_ortho.tsv", show_col_types = FALSE) %>%
  separate_rows(Elgaria, sep = ",\\s*") %>% 
  distinct(Heloderma, Elgaria) %>% 
  left_join(read_tsv("outgroups/transcript_lengths/elgaria_transcript_lengths.tsv", col_names = c("Elgaria_transcript", "Elgaria", "Elgaria_len"), show_col_types = FALSE),
            by = "Elgaria") %>%
  group_by(Heloderma) %>%
  slice_max(order_by = Elgaria_len, n = 1, with_ties = FALSE) %>%
  select(Heloderma, Elgaria_transcript)

furcifer <- read_tsv("outgroups/orthologs/furcifer_ortho.tsv", show_col_types = FALSE) %>%
  separate_rows(Furcifer, sep = ",\\s*") %>% 
  distinct(Heloderma, Furcifer) %>% 
  left_join(read_tsv("outgroups/transcript_lengths/furcifer_transcript_lengths.tsv", col_names = c("Furcifer_transcript", "Furcifer", "Furcifer_len"), show_col_types = FALSE),
            by = "Furcifer") %>%
  group_by(Heloderma) %>%
  slice_max(order_by = Furcifer_len, n = 1, with_ties = FALSE) %>%
  select(Heloderma, Furcifer_transcript)

podarcis <- read_tsv("outgroups/orthologs/podarcis_ortho.tsv", show_col_types = FALSE) %>%
  separate_rows(Podarcis, sep = ",\\s*") %>% 
  distinct(Heloderma, Podarcis) %>% 
  left_join(read_tsv("outgroups/transcript_lengths/podarcis_transcript_lengths.tsv", col_names = c("Podarcis_transcript", "Podarcis", "Podarcis_len"), show_col_types = FALSE),
            by = "Podarcis") %>%
  group_by(Heloderma) %>%
  slice_max(order_by = Podarcis_len, n = 1, with_ties = FALSE) %>%
  select(Heloderma, Podarcis_transcript)

shinisaurus <- read_tsv("outgroups/orthologs/shinisaurus_ortho.tsv", show_col_types = FALSE) %>%
  separate_rows(Shinisaurus, sep = ",\\s*") %>% 
  distinct(Heloderma, Shinisaurus) %>% 
  left_join(read_tsv("outgroups/transcript_lengths/shinisaurus_transcript_lengths.tsv", col_names = c("Shinisaurus_transcript", "Shinisaurus", "Shinisaurus_len"), show_col_types = FALSE),
            by = "Shinisaurus") %>%
  group_by(Heloderma) %>%
  slice_max(order_by = Shinisaurus_len, n = 1, with_ties = FALSE) %>%
  select(Heloderma, Shinisaurus_transcript)

thamnophis <- read_tsv("outgroups/orthologs/thamnophis_ortho.tsv", show_col_types = FALSE) %>%
  separate_rows(Thamnophis, sep = ",\\s*") %>% 
  distinct(Heloderma, Thamnophis) %>% 
  left_join(read_tsv("outgroups/transcript_lengths/thamnophis_transcript_lengths.tsv", col_names = c("Thamnophis_transcript", "Thamnophis", "Thamnophis_len"), show_col_types = FALSE),
            by = "Thamnophis") %>%
  group_by(Heloderma) %>%
  slice_max(order_by = Thamnophis_len, n = 1, with_ties = FALSE) %>%
  select(Heloderma, Thamnophis_transcript)

varanus <- read_tsv("outgroups/orthologs/varanus_ortho.tsv", show_col_types = FALSE) %>%
  distinct(Heloderma, Varanus) %>% 
  left_join(read_tsv("outgroups/transcript_lengths/varanus_transcript_lengths.tsv", col_names = c("Varanus_transcript", "Varanus", "Varanus_len"), show_col_types = FALSE),
            by = "Varanus") %>%
  group_by(Heloderma) %>%
  slice_max(order_by = Varanus_len, n = 1, with_ties = FALSE) %>%
  select(Heloderma, Varanus, Varanus_transcript)

# Read in Varanus W file and format properly 
varanus_w <- read_tsv("varanus_only/varanus_Z_W_gametologs_final.tsv", show_col_types = FALSE) %>% 
  rename(Varanus = Z, Varanus_W = W) %>% 
  mutate(Varanus = str_remove(Varanus, "Z")) %>% 
  distinct(Varanus, Varanus_W) %>% 
  left_join(read_tsv("outgroups/transcript_lengths/varanus_w_transcript_lengths.tsv", col_names = c("Varanus_W_transcript", "Varanus_W", "Varanus_W_len"), show_col_types = FALSE),
            by = "Varanus_W") %>% 
  group_by(Varanus) %>% 
  slice_max(order_by = Varanus_W_len, n = 1, with_ties = FALSE) %>%
  select(Varanus, Varanus_W_transcript) 

# === 2) Build table of orthologs to the shared Heloderma/Varanus Z-W genes === 
# Join everything into one table, only strictly keeping all Gila and Varanus Z/W genes. Any other species can have an N/A if they don't have the ortholog.  
table <- gila_z %>% 
  inner_join(gila_w, by = "Heloderma") %>% 
  inner_join(varanus, by = "Heloderma") %>% 
  inner_join(varanus_w, by = "Varanus") %>% 
  left_join(aniella, by = "Heloderma") %>% 
  left_join(anolis, by = "Heloderma") %>% 
  left_join(candoia, by = "Heloderma") %>% 
  left_join(elgaria, by = "Heloderma") %>% 
  left_join(furcifer, by = "Heloderma") %>% 
  left_join(podarcis, by = "Heloderma") %>% 
  left_join(shinisaurus, by = "Heloderma") %>% 
  left_join(thamnophis, by = "Heloderma") %>% 
  ungroup() %>% 
  select(-Heloderma, -Varanus) %>% 
  rename(Heloderma_Z = Heloderma_transcript, Heloderma_W = Heloderma_W_transcript, Varanus_Z = Varanus_transcript, Varanus_W = Varanus_W_transcript, Anniella = Anniella_transcript, Anolis = Anolis_transcript, Candoia = Candoia_transcript, Elgaria = Elgaria_transcript, Furcifer = Furcifer_transcript, Podarcis = Podarcis_transcript, Shinisaurus = Shinisaurus_transcript, Thamnophis = Thamnophis_transcript) %>% 
  # Rename Varanus_Z gene names to match FASTA headers (LOCZ_ instead of LOC_)
  mutate(Varanus_Z = str_replace(Varanus_Z, "^LOC_", "LOCZ_"))

# Save table
write_tsv(table, "iqtree/combined_table.tsv")

# === 3) Separate combined table into analyses subsets ===
# a) Make a table that has genes that are present in all species
all_species <- table %>%
  filter(if_all(c(Heloderma_Z, Heloderma_W, Varanus_Z, Varanus_W, Anniella, Anolis, Candoia, Elgaria, Furcifer, Podarcis, Shinisaurus, Thamnophis), 
                ~ !is.na(.)))

# Save table 
write_tsv(all_species, "iqtree/combined_table_all_species.tsv")
  
# b) Make a table that has genes that are in all species except Shinisaurus
no_shinisaurus <- table %>%
  filter(if_all(c(Heloderma_Z, Heloderma_W, Varanus_Z, Varanus_W, Anniella, Anolis, Candoia, Elgaria, Furcifer, Podarcis, Thamnophis), 
                ~ !is.na(.))) %>% 
  filter(is.na(Shinisaurus))

# Save table
write_tsv(no_shinisaurus, "iqtree/combined_table_no_shinisaurus.tsv")
