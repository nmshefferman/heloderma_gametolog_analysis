# ------
# Goal
#   This script will use the list of Heloderma Z gametologs to find orthologs in all outgroup Orthofinder tables
# Contents 
#   1) Anniella stebbinsi
#   2) Anolis sagrei
#   3) Candoia aspera
#   4) Elgaria multicarinata
#   5) Furcifer pardalis
#   6) Podarcis raffonei
#   7) Shinisaurus crocodilurus
#   8) Thamnophis elegans
#   9) Varanus acanthurus
# ------

# Load in packages
library(tidyverse)

# Set working directory
setwd("/vf/users/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/")

# Read in list of Heloderma gametologs
gila_z <- read_tsv("heloderma_only/heloderma_Z_W_gametologs_final.tsv", show_col_types = FALSE) %>% 
  distinct(Z)

# === 1) Anniella stebbinsi ===
# Read in Heloderma v. Anniella table
anniella <- read_tsv("../Heloderma_orthofinder/Heloderma_v_Anniella.tsv", show_col_types = FALSE) %>% 
  # Remove -mRNA labels from gene names (so that it matches up with names in Heloderma Z list)
  mutate(Heloderma = str_remove(Heloderma, "-mRNA-.*"),
         Anniella = str_remove(Anniella, "-mRNA-.*"))

# Keep only genes that are in the Gila Z gametolog list
anniella_ortho <- anniella %>%
  semi_join(gila_z, by = c("Heloderma" = "Z")) %>% 
  distinct(Heloderma, Anniella)

# Save list
write_tsv(anniella_ortho, "outgroups/orthologs/anniella_ortho.tsv")

# === 2) Anolis sagrei ===
anolis <- read_tsv("../Heloderma_orthofinder/Heloderma_v_Anolis.tsv", show_col_types = FALSE) %>% 
  # Remove -mRNA labels from gene names (so that it matches up with names in Heloderma Z list)
  mutate(Heloderma = str_remove(Heloderma, "-mRNA-.*"))

# Keep only genes that are in the Gila Z gametolog list
anolis_ortho <- anolis %>%
  semi_join(gila_z, by = c("Heloderma" = "Z")) %>% 
  distinct(Heloderma, Anolis)

# Save list
write_tsv(anolis_ortho, "outgroups/orthologs/anolis_ortho.tsv")

# === 3) Candoia aspera ===
candoia <- read_tsv("../Heloderma_orthofinder/Heloderma_v_Candoia.tsv", show_col_types = FALSE) %>% 
  # Remove -mRNA labels from gene names (so that it matches up with names in Heloderma Z list)
  mutate(Heloderma = str_remove(Heloderma, "-mRNA-.*"))

# Keep only genes that are in the Gila Z gametolog list
candoia_ortho <- candoia %>%
  semi_join(gila_z, by = c("Heloderma" = "Z")) %>% 
  distinct(Heloderma, Candoia)

# Save list
write_tsv(candoia_ortho, "outgroups/orthologs/candoia_ortho.tsv")

# === 4) Elgaria multicarinata ===
# Read in Heloderma v. Elgaria table
elgaria <- read_tsv("../Heloderma_orthofinder/Heloderma_v_Elgaria.tsv", show_col_types = FALSE) %>% 
  # Remove -mRNA labels from gene names (so that it matches up with names in Heloderma Z list)
  mutate(Heloderma = str_remove(Heloderma, "-mRNA-.*"))

# Keep only genes that are in the Gila Z gametolog list
elgaria_ortho <- elgaria %>%
  semi_join(gila_z, by = c("Heloderma" = "Z")) %>% 
  distinct(Heloderma, Elgaria)

# Save list
write_tsv(elgaria_ortho, "outgroups/orthologs/elgaria_ortho.tsv")

# === 5) Furcifer pardalis === 
furcifer <- read_tsv("../Heloderma_orthofinder/Heloderma_v_Furcifer.tsv", show_col_types = FALSE) %>% 
  # Remove -mRNA labels from gene names (so that it matches up with names in Heloderma Z list)
  mutate(Heloderma = str_remove(Heloderma, "-mRNA-.*"))

# Keep only genes that are in the Gila Z gametolog list
furcifer_ortho <- furcifer %>%
  semi_join(gila_z, by = c("Heloderma" = "Z")) %>% 
  distinct(Heloderma, Furcifer)

# Save list
write_tsv(furcifer_ortho, "outgroups/orthologs/furcifer_ortho.tsv")

# === 6) Podarcis raffonei ===
# Read in Heloderma v. Podarcis table
podarcis <- read_tsv("../Heloderma_orthofinder/Heloderma_v_Podarcis.tsv", show_col_types = FALSE) %>% 
  # Remove -mRNA labels from gene names (so that it matches up with names in Heloderma Z list)
  mutate(Heloderma = str_remove(Heloderma, "-mRNA-.*"))

# Keep only genes that are in the Gila Z gametolog list
podarcis_ortho <- podarcis %>%
  semi_join(gila_z, by = c("Heloderma" = "Z")) %>% 
  distinct(Heloderma, Podarcis)

# Save list
write_tsv(podarcis_ortho, "outgroups/orthologs/podarcis_ortho.tsv")

# === 7) Shinisaurus crocodilurus ===
# Read in Heloderma v. Shinisaurus table
shinisaurus <- read_tsv("../Heloderma_orthofinder/Heloderma_v_Shinisaurus.tsv", show_col_types = FALSE) %>% 
  # Remove -mRNA labels from gene names (so that it matches up with names in Heloderma Z list)
  mutate(Heloderma = str_remove(Heloderma, "-mRNA-.*"))

# Keep only genes that are in the Gila Z gametolog list
shinisaurus_ortho <- shinisaurus %>%
  semi_join(gila_z, by = c("Heloderma" = "Z")) %>% 
  distinct(Heloderma, Shinisaurus)

# Save list
write_tsv(shinisaurus_ortho, "outgroups/orthologs/shinisaurus_ortho.tsv")

# === 8) Thamnophis elegans ===
thamnophis <- read_tsv("../Heloderma_orthofinder/Heloderma_v_Thamnophis.tsv", show_col_types = FALSE) %>% 
  # Remove -mRNA labels from gene names (so that it matches up with names in Heloderma Z list)
  mutate(Heloderma = str_remove(Heloderma, "-mRNA-.*"))

# Keep only genes that are in the Gila Z gametolog list
thamnophis_ortho <- thamnophis %>%
  semi_join(gila_z, by = c("Heloderma" = "Z")) %>% 
  distinct(Heloderma, Thamnophis)

# Save list
write_tsv(thamnophis_ortho, "outgroups/orthologs/thamnophis_ortho.tsv")

# === 9) Varanus acanthurus ===
varanus <- read_tsv("../Heloderma_orthofinder/Heloderma_v_Varanus.tsv", show_col_types = FALSE) %>% 
  # Remove -mRNA labels from gene names (so that it matches up with names in Heloderma Z list)
  mutate(Heloderma = str_remove(Heloderma, "-mRNA-.*"),
         Varanus = str_remove(Varanus, "-mRNA-.*"))

# Keep only genes that are in the Gila Z gametolog list
varanus_ortho <- varanus %>%
  semi_join(gila_z, by = c("Heloderma" = "Z")) %>% 
  distinct(Heloderma, Varanus)

# Save list
write_tsv(varanus_ortho, "outgroups/orthologs/varanus_ortho.tsv")
