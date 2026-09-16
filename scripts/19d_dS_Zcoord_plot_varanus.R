# ------
# Goal
#   Plot dS vs. chrZ position (and chrW position) for Varanus gametolog genes.
# Contents
# 1) Make table for plotting
# 2) Plot dS v. chrZ position
#   a) dS vs. chrZ -- Normal
# 3) Plot dS vs. chrW position
#   a) dS vs. chrW -- Normal
# 4) Plot dS v. chrZ position USING GILA CHR Z COORDINATES
# ------

# Load in packages
library(tidyverse)

# Set working directory
setwd("/vf/users/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/varanus_only")

# === 1) Make table for plotting ===
# Read in final gametolog table (has dS values)
final_gametologs <- read_tsv("varanus_Z_W_gametologs_final.tsv", show_col_types = FALSE)

# Add some modifications
plot_tbl_prank <- final_gametologs %>% 
  # Write position in megabases
  mutate(Z_start = Z_start/1e6, Z_end = Z_end/1e6, W_start = W_start/1e6, W_end = W_end/1e6)

# === 2) Plot dS v. chrZ position ===
# ==== a) dS vs. chrZ -- Normal ====
# Plot
plot_prank_z <- ggplot(plot_tbl_prank, aes(x = Z_start, y = dS)) +
  geom_point() +
  coord_cartesian(xlim = c(0, 11.9), ylim = c(0, 1.5)) +
  # Force X axis to start at 0 with no padding buffer
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  # Force Y axis to start at 0 with no padding buffer
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Position on chrZ (Mb)", y = "dS", title = "dS vs Z-chromosome position (Varanus acanthurus)")
plot_prank_z

# Save plot
ggsave("../figures/dS_Zcoord_plot_varanus.png", plot = plot_prank_z, dpi = 300, width = 10, height = 6)
ggsave("../figures/dS_Zcoord_plot_varanus.pdf", plot = plot_prank_z, dpi = 300, width = 10, height = 6)

# ==== 3) dS vs. chrW position ====
# ==== a) dS vs. chrW -- Normal ====
plot_prank_w <- ggplot(plot_tbl_prank, aes(x = W_start, y = dS)) +
  geom_point() +
  coord_cartesian(xlim = c(0, 15.1), ylim = c(0, 1.5)) +
  # Force X axis to start at 0 with no padding buffer
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  # Force Y axis to start at 0 with no padding buffer
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Position on chrW (Mb)", y = "dS", title = "dS vs W-chromosome position (Varanus acanthurus)")
plot_prank_w

# Save plot
ggsave("../figures/dS_Wcoord_plot_varanus.png", plot = plot_prank_w, dpi = 300, width = 10, height = 6)
ggsave("../figures/dS_Wcoord_plot_varanus.pdf", plot = plot_prank_w, dpi = 300, width = 10, height = 6)

# === 4) Plot dS v. chrZ position USING GILA CHR Z COORDINATES ===
# Read Heloderma-Varanus ortholog table
orthologs <- read_tsv("../iqtree/combined_table.tsv", show_col_types = FALSE) %>% 
  select(Heloderma_Z, Varanus_Z) %>% 
  mutate(Heloderma_Z = str_remove(Heloderma_Z, "-mRNA-.*"),
         Varanus_Z = str_remove(Varanus_Z, "-mRNA-.*"),
         # Standardize Varanus Z names between the two tables
         Varanus_Z = str_replace(Varanus_Z, "LOC_", "LOCZ_")) %>% 
  distinct()

# Read in Gila coordinates
gila_coords <- read_tsv("../heloderma_only/heloderma_Z_W_gametologs_final.tsv", show_col_types = FALSE) %>% 
  select(Heloderma_Z = Z, Gila_Z_start = Z_start, Gila_Z_end = Z_end)

# Combine tables
plot_tbl_gila_coords <- final_gametologs %>% 
  # Match Varanus Z to its Heloderma Z ortholog
  left_join(orthologs, by = c("Z" = "Varanus_Z")) %>% 
  # Add the Heloderma chrZ coordinates
  left_join(gila_coords, by = "Heloderma_Z") %>% 
  # Convert Heloderma position to Mb
  mutate(Gila_Z_start = Gila_Z_start/1e6)

# Plot Varanus dS against Heloderma chrZ position
plot_varanus_gila_z <- ggplot(plot_tbl_gila_coords, aes(x = Gila_Z_start, y = dS)) +
  geom_point() +
  coord_cartesian(xlim = c(0, 20.6), ylim = c(0, 1.5)) +
  # Force X axis to start at 0 with no padding buffer
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  # Force Y axis to start at 0 with no padding buffer
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Position on Heloderma chrZ (Mb)", y = "Varanus dS", title = "Varanus dS vs Heloderma Z-chromosome position")
plot_varanus_gila_z

# Save plot
ggsave("../figures/dS_Zcoord_plot_varanus_gilaZ.png", plot = plot_varanus_gila_z, dpi = 300, width = 10, height = 6)
ggsave("../figures/dS_Zcoord_plot_varanus_gilaZ.pdf", plot = plot_varanus_gila_z, dpi = 300, width = 10, height = 6)
