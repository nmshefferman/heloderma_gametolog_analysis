# ------
# Goal
#   Plot dS vs. chrZ position (and chrW position) for Heloderma gametolog genes.
# Contents
# 1) Make table for plotting
# 2) Plot dS v. chrZ position
# ------

# Load in packages
library(tidyverse)

# Set working directory
setwd("/vf/users/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/PAR")

# === 1) Make table for plotting ===
# Read in final PAR table (does not have dS values so we have to merge the dS table before proceeding)
PAR <- read_tsv("Z_W_PARs.tsv", show_col_types = FALSE) %>%
  mutate(W = str_remove(W, "-mRNA-.*"),
         Z = str_remove(Z, "-mRNA-.*"))

# Read in dS value table for PAR genes
PAR_dS <- read_tsv("paml/paml_yn00_summary.tsv", show_col_types = FALSE) %>%
  separate(pair, into = c("W", "Z"), sep = "__")

# Add dS values to PAR genes
PAR <- PAR %>%
  left_join(PAR_dS %>% select(W, Z, dN, dS, omega), by = c("W", "Z")) %>%
  mutate(region = "PAR")

# Read in final gametolog table (has dS values)
final_gametologs <- read_tsv("../heloderma_Z_W_gametologs_final.tsv", show_col_types = FALSE) %>% 
  mutate(region = "non-PAR")

# Merge both tables
plot_tbl <- bind_rows(
  final_gametologs %>%
    select(W, W_start, W_end, Z, Z_start, Z_end, dS) %>%
    mutate(region = "non-PAR"),
  PAR %>%
    select(W, W_start, W_end, Z, Z_start, Z_end, dS) %>%
    mutate(region = "PAR"))

# Convert coordinates to Mb
plot_tbl <- plot_tbl %>%
  mutate(Z_start = Z_start / 1e6,
         Z_end = Z_end / 1e6,
         W_start = W_start / 1e6,
         W_end = W_end / 1e6)

# === 2) Plot dS v. chrZ position ===
# Plot
plot_dS_PAR <- ggplot(plot_tbl, aes(x = Z_start, y = dS, color = region)) +
  geom_point() +
  coord_cartesian(xlim = c(0,20.6), ylim = c(0, 0.45)) +
  # Force X axis to start at 0 with no padding buffer
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  # Force Y axis to start at 0 with no padding buffer
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  geom_vline(xintercept = 2209992 / 1e6, linetype = "dotted", color = "red", linewidth = 1) +
  scale_color_manual(values = c("non-PAR" = "black", "PAR" = "red")) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(
    x = "Position on chrZ (Mb)",
    y = "dS",
    title = "dS vs Z-chromosome position (Heloderma suspectum)")
plot_dS_PAR

# Save plot
ggsave("../../figures/dS_Zcoord_plot_gila_PAR.png", plot = plot_dS_PAR, dpi = 300, width = 10, height = 6)
ggsave("../../figures/dS_Zcoord_plot_gila_PAR.pdf", plot = plot_dS_PAR, dpi = 300, width = 10, height = 6)

