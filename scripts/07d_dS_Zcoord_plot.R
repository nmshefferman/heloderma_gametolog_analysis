# ------
# Goal
#   Plot dS vs. chrZ position (and chrW position) for Heloderma gametolog genes.
# Contents
# 1) Make table for plotting
# 2) Plot dS v. chrZ position
#   a) dS vs. chrZ -- Normal
#   b) dS vs. chrZ -- Colored
#   c) dS vs. chrZ -- Highlighting the 35 shared orthologs along the Z chrom
# 3) Plot dS vs. chrW position
#   a) dS vs. chrW -- Normal
#   b) dS vs. chrW -- Colored
# ------

# Load in packages
library(tidyverse)

# Set working directory
setwd("/vf/users/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only")

# === 1) Make table for plotting ===
# Read in final gametolog table (has dS values)
final_gametologs <- read_tsv("heloderma_Z_W_gametologs_final.tsv", show_col_types = FALSE)

# Read in list of genes that are orthologs (after curating alignments)
final_orthologs <- c("LOC_00019078", "LOC_00019083", "LOC_00019092", "LOC_00019100", "LOC_00019111", "LOC_00019114", "LOC_00019116", 
                     "LOC_00019120", "LOC_00019130", "LOC_00019157", "LOC_00019162", "LOC_00019168", "LOC_00019169", "LOC_00019183",
                     "LOC_00019197", "LOC_00019204", "LOC_00019207", "LOC_00019236", "LOC_00019256", "LOC_00019264", "LOC_00019282",
                     "LOC_00019284", "LOC_00019285", "LOC_00019299", "LOC_00019305", "LOC_00019325", "LOC_00019334")

# Add some modifications
plot_tbl_prank <- final_gametologs %>% 
  # Write position in megabases
  mutate(Z_start = Z_start/1e6, Z_end = Z_end/1e6, W_start = W_start/1e6, W_end = W_end/1e6) %>% 
  # Separate graph into 5 "chunks" and color them 
  mutate(region = case_when(Z_start < 5 ~ "0-5",
                            Z_start < 10 ~ "5-10",
                            Z_start < 15 ~ "10-15",
                            Z_start < 19 ~ "15-19",
                            TRUE   ~ "19+")) %>% 
  mutate(ortholog = if_else(Z %in% final_orthologs, "yes", "no"))

# === 2) Plot dS v. chrZ position ===
# ==== a) dS vs. chrZ -- Normal ====
# Plot
plot_prank_z <- ggplot(plot_tbl_prank, aes(x = Z_start, y = dS)) +
  geom_point() +
  coord_cartesian(xlim = c(0,20.6), ylim = c(0, 0.45)) +
  # Force X axis to start at 0 with no padding buffer
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  # Force Y axis to start at 0 with no padding buffer
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Position on chrZ (Mb)", y = "dS", title = "dS vs Z-chromosome position (Heloderma suspectum)") +
  # Add vertical line to indicate PAR location
  geom_vline(xintercept = 2209992/1e6, linetype = "dotted", color = "red", linewidth = 1)
plot_prank_z

# Save plot
ggsave("../figures/dS_Zcoord_plot_gila.png", plot = plot_prank_z, dpi = 300, width = 10, height = 6)
ggsave("../figures/dS_Zcoord_plot_gila.pdf", plot = plot_prank_z, dpi = 300, width = 10, height = 6)


# ==== b) dS vs. chrZ -- Colored ====
plot_prank_z_colored <- ggplot(plot_tbl_prank, aes(x = Z_start, y = dS, color = region)) +
  geom_point() +
  coord_cartesian(xlim = c(0,20.6), ylim = c(0, 0.45)) +
  # Force X axis to start at 0 with no padding buffer
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  # Force Y axis to start at 0 with no padding buffer
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  scale_color_manual(values = c(
    "0-5"   = "grey60",
    "5-10"  = "#E69F00",
    "10-15" = "#009E73",
    "15-19" = "#CC79A7",
    "19+"   = "#0072B2")) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Position on chrZ (Mb)", y = "dS", title = "dS vs Z-chromosome position (Heloderma suspectum)")
plot_prank_z_colored

# Save plot
ggsave("../figures/dS_Zcoord_plot_gila_colored.png", plot = plot_prank_z_colored, dpi = 300, width = 10, height = 6)
ggsave("../figures/dS_Zcoord_plot_gila_colored.pdf", plot = plot_prank_z_colored, dpi = 300, width = 10, height = 6)

# ==== c) dS vs. chrZ -- Highlighting the 35 shared orthologs along the Z chrom ====
#plot_prank_z_orthologs <- ggplot(plot_tbl_prank, aes(x = Z_start, y = dS, color = ortholog)) +
#  geom_point() +
#  coord_cartesian(xlim = c(0, 20.6), ylim = c(0, 0.45)) +
  # Force X axis to start at 0 with no padding buffer
#  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  # Force Y axis to start at 0 with no padding buffer
#  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
#  scale_color_manual(values = c(
#    "no"   = "grey60",
#    "yes"  = "#E69F00")) +
#  theme_classic() +
#  theme(legend.position = "none") +
#  labs(x = "Position on chrZ (Mb)", y = "dS", title = "dS vs Z-chromosome position (H. suspectum) (Colors represent genes that have orthologs in other species)")
#plot_prank_z_orthologs

# Save plot
#ggsave("../figures/dS_Zcoord_plot_gila_orthologs.png", plot = plot_prank_z_orthologs, dpi = 300, width = 12, height = 6)
#ggsave("../figures/dS_Zcoord_plot_gila_orthologs.pdf", plot = plot_prank_z_orthologs, dpi = 300, width = 12, height = 6)

# ==== 3) dS vs. chrW position ====
# ==== a) dS vs. chrW -- Normal ====
plot_prank_w <- ggplot(plot_tbl_prank, aes(x = W_start, y = dS)) +
  geom_point() +
  coord_cartesian(xlim = c(0, 70), ylim = c(0, 0.45)) +
  # Force X axis to start at 0 with no padding buffer
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  # Force Y axis to start at 0 with no padding buffer
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Position on chrW (Mb)", y = "dS", title = "dS vs W-chromosome position (Heloderma suspectum)") +
  # Add vertical line to indicate PAR location
  geom_vline(xintercept = 2205462/1e6, linetype = "dotted", color = "red", linewidth = 1)
plot_prank_w

# Save plot
ggsave("../figures/dS_Wcoord_plot_gila.png", plot = plot_prank_w, dpi = 300, width = 10, height = 6)
ggsave("../figures/dS_Wcoord_plot_gila.pdf", plot = plot_prank_w, dpi = 300, width = 10, height = 6)

# ==== b) dS vs. chrW -- Colored ====
plot_prank_w_colored <- ggplot(plot_tbl_prank, aes(x = W_start, y = dS, color = region)) +
  geom_point() +
  coord_cartesian(xlim = c(0, 70), ylim = c(0, 0.45)) +
  # Force X axis to start at 0 with no padding buffer
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  # Force Y axis to start at 0 with no padding buffer
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  scale_color_manual(values = c(
    "0-5"   = "grey60",
    "5-10"  = "#E69F00",
    "10-15" = "#009E73",
    "15-19" = "#CC79A7",
    "19+"   = "#0072B2")) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Position on chrW (Mb)", y = "dS", title = "dS vs W-chromosome position (Heloderma suspectum)")
plot_prank_w_colored

# Save plot
ggsave("../figures/dS_Wcoord_plot_gila_colored.png", plot = plot_prank_w_colored, dpi = 300, width = 10, height = 6)
ggsave("../figures/dS_Wcoord_plot_gila_colored.pdf", plot = plot_prank_w_colored, dpi = 300, width = 10, height = 6)
