######
# Goal
#   Make a plot mapping the Z-W gametolog positions along each chromosome. Lines connect gametolog pairs and colors represent the relative distance from the start of the chromosome.
######

# Load in packages
library(tidyverse)

# Set working directory
setwd("/vf/users/Wilson_Lab/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only")

# Read in gametolog table
gametologs <- read_tsv("Z_W_gametologs.tsv", show_col_types = FALSE) %>%
  # Keep track of gametolog pairs by labeling row number
  mutate(pair_id = row_number())

# Build a long table with one point per gene
plot_points <- bind_rows(gametologs %>% 
                           mutate(chrom = "Z", gene = Z, pos = Z_start) %>% 
                           select(pair_id, chrom, gene, pos),
                         gametologs %>% 
                           mutate(chrom = "W", gene = W, pos = W_start) %>% 
                           select(pair_id, chrom, gene, pos)) %>% 
  group_by(chrom) %>% 
  # Calculate distance (as a percentage) for coloring
  mutate(pos_scaled = pos/max(pos, na.rm = TRUE)) %>% 
  ungroup() %>% 
  # Label Z genes as 1 and W genes as 2 (for x-axis)
  mutate(x = if_else(chrom == "Z", 1, 2))

# Build line segments connecting each Z-W pair
plot_segments <- gametologs %>% 
  # Plot y-axis in terms of megabases
  mutate(x = 1, xend = 2, y = Z_start/1e6, yend = W_start/1e6) %>% 
  select(pair_id, x, xend, y, yend)

# Plot 
plot <- ggplot() +
  geom_segment(data = plot_segments,
               aes(x = x, xend = xend, y = y, yend = yend, group = pair_id),
               color = "grey60",
               alpha = 0.35,
               linewidth = 0.4) +
  geom_point(data = plot_points,
             aes(x = x, y = pos/1e6, color = pos_scaled),
             size = 2) +
  scale_x_continuous(breaks = c(1, 2),
                     labels = c("Z", "W")) +
  scale_color_gradient(low = "lightblue",
                       high = "navy") +
  labs(x = "Chromosome",
       y = "Chromosome position (Mb)",
       title = "H. suspectum Z and W Gametolog Positions",
       color = "Relative position") +
  theme_classic() #+
  #theme(legend.position = "none")
plot

# Save plot
ggsave("../figures/heloderma_Z_W_gametolog_positions.png", plot = plot, dpi = 300, width = 10, height = 6)