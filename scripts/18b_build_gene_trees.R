# ------
# Goal
#   Visualize gene trees built using IQ-TREE.
# ------

# Load in packages
library(tidyverse)
library(ape)
library(ggtree)

# Set working directory
setwd("/vf/users/Wilson_Lab/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/iqtree")

# Make a function to plot each gene tree
# === 1) All species ===
plot_tree_all <- function(treefile, outdir = "../figures/tree_plots/all_species") {
  dir.create(outdir)
  
  tree <- read.tree(treefile)
  
  # Root on Podarcis
  tree_rooted <- root(tree, outgroup = "Podarcis", resolve.root = TRUE)
  
  # Build title for each tree from filename
  title <- basename(treefile) %>% 
    str_remove("\\.treefile$") %>% 
    str_replace_all("__", " | ")
  
  p <- ggtree(tree_rooted) +
    geom_tiplab() +
    geom_text2(aes(subset = !isTip, label = label), hjust = -0.3, size = 3) +
    ggtitle(title)
  
  ggsave(filename = file.path(outdir, paste0(basename(treefile), ".pdf")), plot = p, width = 45, height = 6)
  
  return(p)
}

# Run the function on trees
tree_files_all <- list.files("make_trees/all_species", pattern = "\\.treefile$", full.names = TRUE)
plots <- lapply(tree_files_all, plot_tree_all)


# === 2) No Shinisaurus ===
plot_tree_subset <- function(treefile, outdir = "../figures/tree_plots/no_shinisaurus") {
  dir.create(outdir)
  
  tree <- read.tree(treefile)
  
  # Root on Podarcis
  tree_rooted <- root(tree, outgroup = "Podarcis", resolve.root = TRUE)
  
  # Build title for each tree from filename
  title <- basename(treefile) %>% 
    str_remove("\\.treefile$") %>% 
    str_replace_all("__", " | ")
  
  p <- ggtree(tree_rooted) +
    geom_tiplab() +
    geom_text2(aes(subset = !isTip, label = label), hjust = -0.3, size = 3) +
    ggtitle(title)
  
  ggsave(filename = file.path(outdir, paste0(basename(treefile), ".pdf")), plot = p, width = 45, height = 6)
  
  return(p)
}

# Next, run the function on trees that include all species except Shinisaurus (save to the same tree_plot directory)
tree_files_subset <- list.files("make_trees/no_shinisaurus", pattern = "\\.treefile$", full.names = TRUE)
plots_v2 <- lapply(tree_files_subset, plot_tree_subset)