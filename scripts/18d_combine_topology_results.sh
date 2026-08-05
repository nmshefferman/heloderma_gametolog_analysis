#!/bin/bash

# This script combines all the tree topology resutls from each iqtree output into one joint table, with one row per gene.

cd /data/Wilson_Lab/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/iqtree/topology_tests

out="tree_topology_results.tsv"

printf "Gene\tTree\tlogL\tdeltaL\tbp-RELL\tp-KH\tp-SH\tp-WKH\tp-WSH\tc-ELW\tp-AU\tDataset\n" > "$out"

for dataset in all_species no_shinisaurus; do
  for f in "$dataset/stats/"*.iqtree; do
    gene=$(basename "$f" .iqtree)

    awk -v gene="$gene" -v dataset="$dataset" '
      /^USER TREES/ { in_user=1; next }
      in_user && /^[[:space:]]*Tree[[:space:]]+logL/ { next }
      in_user && /^[[:space:]]*-+[[:space:]]*$/ { next }

      in_user && /^[[:space:]]*[123][[:space:]]/ {
        # Fields for IQ-TREE USER TREES lines:
        # $1  tree
        # $2  logL
        # $3  deltaL
        # $4  bp-RELL
        # $6  p-KH
        # $8  p-SH
        # $10 p-WKH
        # $12 p-WSH
        # $14 c-ELW
        # $16 p-AU
        print gene "\t" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $6 "\t" $8 "\t" $10 "\t" $12 "\t" $14 "\t" $16 "\t" dataset
      }
    ' "$f" >> "$out"
  done
done