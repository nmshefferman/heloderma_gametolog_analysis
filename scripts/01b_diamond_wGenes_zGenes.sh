#!/bin/bash
#SBATCH --job-name=gametolog_blast_wGenes_ZGenes
#SBATCH -o /data/Wilson_Lab/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/scripts/gametolog_blast_wGenes_ZGenes.out
#SBATCH -e /data/Wilson_Lab/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/scripts/gametolog_blast_wGenes_ZGenes.err   
#SBATCH --mem=10g
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL

# This script uses DIAMOND to blast W-linked Heloderma suspectum genes against Z-linked genes in the genome to identify potential gametologs.

# Set working directory
cd /data/Wilson_Lab/squamates/Heloderma_sexchr/gametolog_analysis/

# Load in necessary modules
module load diamond

# Make any necessary directories
mkdir -p blast_results

# Make the DIAMOND database from the chrZ-only peptides (to compare Z genes against W genes)
## Need to make DIAMOND database in order to use it for blasting
diamond makedb --in Heloderma_suspectum.chrZ-only.peptides.fasta.gz -d Heloderma_suspectum.chrZ-only.peptides.dmnd

# Blast each W gene against all Z genes 
diamond blastp -d Heloderma_suspectum.chrZ-only.peptides.dmnd -q Heloderma_suspectum.chrW-only_PAR-masked.peptides.fasta.gz -o sheffermannm/heloderma_only/blast_results/wGenes_to_all_zGenes_all_isoforms.tsv --unal 1 --max-target-seqs 5 --header

# Collapse isoforms (only one row per gene pair)
awk 'BEGIN{FS=OFS="\t"}
{
    q=$1; s=$2;
    sub(/-mRNA-.*/, "", q);
    sub(/-mRNA-.*/, "", s);

    key=q OFS s;

    if (!(key in seen)) {
        seen[key]=1;
        print $0;
    }
}' sheffermannm/heloderma_only/blast_results/wGenes_to_all_zGenes_all_isoforms.tsv > sheffermannm/heloderma_only/blast_results/wGenes_zGenes.tsv

# Start: 12:59:27
# End: 12:59:33
# Total run time: 00:00:06