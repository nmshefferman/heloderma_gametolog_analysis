#!/bin/bash
#SBATCH --job-name=gametolog_blast_wGenes_allGenes
#SBATCH -o heloderma_gametolog_blast_wGenes_allGenes.out
#SBATCH -e heloderma_gametolog_blast_wGenes_allGenes.err   
#SBATCH --mem=10g
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL

# This script uses DIAMOND to blast W-linked Heloderma suspectum genes against all genes in the genome to identify potential gametologs.

# Set working directory
cd /data/Wilson_Lab/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only

# Load in necessary modules
module load diamond

# Make any necessary directories
mkdir -p blast_results

# Make the DIAMOND database from the whole genome chrW masked peptides (to compare W genes against all genes)
## Need to make DIAMOND database in order to use it for blasting
diamond makedb --in ../references/Heloderma_suspectum.final.reference.chrW-hardmasked.peptides.fasta.gz -d ../references/Heloderma_suspectum.final.reference.chrW-hardmasked.peptides.dmnd
    
# Blast each W protein coding gene against all genes in the genome
diamond blastp -d ../references/Heloderma_suspectum.final.reference.chrW-hardmasked.peptides.dmnd -q ../references/Heloderma_suspectum.chrW-only_PAR-masked.peptides.fasta.gz -o blast_results/wGenes_to_allGenes_all_isoforms.tsv --unal 1 --max-target-seqs 5 --header 

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
}' blast_results/wGenes_to_allGenes_all_isoforms.tsv > blast_results/wGenes_allGenes.tsv

# Start: 12:59:27
# End: 12:59:42
# Total run time: 00:00:15