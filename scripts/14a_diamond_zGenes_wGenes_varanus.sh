#!/bin/bash
#SBATCH --job-name=gametolog_blast_zGenes_wGenes
#SBATCH -o varanus_gametolog_blast_zGenes_wGenes.out
#SBATCH -e varanus_gametolog_blast_zGenes_wGenes.err   
#SBATCH --mem=10g
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL

# This script uses DIAMOND to blast Z-linked Varanus acanthurus genes against W-linked genes in the genome to identify potential gametologs.

# Set working directory
cd /data/Wilson_Lab/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/varanus_only

# Load in necessary modules
module load diamond

# Make any necessary directories
mkdir -p blast_results

# Make the DIAMOND database from the chrW-only peptides (to compare W genes against Z genes)
## Need to make DIAMOND database in order to use it for blasting
## However, the `Varanus_acanthurus_chrW.peptides.fasta.gz` file contains unexpected characters that cause DIAMOND to fail. I first need to make a clean fasta file:
zcat ../references/Varanus_acanthurus_chrW.peptides.fasta.gz \
| awk '
BEGIN{RS=">"; ORS=""}
NR>1{
    n=split($0,a,"\n")
    hdr=a[1]

    seq=""
    for(i=2;i<=n;i++) seq=seq a[i]

    if(seq !~ /\./)
        print ">" hdr "\n" seq "\n"
}' > Varanus_acanthurus_chrW.peptides.clean.fasta

# Now make database
diamond makedb --in ../references/Varanus_acanthurus_chrW.peptides.clean.fasta -d ../references/Varanus_acanthurus_chrW.peptides.dmnd

# The `Varanus_acanthurus_chrZ.peptides.fasta.gz` file also contains unexpected characters, so I also need to make a clean fasta file:
zcat ../references/Varanus_acanthurus_chrZ.peptides.fasta.gz \
| awk '
BEGIN{RS=">"; ORS=""}
NR>1{
    n=split($0,a,"\n")
    hdr=a[1]

    seq=""
    for(i=2;i<=n;i++) seq=seq a[i]

    if(seq !~ /\./)
        print ">" hdr "\n" seq "\n"
}' > ../references/Varanus_acanthurus_chrZ.peptides.clean.fasta

# Blast each Z gene against all W genes
diamond blastp -d ../references/Varanus_acanthurus_chrW.peptides.dmnd -q ../references/Varanus_acanthurus_chrZ.peptides.clean.fasta -o blast_results/zGenes_to_all_wGenes_all_isoforms.tsv --unal 1 --max-target-seqs 5 --header

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
}' blast_results/zGenes_to_all_wGenes_all_isoforms.tsv > blast_results/zGenes_wGenes.tsv

# Start: 16:30:41 
# End: 16:30:45
# Total run time: 0:00:04