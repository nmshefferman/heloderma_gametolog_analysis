# Evolution of Z-W gametologs in the Gila monster

### Goal: Look at gametolog evolution in the Gila monster
#### What is a gametolog?
Pairs of homologous genes that have retained functional copies on both sex chromosomes, even if they are present in the non-recombining region of the chromosomes. These genes are derived from a common ancestral gene but have diverged because the chromosomes stopped recombining some time ago. 

## *Heloderma suspectum* Z v. W
#### /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only
The first step in this project is to compare gametologs on the *Heloderma suspectum* sex chromosomes.

**`scenario_expectation_blast_results.txt` file explains different expectations of how gene pairs might appear in the BLAST tables depending on how recently they stopped recombining.**

### 1) Identify Z-W gametologs via BLAST database search 
#### a) Z to W 
##### i) BLAST each Z gene to all W genes (using diamond)
https://hpc.nih.gov/apps/diamond.html 
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only
sbatch ../scripts/01a_diamond_zGenes_wGenes.sh
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/blast_results/zGenes_wGenes.tsv`

##### ii) Make a list of Z gene IDs from BLAST output table 
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/blast_results
awk 'BEGIN{FS="\t"} !/^#/ && $1 != "*" {print $1}' zGenes_wGenes.tsv | sort -u > ../gene_lists/Z-genes_transcripts.txt
# Convert trasncript IDs to gene IDs (remove -mRNA flag)
sed 's/-mRNA-.*//' ../gene_lists/Z-genes_transcripts.txt | sort -u > ../gene_lists/Z-genes.txt
rm ../gene_lists/Z-genes_transcripts.txt
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/gene_lists/Z-genes.txt`

#### b) W to Z
##### i) BLAST each W gene to all Z genes (using diamond)
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only
sbatch ../scripts/01b_diamond_wGenes_zGenes.sh
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/blast_results/wGenes_zGenes.tsv`

##### ii) Make a list of W gene IDs from BLAST output table
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/blast_results
awk 'BEGIN{FS="\t"} !/^#/ && $2 != "*" {print $1}' wGenes_zGenes.tsv | sort -u > ../gene_lists/W-genes_transcripts.txt
# Convert trasncript IDs to gene IDs (remove -mRNA flag)
sed 's/-mRNA-.*//' ../gene_lists/W-genes_transcripts.txt | sort -u > ../gene_lists/W-genes.txt
rm ../gene_lists/W-genes_transcripts.txt
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/gene_lists/W-genes.txt`

#### c) W to autosomes
##### i) BLAST each W protein coding gene against all genes in the genome
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/scripts
sbatch 01c_diamond_wGenes_allGenes.sh
```

! NOTE: The *Heloderma* reference genome had the W masked, but not the Z. This means that W genes will get hits to the Z, which we don't care about because we already have a table of W-Z and Z-W hits. We want to only compare W genes to the autosomes. 

This means that we have to filter out the Z genes from the BLAST result:
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/blast_results
awk 'BEGIN{FS=OFS="\t"}
NR==FNR { z[$1]=1; next }
!/^#/ && $2 != "*" {
    subj=$2
    sub(/-mRNA-.*/, "", subj)
    if (!(subj in z)) print
}' ../gene_lists/Z-genes.txt wGenes_allGenes.tsv > wGenes_autosomes.tsv
rm wGenes_allGenes.tsv

# Add this header to wGenes_autosomes.tsv:
# DIAMOND v2.1.12. http://github.com/bbuchfink/diamond
# Invocation: diamond blastp -d Heloderma_suspectum.final.reference.chrW-hardmasked.peptides.dmnd -q Heloderma_suspectum.chrW-only_PAR-masked.peptides.fasta.gz -o sheffermannm/blast_results/wGenes_to_allGenes_all_isoforms.tsv --unal 1 --max-target-seqs 1 --header
# Fields: Query ID (W), Subject ID (autosomes), Percentage of identical matches, Alignment length, Number of mismatches, Number of gap openings, Start of alignment in query, End of alignment in query, Start of alignment in subject, End of alignment in subject, Expected value, Bit score
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/blast_results/heloderma_only/wGenes_autosomes.tsv`

##### ii) Make a list of autosome gene IDs from BLAST output table
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/blast_results
awk 'BEGIN{FS="\t"} !/^#/ && $2 != "*" {print $2}' wGenes_autosomes.tsv | sort -u > ../gene_lists/autosome-genes_transcripts.txt
# Convert trasncript IDs to gene IDs (remove -mRNA flag)
sed 's/-mRNA-.*//' ../gene_lists/autosome-genes_transcripts.txt | sort -u > ../gene_lists/autosome-genes.txt
rm ../gene_lists/autosome-genes_transcripts.txt
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/gene_lists/autosome-genes.txt`

### 2) Find gene coordinates
Z-genes
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/gene_coords

zgrep -F -f /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/gene_lists/Z-genes.txt /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/references/Heloderma_suspectum.final.reference.gff3.gz | awk '$3=="gene"' > Z_gametolog_coords.gff3

# Add this header: #seqid  source  type  start  end  score  strand  phase  attributes
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/gene_coords/Z_gametolog_coords.gff3`

W-genes
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/gene_coords

zgrep -F -f /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/gene_lists/W-genes.txt /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/references/Heloderma_suspectum.final.reference.gff3.gz | awk '$3=="gene"' > W_gametolog_coords.gff3

# Add this header: #seqid  source  type  start  end  score  strand  phase  attributes
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/gene_coords/W_gametolog_coords.gff3`

Autosomes 
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/gene_coords
zgrep -F -f /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/gene_lists/autosome-genes.txt /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/references/Heloderma_suspectum.final.reference.gff3.gz | awk '$3=="gene"' > autosome_gametolog_coords.gff3
# Add this header: #seqid  source  type  start  end  score  strand  phase  attributes
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/gene_coords/autosome_gametolog_coords.gff3`

### 3) Compare Z, W, and autosome genes that returned hits to get a list of Z-W gametologs
Make tables of Z-W-autosome best hits. Run R script:
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/scripts
module load R
Rscript 03_gametologs_best_hits.R
```
**Output:** The first table contains all Z genes that had a hit to a W gene. It also lists if the Z gene had a second hit to another W gene. For each gene, it lists the evalue and alignment length. The table also lists if the Z-W hit is also the strongest W-Z hit (TRUE/FALSE).
            `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/Z_W_summary.tsv`,
            The second table contains all W genes that had a hit to a Z gene. It also lists if the W gene had a second hit to another Z gene and/or if the W gene had a hit to an autosome. For each gene, it lists the evalue and alignment length. The table also lists whether the Z or autosome gene is the best match for the W gene. If the Z gene is the best match, then it lists what the best W-Z hit is. If the best W-Z hit is the same as the W gene in the first column, then they are reciprocal best hits (TRUE/FALSE).
            `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/W_Z_autosome_summary.tsv`, 
            The third table contains only Z and W genes that were reciprocal best hits. These are the most high-confidence gametologs before aligning. The table contains the name and start and end coordinates for each gene.
            `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/Z_W_gametologs.tsv`

! NOTE: The coordinates of the chrZ PAR are 0-2209992 and the coordinate of the chrW PAR are 0-2205462.

### 4) Extract Z and W gametologs from FASTA files
Z genes
```
# Extract Gila Z genes from gametolog table
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only

tail -n +2 Z_W_gametologs.tsv | cut -f5 | sed 's/-mRNA-.*//' | sort -u > gene_lists/gila_gametolog_Z_genes.txt

# Using the list of Gila Z gametologs, extract peptide FASTA files
module load seqkit
seqkit grep -n -r -f gene_lists/gila_gametolog_Z_genes.txt ../references/Heloderma_suspectum.chrZ-only.peptides.fasta.gz -o gene_lists/Z_gametologs_isoforms.fasta

# Keep only one isoform per gene (the longest one)
awk '
BEGIN { RS=">"; ORS="" }
NR > 1 {
    n = split($0, a, "\n")
    hdr = a[1]
    id = hdr
    sub(/ .*/, "", id)

    gene = id
    sub(/-mRNA-.*/, "", gene)

    seq = ""
    for (i = 2; i <= n; i++) seq = seq a[i]

    if (!(gene in maxlen) || length(seq) > maxlen[gene]) {
        maxlen[gene] = length(seq)
        best_hdr[gene] = hdr
        best_seq[gene] = seq
    }
}
END {
    for (g in best_seq)
        print ">" best_hdr[g] "\n" best_seq[g] "\n"
}
' gene_lists/Z_gametologs_isoforms.fasta > gene_lists/Z_gametologs.fasta
rm gene_lists/Z_gametologs_isoforms.fasta
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/gene_lists/Z_gametologs.fasta`

W genes
```
# Extract Gila W genes from gametolog table
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only

tail -n +2 Z_W_gametologs.tsv | cut -f1 | sed 's/-mRNA-.*//' | sort -u > gene_lists/gila_gametolog_W_genes.txt

# Using the list of Gila W gametologs, extract peptide FASTA files
module load seqkit
seqkit grep -n -r -f gene_lists/gila_gametolog_W_genes.txt ../references/Heloderma_suspectum.chrW-only_PAR-masked.peptides.fasta.gz -o gene_lists/W_gametologs_isoforms.fasta

# Keep only one isoform per gene
awk '
BEGIN { RS=">"; ORS="" }
NR > 1 {
    n = split($0, a, "\n")
    hdr = a[1]
    id = hdr
    sub(/ .*/, "", id)

    gene = id
    sub(/-mRNA-.*/, "", gene)

    seq = ""
    for (i = 2; i <= n; i++) seq = seq a[i]

    if (!(gene in maxlen) || length(seq) > maxlen[gene]) {
        maxlen[gene] = length(seq)
        best_hdr[gene] = hdr
        best_seq[gene] = seq
    }
}
END {
    for (g in best_seq)
        print ">" best_hdr[g] "\n" best_seq[g] "\n"
}
' gene_lists/W_gametologs_isoforms.fasta > gene_lists/W_gametologs.fasta
rm gene_lists/W_gametologs_isoforms.fasta
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/gene_lists/W_gametologs.fasta`

### 5) Plot where the genes are along the Z and W (sanity check step)
Again, this is easiest to do in R. Run R script:
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only
module load R
Rscript ../scripts/05_plot_ZW_genes.R
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/figures/heloderma_Z_W_gametolog_positions.png`,
            `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/figures/heloderma_Z_W_gametolog_positions.pdf`

### 6) Align gametolog partners
#### a) Make pairwise FASTA files for alignment (using gametolog table, make a file for each row)
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only
chmod u+x ../scripts/06a_pairwise_fastas_gila.sh
./../scripts/06a_pairwise_fastas_gila.sh
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/align/pairwise_fastas`

#### b) Align
Using PRANK (codon aware) https://hpc.nih.gov/apps/PRANK.html
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/align/
mkdir -p prank

# Align nucleotide gametolog pairs with PRANK
for f in pairwise_fastas/*.fasta; do
  base=$(basename "${f%.fasta}")
  # It is important to have the -codon flag so that the aligner is codon aware
  echo "prank -d=\"$f\" -o=prank/${base}.prank -codon -F"
done > prank/prank.swarm

# Run swarm
module load prank
swarm -f prank/prank.swarm
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/align/prank`

This directory contains an alignment file for each gametolog pair, generated using PRANK.

! NOTE: Check how many files are in the output directory before proceeding, should be the same number of gametolog pairs.

! NOTE: Using the codon-aware version of PRANK removes stop codons, so none of the aligned files have stop codons

#### c) Curate alignments
After running the alignments through PRANK, Brendan curated the alignments to exclude any with large gaps or other problems. This resulted in 100/106 alignments kept and 6/106 discarded. Of the 100 kept alignments, he curated 15 of them. I want to add all the best alignments to the same directory to make downstream anaylsis more straightforward:
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/align/

# Make new directory to keep all the best alignments
mkdir -p Heloderma_pairwise_gametologs_best

# Copy all original alignments to this new directory
cp prank/*.prank.best.fas Heloderma_pairwise_gametologs_best/

# Copy the curated alignments on top of the original prank files
for f in Heloderma_pairwise_gametologs_curated/*.fasta
do
    base=$(basename "$f")
    newname="${base/ (modified)/}"
    newname="${newname%.fasta}.fas"

    cp "$f" "Heloderma_pairwise_gametologs_best/$newname"
done

# Remove the 6 discarded files
cd Heloderma_pairwise_gametologs_best
rm LOCW_00021332* LOCW_00021360* LOCW_00021302* LOCW_00021308* LOCW_00021437* LOCW_00021304*
```
The final directory of high-confidence alignements can be found at `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/align/Heloderma_pairwise_gametologs_best`.

### 7) Compute dS (of non-PAR genes) and plot 
#### a) PAML (https://hpc.nih.gov/apps/PAML.html)
Use PAML to detect natural selection of proteins. Specifically, the yn00 feature of PAML calculates pairwise comparisons between two sequences and outputs dN/dS/ω values.
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/

mkdir -p paml/jobs paml/out
module load paml

# Update user permissions
chmod +x ../scripts/07a_run_yn00.sh

# Make into swarm file
# Find all PRANK codon alignments
find align/Heloderma_pairwise_gametologs_best -name '*.best.fas' | sort \
  | awk '{print "./../scripts/07a_run_yn00.sh \"" $0 "\""}' > paml/paml.swarm

# Run swarm
swarm -f paml/paml.swarm
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/paml/out`

This directory contains individual txt files of yn00 outputs for each gametolog pair, made using the PRANK aligned FASTAs.

! NOTE: The Yang & Nielsen (2000) section contains the reported dN, dS, omega (dN/dS ratio) values. The other methods printed in the txt output files are ignored.

#### b) Combine each pair dN, dS, and omega values into one summary table.
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/paml/

out="paml_yn00_summary.tsv"

printf "pair\tdN\tdS\tomega\n" > "$out"

for f in out/*.yn00.txt; do
  pair=$(basename "$f" .yn00.txt)

  awk -v pair="$pair" '
    /Yang & Nielsen \(2000\)/ {flag=1; next}
    flag && $1 ~ /^[0-9]+$/ && $2 ~ /^[0-9]+$/ && $7 ~ /^[0-9.]+$/ {
      # Fields for the YN00 line:
      # $7 = omega
      # $8 = dN
      # $11 = dS
      print pair "\t" $8 "\t" $11 "\t" $7
      exit
    }
  ' "$f" >> "$out"
done
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/paml/paml_yn00_summary.tsv`

#### c) Finalize *Heloderma suspectum* Z-W gametolog list
We also have to remove the 6 discarded genes from the final gametolog list and add the pairwise dS values.
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only

awk 'BEGIN{FS=OFS="\t"}
# Reads paml summary table and takes dS column (third column)
NR==FNR {
    if (FNR > 1) {
        dS[$1] = $3
    }
    next
}
NR==1 {
    print $0, "dS"
    next
}
$1!="LOCW_00021332" &&
$1!="LOCW_00021360" &&
$1!="LOCW_00021302" &&
$1!="LOCW_00021437" &&
$1!="LOCW_00021308" &&
$1!="LOCW_00021304" {
    pair = $1 "__" $5
    print $0, (pair in dS ? dS[pair] : "dS")
}' paml/paml_yn00_summary.tsv Z_W_gametologs.tsv > heloderma_Z_W_gametologs_final.tsv
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/heloderma_Z_W_gametologs_final.tsv`

#### d) Plot dS vs. position on the Z 
Using the final gametolog table, plot dS vs. position on the Z chromosome in R. Run R script:
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/scripts
module load R
Rscript 07d_dS_Zcoord_plot.R
```
**Output:** Plot dS v. chrZ postion (Mb) with a red dotted line showing the chrZ PAR boundary
            `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/figures/dS_Zcoord_plot_gila.png`,
            `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/figures/dS_Zcoord_plot_gila.pdf`,
            Plot dS v. chrZ position (Mb) with colors for every ~5 Mb (can be used to compare gene position on chrW)
            `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/figures/dS_Zcoord_plot_gila_colored.png`,
            `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/figures/dS_Zcoord_plot_gila_colored.pdf`,
            Plot dS v. chrW position (Mb) with a red dotted line showing the chrW PAR boundary
            `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/figures/dS_Wcoord_plot_gila.png`,
            `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/figures/dS_Wcoord_plot_gila.pdf`,
            Plot dS v. chrW position (Mb) with colors that correspond to chrZ gametologs (can be used to compare gene position on chrZ)
            `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/figures/dS_Wcoord_plot_gila_colored.png`
            `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/figures/dS_Wcoord_plot_gila_colored.pdf`

## *Heloderma suspectum* Z v. W -- PAR
Repeat analysis, this time with PAR genes only to calculate their dS value. 

### 8) Identify PAR homologs via BLAST
#### a) BLAST all Z genes to just W PAR genes
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only
sbatch ../scripts/08a_diamond_zGenes_wPAR.sh
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/PAR/blast_results/PAR_zGenes_wGenes.tsv`

##### i) Make a list of Z gene IDs from BLAST output table 
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/PAR/blast_results
awk 'BEGIN{FS="\t"} !/^#/ && $1 != "*" {print $1}' PAR_zGenes_wGenes.tsv | sort -u > ../gene_lists/PAR-Z-genes_transcripts.txt
# Convert trasncript IDs to gene IDs (remove -mRNA flag)
sed 's/-mRNA-.*//' ../gene_lists/PAR-Z-genes_transcripts.txt | sort -u > ../gene_lists/PAR-Z-genes.txt
rm ../gene_lists/PAR-Z-genes_transcripts.txt
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/PAR/gene_lists/PAR-Z-genes.txt`

#### b) BLAST just W PAR genes to all Z genes 
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only
sbatch ../scripts/08b_diamond_wPAR_zGenes.sh
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/PAR/blast_results/PAR_wGenes_zGenes.tsv`

##### i) Make a list of W gene IDs from BLAST output table 
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/PAR/blast_results
awk 'BEGIN{FS="\t"} !/^#/ && $1 != "*" {print $1}' PAR_wGenes_zGenes.tsv | sort -u > ../gene_lists/PAR-W-genes_transcripts.txt
# Convert trasncript IDs to gene IDs (remove -mRNA flag)
sed 's/-mRNA-.*//' ../gene_lists/PAR-W-genes_transcripts.txt | sort -u > ../gene_lists/PAR-W-genes.txt
rm ../gene_lists/PAR-W-genes_transcripts.txt

# The W genes in the final reference genome are labeled as LOCW (instead of just LOC), so we have to make a list of the W genes with the correct LOC labeling
sed 's/^LOC_/LOCW_/' ../gene_lists/PAR-W-genes.txt > ../gene_lists/PAR-W-genes_revised.txt
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/PAR/gene_lists/PAR-W-genes.txt`

### 9) Find gene coordinates
Z-genes
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/PAR/gene_coords
zgrep -F -f /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/PAR/gene_lists/PAR-Z-genes.txt /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/references/Heloderma_suspectum.final.reference.gff3.gz | awk '$3=="gene"' > Z_PAR_coords.gff3
# Add this header: #seqid  source  type  start  end  score  strand  phase  attributes
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/PAR/gene_coords/Z_PAR_coords.gff3`

W-genes
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/PAR/gene_coords
zgrep -F -f /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/PAR/gene_lists/PAR-W-genes_revised.txt /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/references/Heloderma_suspectum.final.reference.gff3.gz | awk '$3=="gene"' > W_PAR_coords.gff3
# Add this header: #seqid  source  type  start  end  score  strand  phase  attributes
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/PAR/gene_coords/W_PAR_coords.gff3`

### 10) Compare Z and W genes that returned hits to get a list of PAR genes
Make table of Z-W PAR best hits. Run R scrpit:
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/scripts
module load R
Rscript 10_PAR_best_hits.R
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/PAR/Z_W_PARs.tsv`

This table contains only Z and W genes that were reciprocal best hits in the PAR. The table contains the name and start and end coordinates for each gene.

### 11) Extract Z and W gametologs from FASTA files
Z genes
```
# Extract Gila Z genes from PAR table
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/PAR
tail -n +2 Z_W_PARs.tsv | cut -f5 | sed 's/-mRNA-.*//' | sort -u > gene_lists/gila_PAR_Z_genes.txt

# Using the list of Gila Z gametologs, extract peptide FASTA files
module load seqkit
seqkit grep -n -r -f gene_lists/gila_PAR_Z_genes.txt ../../references/Heloderma_suspectum.chrZ-only.peptides.fasta.gz -o gene_lists/Z_PAR_isoforms.fasta

# Keep only one isoform per gene (the longest one)
awk '
BEGIN { RS=">"; ORS="" }
NR > 1 {
    n = split($0, a, "\n")
    hdr = a[1]
    id = hdr
    sub(/ .*/, "", id)

    gene = id
    sub(/-mRNA-.*/, "", gene)

    seq = ""
    for (i = 2; i <= n; i++) seq = seq a[i]

    if (!(gene in maxlen) || length(seq) > maxlen[gene]) {
        maxlen[gene] = length(seq)
        best_hdr[gene] = hdr
        best_seq[gene] = seq
    }
}
END {
    for (g in best_seq)
        print ">" best_hdr[g] "\n" best_seq[g] "\n"
}
' gene_lists/Z_PAR_isoforms.fasta > gene_lists/Z_PARs.fasta
rm gene_lists/Z_PAR_isoforms.fasta
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/PAR/gene_lists/Z_PARs.fasta`

W genes
```
# Extract Gila W genes from gametolog table
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/PAR

tail -n +2 Z_W_PARs.tsv | cut -f1 | sed 's/^LOCW_/LOC_/' | sed 's/-mRNA-.*//' | sort -u > gene_lists/gila_PAR_W_genes.txt

# Using the list of Gila W gametologs, extract peptide FASTA files
module load seqkit
seqkit grep -n -r -f gene_lists/gila_PAR_W_genes.txt ../../references/Heloderma_suspectum.chrW_PAR.peptides.fasta.gz -o gene_lists/W_PAR_isoforms.fasta

# Keep only one isoform per gene
awk '
BEGIN { RS=">"; ORS="" }
NR > 1 {
    n = split($0, a, "\n")
    hdr = a[1]
    id = hdr
    sub(/ .*/, "", id)

    gene = id
    sub(/-mRNA-.*/, "", gene)

    seq = ""
    for (i = 2; i <= n; i++) seq = seq a[i]

    if (!(gene in maxlen) || length(seq) > maxlen[gene]) {
        maxlen[gene] = length(seq)
        best_hdr[gene] = hdr
        best_seq[gene] = seq
    }
}
END {
    for (g in best_seq)
        print ">" best_hdr[g] "\n" best_seq[g] "\n"
}
' gene_lists/W_PAR_isoforms.fasta > gene_lists/W_PARs.fasta
rm gene_lists/W_PAR_isoforms.fasta
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/gene_lists/W_PARs.fasta`

### 12) Align PAR genes
#### a) Make pairwise FASTA files for alignment (using PAR table, make a file for each row)
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only 
chmod u+x ../scripts/12a_pairwise_fastas_gila_PAR.sh
./../scripts/12a_pairwise_fastas_gila_PAR.sh
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/PAR/align/pairwise_fastas`

#### b) Align
Using PRANK (codon aware) 
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/PAR/align/
mkdir -p prank

# Align nucleotide gametolog pairs with PRANK
for f in pairwise_fastas/*.fasta; do
  base=$(basename "${f%.fasta}")
  # It is important to have the -codon flag so that the aligner is codon aware
  echo "prank -d=\"$f\" -o=prank/${base}.prank -codon -F"
done > prank/prank.swarm

# Run swarm
module load prank
swarm -f prank/prank.swarm
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/PAR/heloderma_only/align/prank`

This directory contains an alignment file for each gametolog pair, generated using PRANK.

! NOTE: Check how many files are in the output directory before proceeding, should be the same number of gametolog pairs.

! NOTE: Using the codon-aware version of PRANK removes stop codons, so none of the aligned files have stop codons

### 13) Compute dS and plot
#### a) PAML (https://hpc.nih.gov/apps/PAML.html)
Use PAML to detect natural selection of proteins. Specifically, the yn00 feature of PAML calculates pairwise comparisons between two sequences and outputs dN/dS/ω values.
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/PAR

mkdir -p paml/jobs paml/out
module load paml

# Update user permissions
chmod +x ../../scripts/13a_run_yn00_PAR.sh

# Make swarm file
find align -name '*.best.fas' | sort \
  | awk '{print "bash ../../scripts/13a_run_yn00_PAR.sh \"" $0 "\""}' > paml/paml.swarm

# Run swarm
swarm -f paml/paml.swarm
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/PAR/paml/out`

This directory contains individual txt files of yn00 outputs for each gametolog pair, made using the PRANK aligned FASTAs.

! NOTE: The Yang & Nielsen (2000) section contains the reported dN, dS, omega (dN/dS ratio) values. The other methods printed in the txt output files are ignored.

#### b) Combine each pair dN, dS, and omega values into one summary table.
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/PAR/paml/

out="paml_yn00_summary.tsv"

printf "pair\tdN\tdS\tomega\n" > "$out"

for f in out/*.yn00.txt; do
  pair=$(basename "$f" .yn00.txt)

  awk -v pair="$pair" '
    /Yang & Nielsen \(2000\)/ {flag=1; next}
    flag && $1 ~ /^[0-9]+$/ && $2 ~ /^[0-9]+$/ && $7 ~ /^[0-9.]+$/ {
      # Fields for the YN00 line:
      # $7 = omega
      # $8 = dN
      # $11 = dS
      print pair "\t" $8 "\t" $11 "\t" $7
      exit
    }
  ' "$f" >> "$out"
done
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/PAR/paml/paml_yn00_summary.tsv`

#### c) Plot dS vs. position on the Z, including both PAR and nonPAR values
Using the final PAR table, plot dS vs. position on the Z chromosome in R. Run R script:
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/scripts
module load R
13c_dS_Zcoord_plot_PAR.R
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/figures/dS_Zcoord_plot_gila_PAR.png`,
            `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/figures/dS_Zcoord_plot_gila_PAR.png`

## *Varanus acanthurus* Z v. W
#### /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/varanus_only
Now that we have looked at gametologs within the Gila monster, we want to look at the gametologs in another closely related lizard species to better understand how sex chromosome evolution occured.

### 14) Find *Varanus* Z-W gametologs
#### a) Z to W
##### i) BLAST each Z gene to all W genes (using diamond)
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/varanus_only
sbatch ../scripts/14a_diamond_zGenes_wGenes_varanus.sh
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/varanus_only/blast_results/zGenes_wGenes.tsv`

##### ii) Make a list of Z gene IDs from BLAST output table 
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/varanus_only/blast_results
awk 'BEGIN{FS="\t"} !/^#/ && $1 != "*" {print $1}' zGenes_wGenes.tsv | sort -u > ../gene_lists/Z-genes_transcripts.txt
# Convert trasncript IDs to gene IDs (remove -mRNA flag)
sed 's/-mRNA-.*//' ../gene_lists/Z-genes_transcripts.txt | sort -u > ../gene_lists/Z-genes.txt
rm ../gene_lists/Z-genes_transcripts.txt
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/varanus_only/gene_lists/Z-genes.txt`

#### b) W to Z
##### i) BLAST each W gene to all Z genes (using diamond)
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/varanus_only
sbatch ../scripts/14b_diamond_wGenes_zGenes_varanus.sh
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/varanus_only/blast_results/wGenes_zGenes.tsv`

##### ii) Make a list of W gene IDs from BLAST output table
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/varanus_only/blast_results
awk 'BEGIN{FS="\t"} !/^#/ && $2 != "*" {print $1}' wGenes_zGenes.tsv | sort -u > ../gene_lists/W-genes_transcripts.txt
# Convert trasncript IDs to gene IDs (remove -mRNA flag)
sed 's/-mRNA-.*//' ../gene_lists/W-genes_transcripts.txt | sort -u > ../gene_lists/W-genes.txt
rm ../gene_lists/W-genes_transcripts.txt
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/varanus_only/gene_lists/W-genes.txt`

### 15) Find gene coordinate
Z-genes
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/varanus_only/gene_coords

zgrep -F -f /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/varanus_only/gene_lists/Z-genes.txt /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/references/Varanus_acanthurus_chrZ.gff3.gz | awk '$3=="gene"' > Z_gametolog_coords.gff3

# Add this header: #seqid  source  type  start  end  score  strand  phase  attributes
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/varanus_only/gene_coords/Z_gametolog_coords.gff3`

W-genes
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/varanus_only/gene_coords

zgrep -F -f /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/varanus_only/gene_lists/W-genes.txt /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/references/Varanus_acanthurus_chrW.gff3.gz | awk '$3=="gene"' > W_gametolog_coords.gff3

# Add this header: #seqid  source  type  start  end  score  strand  phase  attributes
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/heloderma_only/gene_coords/W_gametolog_coords.gff3`

### 16) Compare Z and W genes that returned hits to get a list of high-confidence Z-W *Varanus* gametologs
Make tables of Z-W best hits. Run R script:
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/scripts
module load R
Rscript 16_varanus_gametologs_best_hits.R
```
**Output:** The first table contains all Z genes that had a hit to a W gene. It also lists if the Z gene had a second hit to another W gene. For each gene, it lists the evalue and alignment length. The table also lists if the Z-W hit is also the strongest W-Z hit (TRUE/FALSE).
`/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/varanus_only/Z_W_summary.tsv`,
The second table contains all W genes that had a hit to a Z gene. It also lists if the W gene had a second hit to another Z gene. For each gene, it lists the evalue and alignment length. The table also lists if the W-Z hit is also the best Z-W hit (TRUE/FALSE).
`/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/varanus_only/W_Z_summary.tsv`, 
The third table contains only Z and W genes that were reciprocal best hits. These are the most high-confidence gametolog before aligning. The table contains the name and start and end coordinates for each gene.
`/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/varanus_only/Z_W_gametologs.tsv`

### 17) Extract Z and W gametologs from FASTA files
Z genes
```
# Extract Varanus Z genes from gametolog table
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/varanus_only

tail -n +2 Z_W_gametologs.tsv | cut -f5 | sed 's/-mRNA-.*//' | sort -u > gene_lists/varanus_gametolog_Z_genes.txt

# Using the list of Varanus Z gametologs, extract peptide FASTA files
module load seqkit
seqkit grep -n -r -f gene_lists/varanus_gametolog_Z_genes.txt ../references/Varanus_acanthurus_chrZ.peptides.fasta.gz -o gene_lists/Z_gametologs_isoforms.fasta

# Keep only one isoform per gene
awk '
BEGIN { RS=">"; ORS="" }
NR > 1 {
    n = split($0, a, "\n")
    hdr = a[1]
    id = hdr
    sub(/ .*/, "", id)

    gene = id
    sub(/-mRNA-.*/, "", gene)

    seq = ""
    for (i = 2; i <= n; i++) seq = seq a[i]

    if (!(gene in maxlen) || length(seq) > maxlen[gene]) {
        maxlen[gene] = length(seq)
        best_hdr[gene] = hdr
        best_seq[gene] = seq
    }
}
END {
    for (g in best_seq)
        print ">" best_hdr[g] "\n" best_seq[g] "\n"
}
' gene_lists/Z_gametologs_isoforms.fasta > gene_lists/Z_gametologs.fasta
rm gene_lists/Z_gametologs_isoforms.fasta
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/varanus_only/gene_lists/Z_gametologs.fasta`

W genes
```
# Extract Varanus W genes from gametolog table
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/varanus_only

tail -n +2 Z_W_gametologs.tsv | cut -f5 | sed 's/-mRNA-.*//' | sort -u > gene_lists/varanus_gametolog_W_genes.txt

# Using the list of Varanus W gametologs, extract peptide FASTA files
module load seqkit
seqkit grep -n -r -f gene_lists/varanus_gametolog_W_genes.txt ../references/Varanus_acanthurus_chrW.peptides.fasta.gz -o gene_lists/W_gametologs_isoforms.fasta

# Keep only one isoform per gene
awk '
BEGIN { RS=">"; ORS="" }
NR > 1 {
    n = split($0, a, "\n")
    hdr = a[1]
    id = hdr
    sub(/ .*/, "", id)

    gene = id
    sub(/-mRNA-.*/, "", gene)

    seq = ""
    for (i = 2; i <= n; i++) seq = seq a[i]

    if (!(gene in maxlen) || length(seq) > maxlen[gene]) {
        maxlen[gene] = length(seq)
        best_hdr[gene] = hdr
        best_seq[gene] = seq
    }
}
END {
    for (g in best_seq)
        print ">" best_hdr[g] "\n" best_seq[g] "\n"
}
' gene_lists/W_gametologs_isoforms.fasta > gene_lists/W_gametologs.fasta
rm gene_lists/W_gametologs_isoforms.fasta
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/varanus_only/gene_lists/W_gametologs.fasta`

### 18) Align gametolog partners
#### a) Make pairwise FASTA files for alignment (using gametolog table, make a file for each row)
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/varanus_only
chmod u+x ../scripts/18a_pairwise_fastas_varanus.sh
./../scripts/18a_pairwise_fastas_varanus.sh
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/varanus_only/align/pairwise_fastas/`

#### b) Align
Using PRANK (codon aware)
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/varanus_only/align/
mkdir -p prank

# Align ortholog pairs with PRANK
for f in pairwise_fastas/*.fasta; do
  base=$(basename "${f%.fasta}")
  # It is important to have the -codon flag so that the aligner is codon aware
  echo "prank -d=\"$f\" -o=prank/${base}.prank -codon -F"
done > prank/prank.swarm

# Run swarm
module load prank
swarm -f prank/prank.swarm
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/varanus_only/align/prank`

! NOTE: There are 10 pairs that do not have nucleotides in frame (the sequence length is not a multiple of 3). These pairs are listed in `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/varanus_only/align/exclude_pairs.tsv`. I will exclude these pairs in downstream steps. 

#### c) Curate alignments
After running the alignments through PRANK, Brendan curated the alignments to exclude any with large gaps or other problems. This resulted in 49/52 alignments kept and 3/52 discarded. Of the 49 kept alignments, he curated 11 of them. I want to add all the best alignments to the same directory to make downstream anaylsis more straightforward:
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/varanus_only/align/

# Make new directory to keep all the best alignments
mkdir -p Varanus_pairwise_gametologs_curated_best

# Rename all files to remove "(modified)" tag
for f in Varanus_pairwise_gametologs_curated/*.fasta.gz
do
    base=$(basename "$f")
    newname="${base/ (modified)/}"
    newname="${newname%.fasta}"

    cp "$f" "Varanus_pairwise_gametologs_curated_best/$newname"
done
```
The final directroy of high-confidence alignments can be found at `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/varanus_only/align/Varanus_pairwise_gametologs_curated_best`.

### 19) Compute dS and plot
#### a) PAML 
Use PAML to detect natural selection of proteins. Specifically, the yn00 feature of PAML calculates pairwise comparisons between two sequences and outputs dN/dS/ω values.
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/varanus_only/

mkdir -p paml/jobs paml/out
module load paml

# Update user permissions
chmod +x ../scripts/19a_run_yn00_varanus.sh

# Make into swarm file
# Find all PRANK codon alignments
find align/Varanus_pairwise_gametologs_curated_best -name '*.fasta.gz' | sort \
  | awk '{print "./../scripts/19a_run_yn00_varanus.sh \"" $0 "\""}' > paml/paml.swarm

# Run swarm
swarm -f paml/paml.swarm
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/varanus_only/paml/out`

This directory contains individual txt files of yn00 outputs for each gametolog pair, made using the PRANK aligned FASTAs.

! NOTE: The Yang & Nielsen (2000) section contains the reported dN, dS, omega (dN/dS ratio) values. The other methods printed in the txt output files are ignored.

#### b) Combine each pair dN, dS, and omega values into one summary table.
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/varanus_only/paml/

out="paml_yn00_summary.tsv"

printf "pair\tdN\tdS\tomega\n" > "$out"

for f in out/*.yn00.txt; do
  pair=$(basename "$f" .yn00.txt)

  awk -v pair="$pair" '
    /Yang & Nielsen \(2000\)/ {flag=1; next}
    flag && $1 ~ /^[0-9]+$/ && $2 ~ /^[0-9]+$/ && $7 ~ /^[0-9.]+$/ {
      # Fields for the YN00 line:
      # $7 = omega
      # $8 = dN
      # $11 = dS
      print pair "\t" $8 "\t" $11 "\t" $7
      exit
    }
  ' "$f" >> "$out"
done
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/varanus_only/paml/paml_yn00_summary.tsv`

#### c) Finalize *Varanus acanthurus* Z-W gametolog list
We have to remove the 10 failed and 3 discarded alignments and add the pairwise dS values.
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/varanus_only

awk 'BEGIN{FS=OFS="\t"}
# Reads paml summary table and takes dS column (third column)
NR==FNR {
    if (FNR > 1) {
        dS[$1] = $3
    }
    next
}
NR==1 {
    print $0, "dS"
    next
}
$1!="LOC_00019093" &&
$1!="LOC_00019197" &&
$1!="LOC_00019322" &&
$1!="LOC_00019069" &&
$1!="LOC_00019131" &&
$1!="LOC_00019137" &&
$1!="LOC_00019160" &&
$1!="LOC_00019166" &&
$1!="LOC_00019193" &&
$1!="LOC_00019218" &&
$1!="LOC_00019233" &&
$1!="LOC_00019323" &&
$1!="LOC_00019378" {
    pair = $1 "__" $5
    print $0, (pair in dS ? dS[pair] : "dS")
}' paml/paml_yn00_summary.tsv Z_W_gametologs.tsv > varanus_Z_W_gametologs_final.tsv
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/varanus_only/varanus_Z_W_gametologs_final.tsv`

#### d) Plot dS vs. position on the Z
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/scripts
module load R
Rscript 19d_dS_Zcoord_plot_varanus.R
```
**Output:** Plot dS v. chrZ postion (Mb)
        `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/figures/dS_Zcoord_plot_varanus.png`, 
        `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/figures/dS_Zcoord_plot_varanus.pdf`,
        Plot dS v. chrW position (Mb)
        `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/figures/dS_Wcoord_plot_varanus.png`,
        `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/figures/dS_Wcoord_plot_varanus.pdf`,
        Plot dS v. **Heloderma chrZ position (Mb)**
        `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/figures/dS_Zcoord_plot_varanus_gilaZ.png`,
        `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/figures/dS_Zcoord_plot_varanus_gilaZ.pdf`

## Find orthologs of outgroups (using OrthoFinder outputs) & make gene trees  
#### /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/outgroups
Now that we have a list of gametologs for *Heloderma suspectum* and a close relative, *Varanus acanthurus* (they are both anguimorphs), we want to see how other orthologous genes from outgroup species might be related. The end goal of this section is to create gene trees to see how orthologs group together (might follow species trees, might not).

### 20) Make a list of orthologs of *Heloderma* Z gametolog genes for each species
First, we want to use the list of *Heloderma* Z gametologs to find orthologs in all outgroups (using Orthofinder tables). In R:
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/scripts
module load R
Rscript 20_orthofinder_gametologs.R
```
**Output**: `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/outgroups/orthologs`
This script outputs tables for each outgroup that list any orthologs from that species to any Gila Z gametologs.

### 21) Make table grouping orthologs by *Heloderma* Z, *Heloderma* W, *Varanus* Z, and *Varanus* W genes
Next, the goal is to make a table containing all species (and *Heloderma* and *Varanus* W genes) of each *Heloderma* Z gametolog that has an ortholog in each species. 

! NOTE: The outputs of Orthofinder often yield multiple genes for the same *Heloderma* gene, as either isoforms of the same gene or paralogs. We decided to use the longest transcript when deciding which ortholog to keep, so first we must find the length of each transcript: 
```
# Make tables for each species listing the trancript length of each gene (to be used later)
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/outgroups/ 

module load seqkit

mkdir -p transcript_lengths

for species in \
    Anniella_stebbinsi_HiFi_2024.asm.hic.hap2 \
    Furcifer_GCA_030440675.1_ASM3044067v1 \
    GCA_021292165.1_IOZ_Scro_1.0 \
    GCF_009769535.1_rThaEle1.pri \
    GCF_023053635.1_rElgMul1.1.pri \
    Heloderma_suspectum.final \
    GCF_027172205.1_rPodRaf1.pri \
    GCF_035149785.1_rCanAsp1.hap2 \
    GCF_037176765.1_rAnoSag1.mat \
    Varanus_acanthurus_GCA_050042745.1
do
    out=$(basename "$species")

    seqkit fx2tab -n -l "../../Heloderma_orthofinder/${species}_transcripts.fasta.gz" |
    awk 'BEGIN{FS=OFS="\t"}
    {
        id=$1
        sub(/ .*/, "", id)

        gene=id
        sub(/-mRNA-.*/, "", gene)

        print id, gene, $2
    }' > "transcript_lengths/${out}_transcript_lengths.tsv"
done

# Separate steps for Heloderma and Varanus W gametologs (their transcript fastas are in a different directory and follow a slightly different naming scheme)
for species in \
    Heloderma_suspectum.chrW-only_PAR-masked \
    Varanus_acanthurus/Varanus_acanthurus_chrW
do
    out=$(basename "$species")

    seqkit fx2tab -n -l "../../${species}.transcripts.fasta.gz" |
    awk 'BEGIN{FS=OFS="\t"}
    {
        id=$1
        sub(/ .*/, "", id)

        gene=id
        sub(/-mRNA-.*/, "", gene)

        print id, gene, $2
    }' > "transcript_lengths/${out}_transcript_lengths.tsv"

done

# Will need to rename files in `transcript_lengths` to "species"_transcript_lengths.tsv manually
```

Then, we will build a table that lists the names of the genes (specifically, the longest isoform) that are orthologs to the same *Heloderma* Z-W and *Varanus* Z-W genes. The table will be used to pull out the FASTA sequences in next steps. In R:
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/scripts
module load R
Rscript 21_merge_gametologs.R
```
**Output:** This table lists all orthologs of the shared Heloderma/Varanus Z/W genes (some species will have N/A).
            `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/iqtree/combined_table.tsv`,
            This table only lists orthologs that are present in ALL species (no N/As).
            `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/iqtree/combined_table_all_species.tsv`,
            This table only lists orthologs that are present in all species EXCEPT Shinisaurus (only Shinisaurus will have N/As).
            `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/iqtree/combined_table_no_shinisaurus.tsv`

! NOTE: Some species will have an N/A as their gene name if there is not an ortholog for the *Heloderma*/*Varanus* genes in that species. That is ok.

**IMPORTANT:** We are going to keep genes under two conditions; we will take all the genes that have orthologs in __all__ species. Then, we will also take all the genes that have orthologs in all species __except__ for *Shinisaurus*. *Shinisaurus* is not well annotated, meaning that a lot of its genes are missing, but we don't want to exclude these genes completely from the analysis if they will be informative.

### 22) Extract one FASTA per row
Using the table that has genes in all species and genes in all species but *Shinisaurus*, build one FASTA file per row.

#### i) All species
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/iqtree
mkdir -p iqtree_fastas/all_species
tail -n +2 combined_table_all_species.tsv | awk -F'\t' -v outdir="iqtree_fastas/all_species" '{
  print "bash /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/scripts/22_merge_fastas.sh \"" $1 "\" \"" $2 "\" \"" $3 "\" \"" $4 "\" \"" $5 "\" \"" $6 "\" \"" $7 "\" \"" $8 "\" \"" $9 "\" \"" $10 "\" \"" $11 "\" \"" $12 "\" \"" outdir "\"" 
}' > iqtree_fastas/all_species/merge_fastas.swarm

# Run swarm
module load seqkit
swarm -f iqtree_fastas/all_species/merge_fastas.swarm
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/iqtree/iqtree_fastas/all_species`

#### ii) All species except *Shinisaurus*
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/iqtree
mkdir -p iqtree_fastas/no_shinisaurus
tail -n +2 combined_table_no_shinisaurus.tsv | awk -F'\t' -v outdir="iqtree_fastas/no_shinisaurus" '{
  print "bash /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/scripts/22_merge_fastas.sh \"" $1 "\" \"" $2 "\" \"" $3 "\" \"" $4 "\" \"" $5 "\" \"" $6 "\" \"" $7 "\" \"" $8 "\" \"" $9 "\" \"" $10 "\" \"" $11 "\" \"" $12 "\" \"" outdir "\"" 
}' > iqtree_fastas/no_shinisaurus/merge_fastas.swarm

# Run swarm
module load seqkit
swarm -f iqtree_fastas/no_shinisaurus/merge_fastas.swarm
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/iqtree/iqtree_fastas/no_shinisaurus`

! NOTE: Each file is named by the Gila monster Z gametolog that is orthologous to all other genes in the file.

### 23) Align
#### a) Align using PRANK (codon aware)
#### i) All species
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/iqtree/
mkdir -p align/all_species

# Align ortholog pairs with PRANK
for f in iqtree_fastas/all_species/*.fasta; do
  base=$(basename "${f%.fasta}")
  # It is important to have the -codon flag so that the aligner is codon aware
  echo "prank -d=\"$f\" -o=align/all_species/${base}.prank -codon -F"
done > align/all_species/align.swarm

# Run swarm
module load prank
swarm -f align/all_species/align.swarm
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/iqtree/align/all_species`

#### ii) No *Shinisaurus*
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/iqtree/
mkdir -p align/no_shinisaurus

# Align ortholog pairs with PRANK
for f in iqtree_fastas/no_shinisaurus/*.fasta; do
  base=$(basename "${f%.fasta}")
  # It is important to have the -codon flag so that the aligner is codon aware
  echo "prank -d=\"$f\" -o=align/no_shinisaurus/${base}.prank -codon -F"
done > align/no_shinisaurus/align.swarm

# Run swarm
module load prank
swarm -f align/no_shinisaurus/align.swarm
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/iqtree/align/no_shinisaurus`

#### b) Alignment curation
After aligning, I sent the files to Brendan for alignment curation. 
##### i) All species
There were 23 original gene files. He sent back 18 files for informative genes, removing 5 from further curation. The informative alignment files can be found in the directory: `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/shefferman/iqtree/align/all_species/Heloderma_gametologs_alignments-curated-bjp`

##### ii) No *Shinisaurus*
There were 9 original gene files. He sent back 8 files for informative genes, removing 1 from further curation. The informative alignment files can be found in the directory: `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/shefferman/iqtree/align/no_shinisaurus/No_Shinisuaurs_subset`

### 24) Build gene trees for each alignment
IQ-TREE will generate a gene tree given MSA files to best predict the phylogeny of each of the gene pairs.
#### a) Run IQ-TREE on each alignment
It is important to make sure that *Podarcis* is set as the root of the tree in the iqtree command (specified by -o flag)
##### i) All species
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/iqtree

# Clean up species names (remove the gene names, only leave the species name)
mkdir -p align/all_species/clean_fastas
for f in align/all_species/Heloderma_gametologs_alignments-curated-bjp/LOC_*
do
    base=$(basename "$f" .fasta)

    cat "$f" |
    awk '
    /^>/{
        sub(/^>/,"")
        split($0,a,"|")
        print ">"a[1]
        next
    }
    {print}
    ' > align/all_species/clean_fastas/${base}.fasta
done

# Make swarm file that runs IQ-TREE on each set of gene groups. 
mkdir -p make_trees/all_species
for f in align/all_species/clean_fastas/LOC_*
do
    base=$(basename "$f" .fasta)
    echo "iqtree2 -s $f -B 10000 -o Podarcis --prefix make_trees/all_species/$base -T 1"
done > make_trees/all_species/make_trees.swarm

module load iqtree
swarm -f make_trees/all_species/make_trees.swarm
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/iqtree/make_trees/all_species`

##### ii) No *Shinisaurus*
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/iqtree

# Clean up species names (remove the gene names, only leave the species name)
mkdir -p align/no_shinisaurus/clean_fastas
for f in align/no_shinisaurus/No_Shinisaurus_subset/LOC_*
do
    base=$(basename "$f" .fasta.gz)

    gzip -cd "$f" |
    awk '
    /^>/{
        sub(/^>/,"")
        split($0,a,"|")
        print ">"a[1]
        next
    }
    {print}
    ' > align/no_shinisaurus/clean_fastas/${base}.fasta
done

# Make swarm file that runs IQ-TREE on each set of gene groups. 
mkdir -p make_trees/no_shinisaurus
for f in align/no_shinisaurus/clean_fastas/LOC_*
do
    base=$(basename "$f" .fasta)
    echo "iqtree2 -s $f -B 10000 -o Podarcis --prefix make_trees/no_shinisaurus/$base -T 1"
done > make_trees/no_shinisaurus/make_trees.swarm

module load iqtree
swarm -f make_trees/no_shinisaurus/make_trees.swarm
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/iqtree/make_trees/no_shinisaurus`

#### b) Visualize trees for each individual gene group
The best way to build trees is on FigTree (v.1.5.0). We can put all the .treefile files in one folder for easy download.
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/iqtree

mkdir -p gene_trees

find make_trees -name "*.treefile" -exec cp {} gene_trees/ \;
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/iqtree/gene_trees`

Alternatively, we can add all the trees to one .treefile file so that we can edit within the same file on FigTree.
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/iqtree/gene_trees

out="all_gene_trees.trees"

echo "#NEXUS" > "$out"
echo "begin trees;" >> "$out"

for f in *.treefile; do
    gene=$(basename "$f" .treefile)
    tree=$(cat "$f")
    echo "tree $gene = $tree" >> "$out"
done

echo "end;" >> "$out"
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/iqtree/gene_trees/all_gene_trees.nex`

#### c) Run tree topology tests on each alignment (https://iqtree.github.io/doc/Advanced-Tutorial)
Once the trees are made, we want to make sure that we can support our species-tree hypothesis with the proper statistics. 

##### i) All species
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/iqtree/topology_tests

# Make swarm file that runs a topology test on each tree (will output many files, only need .iqtree files)
mkdir -p all_species/stats

for f in ../align/all_species/clean_fastas/LOC_*.fasta
do
    base=$(basename "$f" .fasta)

    echo "iqtree2 -s $f -B 10000 -z all_species/topologies.treels -au -zb 10000 -zw -o Podarcis --prefix all_species/stats/${base} -T 1"
done > all_species/topology.swarm

module load iqtree
swarm -f all_species/topology.swarm

# After swarm finishes, run:
find all_species/stats -type f ! -name "*.iqtree" ! -name "*.trees" -delete
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/iqtree/topology_tests/all_species/stats`

##### ii) No *Shinisaurus*
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/iqtree/topology_tests

# Make swarm file that runs a topology test on each tree (will output many files, only need .iqtree files)
mkdir -p no_shinisaurus/stats

for f in ../align/no_shinisaurus/clean_fastas/LOC_*.fasta
do
    base=$(basename "$f" .fasta)

    echo "iqtree2 -s $f -B 10000 -z no_shinisaurus/topologies_no_shin.treels -au -zb 10000 -zw -o Podarcis --prefix no_shinisaurus/stats/${base} -T 1"
done > no_shinisaurus/topology.swarm

module load iqtree
swarm -f no_shinisaurus/topology.swarm

# After swarm finishes, run:
find no_shinisaurus/stats -type f ! -name "*.iqtree" ! -name "*.trees" -delete
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/iqtree/topology_tests/no_shinisaurus/stats`

! NOTE: The tree topolgy results can be found in the .iqtree files, under the USER TREES section.

#### d) Make table with topology test results for each gene
These files can be used to make a table of p-values to determine which hypothesis is failed to be rejected.
```
cd /data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/scripts
chmod u+x 24d_combine_topology_results.sh
./24d_combine_topology_results.sh
```
**Output:** `/data/Wilson_Lab/projects/squamates/Heloderma_sexchr/gametolog_analysis/sheffermannm/iqtree/topology_tests/tree_topology_results.tsv`

! NOTE: a p-AU value closer to 1 means the test **fails to reject** that tree as a true topology

! NOTE: In `tree_topology_results.tsv` tree 1 refers to the multiple origin tree and tree 2 referes to the independent origin tree.