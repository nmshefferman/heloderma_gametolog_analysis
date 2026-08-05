#!/bin/bash

# This script takes the list of Varanus Z-W gametolog pairs and creates a pairwise FASTA file for each pair, using transcript (nucleotide) sequences.

module load seqkit

pair_tsv="../varanus_only/Z_W_gametologs.tsv"
w_ref="../../Varanus_acanthurus/Varanus_acanthurus_chrW.transcripts.fasta.gz"
z_ref="../../Varanus_acanthurus/Varanus_acanthurus_chrZ.transcripts.fasta.gz"
outdir="../varanus_only/align/pairwise_fastas/"

mkdir -p "$outdir"

# Keep the longest transcript isoform for each gene
make_longest_per_gene () {
  local in_fa="$1"
  local out_fa="$2"

  gzip -cd "$in_fa" | awk '
    BEGIN { RS=">"; ORS="" }
    NR > 1 {
      n = split($0, a, "\n")
      id = a[1]
      sub(/ .*/, "", id)

      gene = id
      sub(/-mRNA-.*/, "", gene)

      seq = ""
      for (i = 2; i <= n; i++) {
        if (a[i] != "") seq = seq a[i]
      }

      if (!(gene in maxlen) || length(seq) > maxlen[gene]) {
        maxlen[gene] = length(seq)
        bestid[gene] = id
        bestseq[gene] = seq
      }
    }
    END {
      for (g in bestseq) {
        print ">" bestid[g] "\n" bestseq[g] "\n"
      }
    }
  ' > "$out_fa"
}

make_longest_per_gene "$w_ref" "W_longest.fa"
make_longest_per_gene "$z_ref" "Z_longest.fa"

# Make one pairwise FASTA per row in the gametolog table
tail -n +2 "$pair_tsv" | while IFS=$'\t' read -r W W_name W_start W_end Z Z_name Z_start Z_end; do
  out="${outdir}/${W}__${Z}.fasta"

  w_rec=$(seqkit grep -nrp "^${W}-mRNA-" W_longest.fa || true)
  z_rec=$(seqkit grep -nrp "^${Z}-mRNA-" Z_longest.fa || true)

  if [[ -z "$w_rec" || -z "$z_rec" ]]; then
    echo "Skipping ${W} / ${Z} because one sequence was missing" >&2
    continue
  fi

  {
    printf "%s\n" "$w_rec"
    printf "%s\n" "$z_rec"
  } > "$out"
done

# Clean up intermediate files
rm W_longest.fa Z_longest.fa