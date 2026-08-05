#!/bin/bash

# This script builds one combined FASTA per Heloderma Z gene family, containing a sequence from Heloderma and each ortholog from the outgroups species.

module load seqkit

# Define for building swarm script
Heloderma_Z="$1"
Heloderma_W="$2"
Varanus_Z="$3"
Varanus_W="$4"
Anniella="$5"
Anolis="$6"
Candoia="$7"
Elgaria="$8"
Furcifer="$9"
Podarcis="${10}"
Shinisaurus="${11}"
Thamnophis="${12}"
outdir="${13}"

mkdir -p "$outdir"

# Transcript FASTAs
hel_fa="../../Heloderma_orthofinder/Heloderma_suspectum.final_transcripts.fasta.gz"
hel_w_fa="../../Heloderma_suspectum.chrW-only_PAR-masked.transcripts.fasta.gz"
varanus_fa="../../Varanus_acanthurus/Varanus_acanthurus_chrZ.transcripts.fasta.gz"
varanus_w_fa="../../Varanus_acanthurus/Varanus_acanthurus_chrW.transcripts.fasta.gz"
anniella_fa="../../Heloderma_orthofinder/Anniella_stebbinsi_HiFi_2024.asm.hic.hap2_transcripts.fasta.gz"
anolis_fa="../../Heloderma_orthofinder/GCF_037176765.1_rAnoSag1.mat_transcripts.fasta.gz"
candoia_fa="../../Heloderma_orthofinder/GCF_035149785.1_rCanAsp1.hap2_transcripts.fasta.gz"
elgaria_fa="../../Heloderma_orthofinder/GCF_023053635.1_rElgMul1.1.pri_transcripts.fasta.gz"
furcifer_fa="../../Heloderma_orthofinder/Furcifer_GCA_030440675.1_ASM3044067v1_transcripts.fasta.gz"
podarcis_fa="../../Heloderma_orthofinder/GCF_027172205.1_rPodRaf1.pri_transcripts.fasta.gz"
shinisaurus_fa="../../Heloderma_orthofinder/GCA_021292165.1_IOZ_Scro_1.0_transcripts.fasta.gz"
thamnophis_fa="../../Heloderma_orthofinder/GCF_009769535.1_rThaEle1.pri_transcripts.fasta.gz"

# The table already gives a single transcript ID. This section finds the FASTA sequence whose header matches the ID exactly.
extract_exact () {
  local fasta="$1"
  local id="$2"
  local label="$3"

  gzip -cd "$fasta" | awk -v id="$id" -v label="$label" '
  BEGIN { RS=">"; ORS="" }
  NR > 1 {
      n = split($0, a, "\n")
      hdr = a[1]

      header_id = hdr
      sub(/ .*/, "", header_id)

      if (header_id == id) {
          seq = ""
          for (i = 2; i <= n; i++) seq = seq a[i]
          gsub(/[[:space:]]/, "", seq)
          print ">" label "|" header_id "\n" seq "\n"
          found = 1
          exit
      }
  }
  END {
      if (!found) {
          print "NO_MATCH\t" label "\t" id > "/dev/stderr"
      }
  }'
}

# Make the name of each combined FASTA file the name of the Heloderma Z gene.
hel_gene="${Heloderma_Z%-mRNA-*}"
out="$outdir/${hel_gene}.fasta"

add_if_present () {
  local fasta="$1"
  local id="$2"
  local label="$3"

  [[ "$id" == "NA" || -z "$id" ]] && return 0
  extract_exact "$fasta" "$id" "$label"
}

# Combine sequences into one file.
{
  add_if_present "$hel_fa"        "$Heloderma_Z"   "Heloderma_Z"
  add_if_present "$hel_w_fa"      "$Heloderma_W"   "Heloderma_W"
  add_if_present "$varanus_fa"    "$Varanus_Z"     "Varanus_Z"
  add_if_present "$varanus_w_fa"  "$Varanus_W"     "Varanus_W"
  add_if_present "$anniella_fa"   "$Anniella"      "Anniella"
  add_if_present "$anolis_fa"     "$Anolis"        "Anolis"
  add_if_present "$candoia_fa"    "$Candoia"       "Candoia"
  add_if_present "$elgaria_fa"    "$Elgaria"       "Elgaria"
  add_if_present "$furcifer_fa"   "$Furcifer"      "Furcifer"
  add_if_present "$podarcis_fa"   "$Podarcis"      "Podarcis"
  add_if_present "$shinisaurus_fa" "$Shinisaurus"  "Shinisaurus"
  add_if_present "$thamnophis_fa" "$Thamnophis"    "Thamnophis"
} > "$out"