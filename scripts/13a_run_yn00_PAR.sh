#!/bin/bash

# This code runs PAML yn00 on a pairwise codon alignment. yn00 is used to estimate synonymous and nonsynonymous substitution rates betweent two sequences. 
# *Note* PAML yn00 requires PHYLIP format input, so the script converts the aligned FASTAs to PHYLIP format before running yn00.

# Alignment file supplied as argument 
aln="$1"

# Extract pair name from PRANK file name
base=$(basename "$aln" .prank.best.fas)

# Create working directory for pairwise analysis
work="paml/jobs/$base"
mkdir -p "$work" paml/out

# Define output files
phy="$work/$base.phy"
ctl="$work/yn00.ctl"
log="$work/yn00.log"
out="$work/yn00.out"

# Convert 2-sequence PRANK codon alignment to PHYLIP (necessary format for input to yn00)
awk '
BEGIN { RS=">"; ORS="" }

NR > 1 {
    n = split($0, a, "\n")
    seq = ""

    for (i = 2; i <= n; i++) {
        gsub(/[[:space:]]/, "", a[i])
        seq = seq a[i]
    }

    count++

    if (count == 1) {
        s1 = seq
    } else if (count == 2) {
        s2 = seq
    }
}

END {
    if (count != 2) {
        print "ERROR: expected exactly 2 sequences, found " count > "/dev/stderr"
        exit 1
    }

    if (length(s1) != length(s2)) {
        print "ERROR: aligned sequences are different lengths" > "/dev/stderr"
        exit 1
    }

    print "2 " length(s1) "\n"
    print "SEQ1      " s1 "\n"
    print "SEQ2      " s2 "\n"
}
' "$aln" > "$phy"

# Create yn00 control file (needs specific parameters for pairwise analysis)
# seqfile: input PHYLIP file with 2 sequences
# outfile: output file for yn00 results
# Use all other standard parameters
cat > "$ctl" <<CTL
seqfile = $phy
outfile = $out
verbose = 1
icode = 0
weighting = 0
commonf3x4 = 0
CTL

# Run yn00 analysis
yn00 "$ctl" > "$log" 2>&1

# Save yn00 output
cp "$out" "paml/out/$base.yn00.txt"