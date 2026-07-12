#!/usr/bin/env zsh

# run_clustalo_texshade_analysis.zsh
# Generate Clustal-format alignments for TeXShade/inspection.
# This script does not perform HMM alignment or Python-based conservation counting.

set -euo pipefail

# =========================
# Move to script directory
# =========================

SCRIPT_DIR="${0:A:h}"
cd "$SCRIPT_DIR"

# =========================
# Directories
# =========================

DATASET_DIR="dataset"
OUTPUT_DIR="output"

mkdir -p "$OUTPUT_DIR"

# =========================
# Input files
# =========================

PAN2_GENOMIC_FILE="${DATASET_DIR}/pan-2_genomic.fasta"
PTHR_FILE="${DATASET_DIR}/PTHR20881.fasta"

# =========================
# Output files
# =========================

PAN2_CLUSTAL="${OUTPUT_DIR}/results_pan-2.aln"
PTHR_CLUSTAL="${OUTPUT_DIR}/PTHR20881.aln"

# =========================
# Check required commands
# =========================

print "Checking required commands..."

for cmd in clustalo; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        print -u2 "Error: $cmd was not found in PATH."
        exit 1
    fi
    print "$cmd: $(command -v "$cmd")"
done

print ""

# =========================
# Check required input files
# =========================

print "Checking required input files..."

for file in "$PAN2_GENOMIC_FILE" "$PTHR_FILE"; do
    if [[ ! -f "$file" ]]; then
        print -u2 "Error: required file not found: $file"
        exit 1
    fi
    print "Found: $file"
done

print ""

# =========================
# Show input summary
# =========================

print "Sequence names in $PAN2_GENOMIC_FILE:"
grep "^>" "$PAN2_GENOMIC_FILE" || true
print ""

print "Number of PTHR20881 family sequences:"
grep -c "^>" "$PTHR_FILE"
print ""

# =========================
# Generate Clustal alignments with Clustal Omega
# =========================

print "Running Clustal Omega for pan-2 genomic sequences..."

clustalo \
    --infile "$PAN2_GENOMIC_FILE" \
    --verbose \
    --outfmt clustal \
    --outfile "$PAN2_CLUSTAL" \
    --force

print "Created: $PAN2_CLUSTAL"
print ""

print "Running Clustal Omega for PTHR20881 family sequences..."

clustalo \
    --infile "$PTHR_FILE" \
    --verbose \
    --outfmt clustal \
    --outfile "$PTHR_CLUSTAL" \
    --force

print "Created: $PTHR_CLUSTAL"
print ""

print "Analysis complete."
print "PAN-2 Clustal alignment: $PAN2_CLUSTAL"
print "PTHR20881 Clustal alignment: $PTHR_CLUSTAL"
