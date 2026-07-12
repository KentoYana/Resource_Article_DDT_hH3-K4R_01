#!/bin/zsh
# run_palindrome_OGW1.zsh
#
# 1. Detect inverted repeats with EMBOSS palindrome using ungapped FASTA sequences.
# 2. Parse every OGW1 IR whose complete span falls within alignment positions 1-80.
# 3. Convert ungapped OGW1 coordinates to the gapped alignment coordinate system
#    used by the R analysis (the OGW1 deletion is represented by "-" at position 56).
# 4. Cluster strongly overlapping candidates by full-span overlap.
# 5. Select one representative IR per cluster reproducibly.

set -euo pipefail

SCRIPT_DIR="${0:A:h}"
cd "$SCRIPT_DIR"

INPUT_FASTA="dataset/pan-2_genomic.fasta"
PALINDROME_OUT="output/palindrome.txt"
IR_CANDIDATES_OUT="output/ir_candidates.csv"
IR_PAIRS_OUT="output/ir_pairs.csv"

# EMBOSS palindrome parameters
MIN_PAL_LEN=5
MAX_PAL_LEN=15
GAP_LIMIT=10
NUM_MISMATCHES=1

# Coordinate conversion and analysis region
OGW1_GAP_POSITION=56
REGION_START=1
REGION_END=80

# Two candidates are assigned to the same cluster when the overlap between
# their complete IR spans is at least this fraction of the shorter span.
# Example: 0.80 means >=80% overlap of the shorter full span.
SPAN_OVERLAP_THRESHOLD=0.80

if ! command -v palindrome >/dev/null 2>&1; then
  print -u2 "Error: EMBOSS palindrome is not installed or not in PATH."
  print -u2 "Install it with: brew install emboss"
  exit 1
fi

if [[ ! -f "$INPUT_FASTA" ]]; then
  print -u2 "Error: input FASTA file not found: $INPUT_FASTA"
  exit 1
fi

mkdir -p "${PALINDROME_OUT:h}"

print "Running EMBOSS palindrome..."
print "Input : $INPUT_FASTA"
print "Output: $PALINDROME_OUT"

palindrome \
  -sequence "$INPUT_FASTA" \
  -outfile "$PALINDROME_OUT" \
  -minpallen "$MIN_PAL_LEN" \
  -maxpallen "$MAX_PAL_LEN" \
  -gaplimit "$GAP_LIMIT" \
  -nummismatches "$NUM_MISMATCHES"

candidate_tmp="$(mktemp -t ogw1_ir_candidates.XXXXXX)"
selected_tmp="$(mktemp -t ogw1_ir_selected.XXXXXX)"
trap 'rm -f "$candidate_tmp" "$selected_tmp"' EXIT

print "Extracting OGW1 IR candidates within alignment positions ${REGION_START}-${REGION_END}..."

# Parse the OGW1 section. EMBOSS prints each hit as:
#   left-coordinate  left-sequence  left-coordinate
#                    match symbols
#   right-coordinate right-sequence right-coordinate
#
# The candidate table retains scoring information so representative selection
# is auditable.
awk \
  -v gap_position="$OGW1_GAP_POSITION" \
  -v region_start="$REGION_START" \
  -v region_end="$REGION_END" '
function to_alignment(x) {
  return (x >= gap_position) ? x + 1 : x
}

function emit_candidate(    ls, le, rs, re, tmp, stem_len, matches,
                            mismatches, loop_len, span_len, symbols) {
  ls = to_alignment(left_start)
  le = to_alignment(left_end)
  rs = to_alignment(right_coord_1)
  re = to_alignment(right_coord_2)

  if (ls > le) {
    tmp = ls; ls = le; le = tmp
  }
  if (rs > re) {
    tmp = rs; rs = re; re = tmp
  }

  if (ls < region_start || le > region_end ||
      rs < region_start || re > region_end) {
    return
  }

  symbols = match_symbols
  gsub(/[[:space:]]/, "", symbols)
  matches = gsub(/\|/, "", symbols)
  stem_len = length(left_sequence)
  mismatches = stem_len - matches
  loop_len = rs - le - 1
  span_len = re - ls + 1

  candidate_count++
  printf "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n", \
         candidate_count, ls, le, rs, re, stem_len, matches, mismatches, \
         loop_len, span_len
}

BEGIN {
  in_ogw1 = 0
  state = 0
  candidate_count = 0
}

/^Palindromes of:[[:space:]]+OGW1[[:space:]]*$/ {
  in_ogw1 = 1
  state = 0
  next
}

/^Palindromes of:/ {
  if (in_ogw1) exit
  next
}

!in_ogw1 { next }

$1 ~ /^[0-9]+$/ && $2 ~ /^[[:alpha:]]+$/ && $3 ~ /^[0-9]+$/ {
  if (state == 0) {
    left_start = $1 + 0
    left_sequence = $2
    left_end = $3 + 0
    state = 1
  } else if (state == 2) {
    right_coord_1 = $1 + 0
    right_sequence = $2
    right_coord_2 = $3 + 0
    emit_candidate()
    state = 0
  }
  next
}

state == 1 && index($0, "|") > 0 {
  match_symbols = $0
  state = 2
  next
}
' "$PALINDROME_OUT" > "$candidate_tmp"

candidate_count="$(wc -l < "$candidate_tmp" | tr -d '[:space:]')"
if (( candidate_count == 0 )); then
  print -u2 "Error: no OGW1 IR candidates were found within alignment positions ${REGION_START}-${REGION_END}."
  exit 1
fi

# Cluster candidates using the overlap coefficient of their complete spans:
#
#   overlap length / length of the shorter span
#
# Clustering is transitive through union-find. Within each cluster, select the
# representative using the following deterministic priority:
#   1. more matched bases
#   2. fewer mismatches
#   3. longer stem
#   4. longer loop
#   5. longer complete span
#   6. earlier coordinates (final tie-break)
awk -F, -v threshold="$SPAN_OVERLAP_THRESHOLD" '
function find_root(x) {
  while (parent[x] != x) {
    parent[x] = parent[parent[x]]
    x = parent[x]
  }
  return x
}

function unite(a, b,    ra, rb) {
  ra = find_root(a)
  rb = find_root(b)
  if (ra != rb) parent[rb] = ra
}

function min(a, b) { return (a < b) ? a : b }
function max(a, b) { return (a > b) ? a : b }

function is_better(a, b) {
  if (match_n[a] != match_n[b]) return match_n[a] > match_n[b]
  if (mismatch_n[a] != mismatch_n[b]) return mismatch_n[a] < mismatch_n[b]
  if (stem_len[a] != stem_len[b]) return stem_len[a] > stem_len[b]
  if (loop_len[a] != loop_len[b]) return loop_len[a] > loop_len[b]
  if (span_len[a] != span_len[b]) return span_len[a] > span_len[b]
  if (ls[a] != ls[b]) return ls[a] < ls[b]
  if (le[a] != le[b]) return le[a] < le[b]
  if (rs[a] != rs[b]) return rs[a] < rs[b]
  return re[a] < re[b]
}

{
  n++
  source_id[n] = $1 + 0
  ls[n] = $2 + 0
  le[n] = $3 + 0
  rs[n] = $4 + 0
  re[n] = $5 + 0
  stem_len[n] = $6 + 0
  match_n[n] = $7 + 0
  mismatch_n[n] = $8 + 0
  loop_len[n] = $9 + 0
  span_len[n] = $10 + 0
  parent[n] = n
}

END {
  for (i = 1; i <= n; i++) {
    for (j = i + 1; j <= n; j++) {
      overlap = min(re[i], re[j]) - max(ls[i], ls[j]) + 1
      if (overlap < 0) overlap = 0
      shorter = min(span_len[i], span_len[j])
      coefficient = (shorter > 0) ? overlap / shorter : 0
      if (coefficient >= threshold) unite(i, j)
    }
  }

  for (i = 1; i <= n; i++) {
    root = find_root(i)
    cluster_id[i] = root
    if (!(root in representative) || is_better(i, representative[root])) {
      representative[root] = i
    }
  }

  for (root in representative) {
    i = representative[root]
    printf "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n", \
           ls[i], le[i], rs[i], re[i], stem_len[i], match_n[i], \
           mismatch_n[i], loop_len[i], span_len[i], root, source_id[i]
  }
}
' "$candidate_tmp" \
  | sort -t, -k1,1n -k2,2n -k3,3n -k4,4n > "$selected_tmp"

# Write an audit table containing every candidate and the scoring fields used
# for representative selection.
{
  print "candidate_id,left_start,left_end,right_start,right_end,stem_length,matches,mismatches,loop_length,span_length"
  cat "$candidate_tmp"
} > "$IR_CANDIDATES_OUT"

# Final table consumed by the R script.
awk -F, 'BEGIN {
  OFS=","
  print "IR_ID,left_start,left_end,right_start,right_end"
}
{
  printf "IR%d,%d,%d,%d,%d\n", NR, $1, $2, $3, $4
}' "$selected_tmp" > "$IR_PAIRS_OUT"

selected_count="$(( $(wc -l < "$IR_PAIRS_OUT") - 1 ))"

print "Done."
print "  Candidates in region : $candidate_count"
print "  Representative IRs   : $selected_count"
print "Generated:"
print "  $PALINDROME_OUT"
print "  $IR_CANDIDATES_OUT"
print "  $IR_PAIRS_OUT"
