#!/usr/bin/env zsh

set -euo pipefail

script_dir="${0:A:h}"
input_file="${script_dir}/output/expData_KNT-Op47-1.fasta"
output_file="${script_dir}/output/results_B36.aln"

if [[ ! -f "${input_file}" ]]; then
  print -u2 "Error: input file not found: ${input_file}"
  exit 1
fi

mkdir -p "${script_dir}/output"

clustalo \
  --infile "${input_file}" \
  --verbose \
  --outfmt clustal \
  --outfile "${output_file}" \
  --force

print "Generated: ${output_file}"
