#!/bin/bash
# bowtiemap - map paired-end R1/R2 fastq.gz files with bowtie2
set -euo pipefail

input_dir="${1:-}"
if [[ -z "$input_dir" ]]; then
	echo "Usage: $0 <input_dir>"
	exit 1
fi
if [[ ! -d "$input_dir" ]]; then
	echo "Input directory not found: $input_dir"
	exit 1
fi

index="/Users/JK/Desktop/Bow/IndexBow/Aae_index/Aae_index"
outbase="/Users/JK/Desktop/Bow/MapBow"

# Collect unique sample prefixes from files named *_R1_001.fastq.gz
mapfile -t samples < <(find "$input_dir" -type f -name "*_R1_001.fastq.gz" -print0 \
	| xargs -0 -n1 basename \
	| cut -d'_' -f1-2 \
	| sort -u)

for sample in "${samples[@]}"; do
	echo "mapping $sample"
	# Find full paths for R1 and R2 (handles nested dirs)
	read1_path=$(find "$input_dir" -type f -name "${sample}_R1_001.fastq.gz" -print -quit || true)
	read2_path=$(find "$input_dir" -type f -name "${sample}_R2_001.fastq.gz" -print -quit || true)

	outdir="${outbase}/${sample}"
	mkdir -p "$outdir"

	if [[ -n "$read1_path" && -n "$read2_path" && -f "$read1_path" && -f "$read2_path" ]]; then
		sam_out="${outdir}/${sample}.sam"
		bowtie2 -x "$index" -1 "$read1_path" -2 "$read2_path" -S "$sam_out"
	else
		echo "Missing files for $sample: R1='${read1_path:-}' R2='${read2_path:-}'"
	fi
done