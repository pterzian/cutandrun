#!/bin/bash

# Input arguments
sample_name="$1"
assay_name="$2"
fastq_r1="$3"
fastq_r2="$4"
THREADS=4

# Create output directory
output_dir="output/${sample_name}/trimmed"
mkdir -p "$output_dir"

# Define output file names
trimmed_r1="${output_dir}/${sample_name}_trimmed_R1.fq.gz"
trimmed_r2="${output_dir}/${sample_name}_trimmed_R2.fq.gz"

# Run fastp
fastp --in1 "$fastq_r1" --in2 "$fastq_r2" \
      --out1 "$trimmed_r1" --out2 "$trimmed_r2" \
      --thread $THREADS \
      --detect_adapter_for_pe \
      --qualified_quality_phred 20 \
      --length_required 20 \
      --json ${output_dir}/fastp_report.json \
      --html ${output_dir}/fastp_report.html
