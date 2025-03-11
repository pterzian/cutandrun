#!/bin/bash

# Input arguments
sample_name="$1"
bam_file="$2"
output_dir="output/${sample_name}/peaks"

# Create output directory for peaks
mkdir -p "$output_dir"

# Broad peak calling
broad_peaks_output="${output_dir}/${sample_name}_broad_peaks"
macs2 callpeak -t "$bam_file" -f BAM -g hs -n "${sample_name}_broad" --broad --outdir "$output_dir"

# Narrow peak calling
narrow_peaks_output="${output_dir}/${sample_name}_narrow_peaks"
macs2 callpeak -t "$bam_file" -f BAM -g hs -n "${sample_name}_narrow" --outdir "$output_dir"
