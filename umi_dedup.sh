#!/bin/bash

# umi_dedup.sh
# Deduplicate BAM files using UMI-tools

# Check arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <sample_name> <input.bam> <output_dedup.bam>"
    exit 1
fi

sample_name="$1"
input_bam="$2"
output_bam="$3"

echo "[$sample_name] Starting UMI deduplication..."

umi_tools dedup \
    --stdin="$input_bam" \
    --stdout="$output_bam" \
    --log="${output_bam%.bam}_dedup.log"

echo "[$sample_name] UMI deduplication completed: $output_bam"
