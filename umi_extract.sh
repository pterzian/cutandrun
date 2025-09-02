#!/bin/bash

# umi_extract.sh
# Extract UMIs from paired-end reads using UMI-tools

# Check arguments
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <sample_name> <trimmed_r1.fq.gz> <trimmed_r2.fq.gz> <umi_r1_out.fq.gz> <umi_r2_out.fq.gz>"
    exit 1
fi

sample_name="$1"
trimmed_r1="$2"
trimmed_r2="$3"
umi_r1_out="$4"
umi_r2_out="$5"

echo "[$sample_name] Starting UMI extraction..."

umi_tools extract \
    --bc-pattern=NNNNNNNN \
    --stdin="$trimmed_r1" \
    --stdout="$umi_r1_out" \
    --read2-in="$trimmed_r2" \
    --read2-out="$umi_r2_out" \
    --log="output/${sample_name}/umi/${sample_name}_umi_extract.log"

echo "[$sample_name] UMI extraction completed."
