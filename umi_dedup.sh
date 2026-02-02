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
DIR=$(dirname "${output_bam}")
flagstat="${DIR}/${sample_name}_dedup_flagstat.txt"
summary_tsv="output/flagstat_summary.tsv"

echo "[$sample_name] Starting UMI deduplication..."

umi_tools dedup \
    --stdin="$input_bam" \
    --stdout="$output_bam" \
    --log="${output_bam%.bam}_dedup.log"


# Generate flagstat for deduplicated spike-in BAM
samtools flagstat "$output_bam" > "$flagstat"
total_reads=$(grep -m 1 "in total" "$flagstat" | awk '{print $1}')
mapped_reads=$(grep -m 1 "mapped (" "$flagstat" | awk '{print $1}')
pct_mapped=$(grep -m 1 "mapped (" "$flagstat" | awk -F '[()%]' '{print $2}')
echo -e "${sample_name}\tdedup\t${total_reads}\t${mapped_reads}\t${pct_mapped}" >> "$summary_tsv"


echo "[$sample_name] UMI deduplication completed: $output_bam"
