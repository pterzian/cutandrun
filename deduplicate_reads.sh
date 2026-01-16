#!/bin/bash

# Input arguments
sample_name="$1"
primary_bam="$2"
spikein_bam="$3"

# Output directory for deduplicated BAM files
output_dir="output/${sample_name}/deduplicated"
mkdir -p "$output_dir"

# Deduplicate Primary Genome BAM
dedup_primary_bam="${output_dir}/${sample_name}_primary_dedup.bam"
metrics_primary="${output_dir}/${sample_name}_primary_dedup_metrics.txt"
flagstat_primary="${output_dir}/${sample_name}_primary_dedup_flagstat.txt"
summary_tsv="output/flagstat_summary.tsv"

# Run Picard to remove duplicates
picard MarkDuplicatesWithMateCigar \
    I="$primary_bam" \
    O="$dedup_primary_bam" \
    M="$metrics_primary" \
    REMOVE_DUPLICATES=true

# Generate flagstat for deduplicated BAM
samtools flagstat "$dedup_primary_bam" > "$flagstat_primary"
total_reads=$(grep -m 1 "in total" "$flagstat_primary" | awk '{print $1}')
mapped_reads=$(grep -m 1 "mapped (" "$flagstat_primary" | awk '{print $1}')
pct_mapped=$(grep -m 1 "mapped (" "$flagstat_primary" | awk -F '[()%]' '{print $2}')
echo -e "${sample_name}\tprimary_dedup\t${total_reads}\t${mapped_reads}\t${pct_mapped}" >> "$summary_tsv"

# Deduplicate Spike-in Genome BAM
dedup_spikein_bam="${output_dir}/${sample_name}_spikein_dedup.bam"
metrics_spikein="${output_dir}/${sample_name}_spikein_dedup_metrics.txt"
flagstat_spikein="${output_dir}/${sample_name}_spikein_dedup_flagstat.txt"

# Run Picard to remove duplicates
picard MarkDuplicatesWithMateCigar \
    I="$spikein_bam" \
    O="$dedup_spikein_bam" \
    M="$metrics_spikein" \
    REMOVE_DUPLICATES=true

# Generate flagstat for deduplicated spike-in BAM
samtools flagstat "$dedup_spikein_bam" > "$flagstat_spikein"
total_reads=$(grep -m 1 "in total" "$flagstat_spikein" | awk '{print $1}')
mapped_reads=$(grep -m 1 "mapped (" "$flagstat_spikein" | awk '{print $1}')
pct_mapped=$(grep -m 1 "mapped (" "$flagstat_spikein" | awk -F '[()%]' '{print $2}')
echo -e "${sample_name}\tspikein_dedup\t${total_reads}\t${mapped_reads}\t${pct_mapped}" >> "$summary_tsv"
