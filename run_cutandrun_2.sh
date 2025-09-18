#!/bin/bash

# run_cutandrun2.sh
# CUT&RUN pipeline extended with UMI-tools steps for both primary and spike-in genomes

# Check if required arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <sample_sheet.csv> <primary_genome_index> <spikein_genome_index>"
    exit 1
fi

# Input parameters
sample_sheet="$1"
primary_genome="$2"
spikein_genome="$3"

# Initialize summary and timing report files
summary_tsv="output/flagstat_summary.tsv"
echo -e "SampleName\tCondition\tTotalReads\tMappedReads\tPercentMapped" > "$summary_tsv"

timing_report="output/timing_report.tsv"
echo -e "SampleName\tStep\tTimeTaken(s)" > "$timing_report"

# Dictionaries to hold sample info
declare -A fastq_r1_files
declare -A fastq_r2_files
declare -A assay_names

# Read the sample sheet and group FASTQ files by sample
echo "Reading sample sheet: $sample_sheet"
while IFS="," read -r sample_name assay_name fastq_r1 fastq_r2; do
    if [[ "$sample_name" != "SampleName" ]]; then
        echo "Processing sample entry: $sample_name"
        assay_names["$sample_name"]="$assay_name"

        if [[ -n "$fastq_r1" && -n "$fastq_r2" ]]; then
            fastq_r1_files["$sample_name"]+="$fastq_r1 "
            fastq_r2_files["$sample_name"]+="$fastq_r2 "
        else
            echo "Warning: Empty or invalid FASTQ paths for sample $sample_name"
        fi
    fi
done < "$sample_sheet"

# Function to measure time for each step
measure_time() {
    local step_name="$1"
    local sample_name="$2"
    shift 2
    local start_time=$(date +%s)
    "$@"
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    echo -e "${sample_name}\t${step_name}\t${duration}" >> "$timing_report"
}

# Process each sample
for sample_name in "${!fastq_r1_files[@]}"; do
    echo "=== Processing sample: $sample_name ==="

    # Create output directories for this sample
    sample_output_dir="output/${sample_name}"
    mkdir -p "$sample_output_dir/merged_fastq" \
             "$sample_output_dir/trimmed" \
             "$sample_output_dir/umi" \
             "$sample_output_dir/mapped_primary" \
             "$sample_output_dir/mapped_spikein" \
             "$sample_output_dir/deduplicated/primary" \
             "$sample_output_dir/deduplicated/spikein" \
             "$sample_output_dir/deduplicated_umi/primary" \
             "$sample_output_dir/deduplicated_umi/spikein" \
             "$sample_output_dir/peaks"

    # Define merged file paths
    merged_r1="${sample_output_dir}/merged_fastq/${sample_name}_merged_R1.fq.gz"
    merged_r2="${sample_output_dir}/merged_fastq/${sample_name}_merged_R2.fq.gz"

    # Step 1: Merge multiple FASTQs
    echo "Merging FASTQ files for $sample_name"
    measure_time "Merging" "$sample_name" cat ${fastq_r1_files["$sample_name"]} > "$merged_r1"
    measure_time "Merging" "$sample_name" cat ${fastq_r2_files["$sample_name"]} > "$merged_r2"

    # Step 2: Trimming
    trimmed_r1="${sample_output_dir}/trimmed/${sample_name}_trimmed_R1.fq.gz"
    trimmed_r2="${sample_output_dir}/trimmed/${sample_name}_trimmed_R2.fq.gz"
    measure_time "Trimming" "$sample_name" ./trim_reads.sh "$sample_name" "${assay_names[$sample_name]}" "$merged_r1" "$merged_r2"

    # Step 3: UMI extraction
    umi_r1="${sample_output_dir}/umi/${sample_name}_umi_R1.fq.gz"
    umi_r2="${sample_output_dir}/umi/${sample_name}_umi_R2.fq.gz"
    measure_time "UMI Extraction" "$sample_name" ./umi_extract.sh "$sample_name" "$trimmed_r1" "$trimmed_r2" "$umi_r1" "$umi_r2"

    # Step 4: Primary Genome Mapping (use UMI-cleaned reads)
    primary_bam="${sample_output_dir}/mapped_primary/${sample_name}_primary_aligned.bam"
    measure_time "Primary Mapping" "$sample_name" ./map_reads_primary.sh "$sample_name" "$umi_r1" "$umi_r2" "$primary_genome"

    # Step 5: Spike-in Genome Mapping (use UMI-cleaned reads)
    spikein_bam="${sample_output_dir}/mapped_spikein/${sample_name}_spikein_aligned.bam"
    measure_time "Spike-in Mapping" "$sample_name" ./map_reads_spikein.sh "$sample_name" "$umi_r1" "$umi_r2" "$spikein_genome"

    # Step 6a: UMI Deduplication (Primary)
    dedup_umi_primary_bam="${sample_output_dir}/deduplicated_umi/primary/${sample_name}_primary_dedup_UMI.bam"
    measure_time "UMI Deduplication Primary" "$sample_name" ./umi_dedup.sh "$sample_name" "$primary_bam" "$dedup_umi_primary_bam"

    # Step 6b: UMI Deduplication (Spike-in)
    dedup_umi_spikein_bam="${sample_output_dir}/deduplicated_umi/spikein/${sample_name}_spikein_dedup_UMI.bam"
    measure_time "UMI Deduplication Spikein" "$sample_name" ./umi_dedup.sh "$sample_name" "$spikein_bam" "$dedup_umi_spikein_bam"

    # Step 7a: Standard Deduplication (Primary)
    #dedup_primary_bam="${sample_output_dir}/deduplicated/primary/${sample_name}_primary_dedup.bam"
    #measure_time "Deduplication Primary" "$sample_name" ./deduplicate_reads.sh "$sample_name" "$primary_bam" "$dedup_primary_bam"

    # Step 7b: Standard Deduplication (Spike-in)
    #dedup_spikein_bam="${sample_output_dir}/deduplicated/spikein/${sample_name}_spikein_dedup.bam"
    #measure_time "Deduplication Spikein" "$sample_name" ./deduplicate_reads.sh "$sample_name" "$spikein_bam" "$dedup_spikein_bam"

    # Step 8: Peak Calling (using standard dedup BAM by default)
    #measure_time "Peak Calling" "$sample_name" ./call_peaks.sh "$sample_name" "$dedup_primary_bam"`
done

echo "Pipeline complete. Summary file: $summary_tsv"
echo "Timing report generated: $timing_report"