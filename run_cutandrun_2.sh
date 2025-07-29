#!/bin/bash

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


# Dictionary to hold lists of FASTQ files per sample
declare -A fastq_r1_files
declare -A fastq_r2_files

# Check if the sample sheet is being read correctly
echo "Reading sample sheet: $sample_sheet"
# cat "$sample_sheet"  # Print the contents of the sample sheet

# Read the sample sheet and group FASTQ files by sample
while IFS="," read -r sample_name assay_name fastq_r1 fastq_r2; do
    # Skip the header line
    if [[ "$sample_name" != "SampleName" ]]; then
        echo "Processing sample: $sample_name"
        # echo "  R1: $fastq_r1"
        # echo "  R2: $fastq_r2"

        # Ensure there are no issues with whitespace or empty lines
        if [[ -n "$fastq_r1" && -n "$fastq_r2" ]]; then
            # echo "  Adding to arrays..."
            fastq_r1_files["$sample_name"]+="$fastq_r1 "
            fastq_r2_files["$sample_name"]+="$fastq_r2 "
            # echo "  Current R1 files for $sample_name: ${fastq_r1_files[$sample_name]}"
            # echo "  Current R2 files for $sample_name: ${fastq_r2_files[$sample_name]}"
        else
            echo "Warning: Empty or invalid FASTQ paths for sample $sample_name"
        fi
    fi
done < "$sample_sheet"  # Read directly from the sample sheet (no tail or piping)

# # Debug print to check if the arrays are populated correctly
# echo "FASTQ R1 Files for each sample:"
# for key in "${!fastq_r1_files[@]}"; do
#     echo "${key}: ${fastq_r1_files[$key]}"
# done

# echo "FASTQ R2 Files for each sample:"
# for key in "${!fastq_r2_files[@]}"; do
#     echo "${key}: ${fastq_r2_files[$key]}"
# done



# Process each sample
for sample_name in "${!fastq_r1_files[@]}"; do
    echo "Processing sample: $sample_name"

    # Create output directory for the current sample
    sample_output_dir="output/${sample_name}"
    mkdir -p "$sample_output_dir/merged_fastq" "$sample_output_dir/trimmed" \
             "$sample_output_dir/mapped_primary" "$sample_output_dir/mapped_spikein" \
             "$sample_output_dir/deduplicated" "$sample_output_dir/peaks"

    # Define merged file paths
    merged_r1="${sample_output_dir}/merged_fastq/${sample_name}_merged_R1.fq.gz"
    merged_r2="${sample_output_dir}/merged_fastq/${sample_name}_merged_R2.fq.gz"

    # Function to measure time for each step
    measure_time() {
        local step_name="$1"
        local start_time=$(date +%s)
        shift
        "$@"  # Run the command
        local end_time=$(date +%s)
        local duration=$((end_time - start_time))
        echo -e "${sample_name}\t${step_name}\t${duration}" >> "$timing_report"
    }

    # Step 1: Merging multiple libraries for the same sample
    echo "Merging FASTQ files for sample: $sample_name"
    measure_time "Merging" cat ${fastq_r1_files["$sample_name"]} > "$merged_r1"
    measure_time "Merging" cat ${fastq_r2_files["$sample_name"]} > "$merged_r2"

    # Step 2: Trimming
    measure_time "Trimming" ./trim_reads.sh "$sample_name" "$assay_name" "$merged_r1" "$merged_r2"
    
    # Define trimmed file paths
    trimmed_r1="${sample_output_dir}/trimmed/${sample_name}_trimmed_R1.fq.gz"
    trimmed_r2="${sample_output_dir}/trimmed/${sample_name}_trimmed_R2.fq.gz"

    # Step 3: Primary Genome Mapping
    measure_time "Primary Mapping" ./map_reads_primary.sh "$sample_name" "$trimmed_r1" "$trimmed_r2" "$primary_genome"
    
    # Define primary BAM path
    primary_bam="${sample_output_dir}/mapped_primary/${sample_name}_primary_aligned.bam"

    # Step 4: Spike-in Genome Mapping
    measure_time "Spike-in Mapping" ./map_reads_spikein.sh "$sample_name" "$trimmed_r1" "$trimmed_r2" "$spikein_genome"

    # Define spike-in BAM path
    spikein_bam="${sample_output_dir}/mapped_spikein/${sample_name}_spikein_aligned.bam"

    # Step 5: Deduplication
    measure_time "Deduplication" ./deduplicate_reads.sh "$sample_name" "$primary_bam" "$spikein_bam"
    
    # Define deduplicated BAM path
    dedup_primary_bam="${sample_output_dir}/deduplicated/${sample_name}_primary_dedup.bam"

    # Step 6: Peak Calling
    measure_time "Peak Calling" ./call_peaks.sh "$sample_name" "$dedup_primary_bam"

done

echo "Pipeline complete. Summary file: $summary_tsv"
echo "Timing report generated: $timing_report"