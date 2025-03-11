#!/bin/bash

# Input arguments
sample_name="$1"
trimmed_r1="$2"
trimmed_r2="$3"
spikein_genome="$4"  # Path to the spike-in reference genome

# Create output directory
output_dir="output/${sample_name}/mapped_spikein"
mkdir -p "$output_dir"

# Define output file names
sam_file="${output_dir}/${sample_name}_spikein_aligned.sam"
bam_file="${output_dir}/${sample_name}_spikein_aligned.bam"
flagstat_file="${output_dir}/${sample_name}_spikein_flagstat.txt"
summary_tsv="output/flagstat_summary.tsv"

# Run bowtie2 to map to the spike-in genome
bowtie2 -x "$spikein_genome" -1 "$trimmed_r1" -2 "$trimmed_r2" --threads 4 --local --very-sensitive-local --no-mixed --no-discordant --phred33 --minins 10 --maxins 700 --dovetail -S "$sam_file"

# if [[ $? -ne 0 ]]; then
#     echo "Error: Bowtie2 alignment failed."
#     exit 1
# fi


# Convert SAM to sorted BAM
samtools view -bS "$sam_file" | samtools sort -o "$bam_file"

# Remove the SAM file after conversion
rm "$sam_file"

# Run flagstat and reformat results into summary file
samtools flagstat "$bam_file" > "$flagstat_file"
total_reads=$(grep -m 1 "in total" "$flagstat_file" | awk '{print $1}')
mapped_reads=$(grep "mapped (" "$flagstat_file" | awk '{print $1}')
pct_mapped=$(grep "mapped (" "$flagstat_file" | awk -F '[()%]' '{print $2}')

# Append formatted stats to summary TSV
echo -e "${sample_name}\tspikein\t${total_reads}\t${mapped_reads}\t${pct_mapped}" >> "$summary_tsv"
