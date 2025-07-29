#!/bin/bash

# Check if the user provided a directory
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <path_to_data_directory>"
    exit 1
fi

# Assign the input directory to a variable
data_dir="$1"

# Check if the directory exists
if [ ! -d "$data_dir" ]; then
    echo "Error: Directory '$data_dir' not found!"
    exit 1
fi

# Set the output sample sheet file
sample_sheet="sample_sheet.csv"

# Write the header
echo "SampleName,AssayName,FASTQ_R1,FASTQ_R2" > "$sample_sheet"

# Loop through each sample directory inside the given data directory
for sample_dir in "$data_dir"/*/; do
    sample_dir=${sample_dir%/}  # Remove trailing slash
    sample_name=$(basename "$sample_dir")
    echo "Processing sample directory: $sample_name"

    # Find all _1.fq.gz files (paired-end read 1)
    find "$sample_dir" -type f -name "*_1.fq.gz" | while read -r fastq_r1; do
        # Derive the corresponding _2.fq.gz file
        fastq_r2="${fastq_r1/_1.fq.gz/_2.fq.gz}"

        # Ensure both pairs exist
        if [[ -f "$fastq_r1" && -f "$fastq_r2" ]]; then
            # Assay name placeholder (customize if needed)
            assay_name="cutandrun"

            # Write to the sample sheet
            echo "$sample_name,$assay_name,$fastq_r1,$fastq_r2" >> "$sample_sheet"
        fi
    done
done

echo "Sample sheet generated: $sample_sheet"

