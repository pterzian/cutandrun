#!/usr/bin/env bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <path_to_data_directory>"
    exit 1
fi

data_dir="$1"

[ -d "$data_dir" ] || { echo "Error: Directory '$data_dir' not found!"; exit 1; }

sample_sheet="sample_sheet.csv"
echo "SampleName,AssayName,FASTQ_R1,FASTQ_R2" > "$sample_sheet"

for sample_dir in "$data_dir"/*/; do
    sample_name=$(basename "$sample_dir")
    echo "Processing: $sample_name"

    find "$sample_dir" -type f -regextype posix-extended \
        -regex '.*(_R?1|\.R1)\.(fq|fastq)\.gz' |
    while read -r fastq_r1; do

        fastq_r2=$(echo "$fastq_r1" | sed -E 's/(_R?1|\.R1)\.(fq|fastq)\.gz/\1/; s/1/2/' )
        fastq_r2="${fastq_r1/_R1/_R2}"
        fastq_r2="${fastq_r2/_1/_2}"
        fastq_r2="${fastq_r2/.R1/.R2}"

        if [[ -f "$fastq_r2" ]]; then
            echo "$sample_name,cutandrun,$fastq_r1,$fastq_r2" >> "$sample_sheet"
        fi
    done

done

echo "Sample sheet generated: $sample_sheet"