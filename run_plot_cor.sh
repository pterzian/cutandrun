#!/bin/bash

# Usage: ./deeptools_analysis.sh <study_name> <bam_dir> <bin_size> [bam_pattern]

# Check for at least 3 arguments
if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <study_name> <bam_dir> <bin_size> [bam_pattern]"
    exit 1
fi

# Input arguments
STUDY_NAME=$1
BAM_DIR=$2
BIN_SIZE=$3
BAM_PATTERN=${4:-*.bam}  # Default to *.bam if not provided

# Full path pattern to match BAM files
BAM_PATH_PATTERN="${BAM_DIR}/${BAM_PATTERN}"

# Collect BAM files
BAM_FILES=($BAM_PATH_PATTERN)
if [ ${#BAM_FILES[@]} -eq 0 ]; then
    echo "No BAM files found matching pattern '${BAM_PATH_PATTERN}'"
    exit 1
fi

# Create labels from BAM file names
LABELS=()
for bam in "${BAM_FILES[@]}"; do
    LABELS+=("$(basename "$bam" .bam)")
done

# Convert arrays to space-separated strings
BAM_FILES_STR="${BAM_FILES[@]}"
LABELS_STR="${LABELS[@]}"

echo "Running multiBamSummary..."
multiBamSummary bins \
    --bamfiles $BAM_FILES_STR \
    --binSize $BIN_SIZE \
    --labels $LABELS_STR \
    -out ${STUDY_NAME}.npz \
    --outRawCounts ${STUDY_NAME}.tab

echo "Generating PCA plot..."
plotPCA -in ${STUDY_NAME}.npz \
    -o PCA_${STUDY_NAME}.png \
    -T "PCA of read counts"

echo "Generating Pearson correlation scatterplot..."
plotCorrelation -in ${STUDY_NAME}.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation of Read Counts" \
    --whatToPlot scatterplot --colorMap RdYlBu --plotNumbers \
    -o scatterplot_${STUDY_NAME}.png \
    --outFileCorMatrix scatterplot_${STUDY_NAME}_readCounts.tab

echo "Generating Pearson correlation heatmap..."
plotCorrelation -in ${STUDY_NAME}.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation of Read Counts" \
    --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
    -o heatmap_${STUDY_NAME}.png \
    --outFileCorMatrix heatmap_${STUDY_NAME}_readCounts.tab

echo "All analyses completed for study: ${STUDY_NAME}"
