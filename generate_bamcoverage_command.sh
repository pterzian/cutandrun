#!/bin/bash

# Usage:
# ./generate_bigwig_commands.sh <bam_dir> <bam_pattern> <normalize_method: CPM|RPKM> <blacklist.bed> <stranded: yes|no> <output_script.sh> [output_dir] [--binSize N] [--nproc N]
# Run with --help or -h to show this message.

show_help() {
    cat << EOF
Usage: $0 <bam_dir> <bam_pattern> <normalize_method: CPM|RPKM> <blacklist.bed> <stranded: stranded|unstranded> <output_script.sh> [output_dir] [--binSize N] [--nproc N]

Parameters:
  bam_dir          Directory containing BAM files.
  bam_pattern      Pattern to match BAM files (e.g., '*.bam').
  normalize_method Normalization method: CPM or RPKM.
  blacklist.bed    Path to blacklist BED file.
  stranded         'stranded' to generate stranded bigWigs (forward & reverse), 'unstranded' for unstranded.
  output_script.sh Output bash script file name to write commands.
  output_dir       (Optional) Directory to save bigWig files. Defaults to current directory.

Optional flags:
  --binSize N      Bin size for bamCoverage (default: 5).
  --nproc N        Number of processors/threads for bamCoverage (default: 2).

Example:
  $0 /data/bams '*.bam' CPM /refs/blacklist.bed stranded bigwig_commands.sh /results --binSize 10 --nproc 4

EOF
}

# Show help if no args or help flag given
if [ "$#" -eq 0 ] || [[ "$1" == "--help" ]] || [[ "$1" == "-h" ]]; then
    show_help
    exit 0
fi

set -e

# Positional arguments
BAM_DIR=$1
BAM_PATTERN=$2
NORM_METHOD=$3
BLACKLIST=$4
STRANDED=$5
OUTPUT_SCRIPT=$6
OUTPUT_DIR=${7:-.}

# Defaults
BIN_SIZE=5
NPROC=2

shift 7
while [[ $# -gt 0 ]]; do
    case "$1" in
        --binSize)
            BIN_SIZE=$2
            shift 2
            ;;
        --nproc)
            NPROC=$2
            shift 2
            ;;
        *)
            echo "Unknown parameter: $1"
            show_help
            exit 1
            ;;
    esac
done

# Validate inputs as before...
if [[ "$NORM_METHOD" != "CPM" && "$NORM_METHOD" != "RPKM" ]]; then
    echo "Error: Normalization method must be CPM or RPKM"
    exit 1
fi

if [[ "$STRANDED" != "stranded" && "$STRANDED" != "unstranded" ]]; then
    echo "Error: Stranded must be 'stranded' or 'unstranded'"
    exit 1
fi

if [ ! -f "$BLACKLIST" ]; then
    echo "Error: Blacklist file '$BLACKLIST' does not exist"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

echo "#!/bin/bash" > "$OUTPUT_SCRIPT"
echo "" >> "$OUTPUT_SCRIPT"

shopt -s nullglob
BAM_FILES=( "$BAM_DIR"/$BAM_PATTERN )

if [ ${#BAM_FILES[@]} -eq 0 ]; then
    echo "No BAM files found matching pattern '$BAM_PATTERN' in '$BAM_DIR'"
    exit 1
fi

for BAM in "${BAM_FILES[@]}"; do
    if [ ! -f "$BAM" ]; then
        echo "Skipping invalid file: $BAM"
        continue
    fi

    BASENAME=$(basename "$BAM" .bam)

    echo "# Commands for: $(basename "$BAM")" >> "$OUTPUT_SCRIPT"

    if [ "$STRANDED" == "stranded" ]; then
        echo "bamCoverage -b \"$BAM\" -p $NPROC --binSize $BIN_SIZE --ignoreForNormalization KJ947872.2 \\" >> "$OUTPUT_SCRIPT"
        echo "  --normalizeUsing $NORM_METHOD --exactScaling --extendReads --skipNonCoveredRegions \\" >> "$OUTPUT_SCRIPT"
        echo "  --blackListFileName \"$BLACKLIST\" --filterRNAstrand reverse -of bigwig \\" >> "$OUTPUT_SCRIPT"
        echo "  -o \"$OUTPUT_DIR/${BASENAME}.rv.bigwig\" &" >> "$OUTPUT_SCRIPT"
        echo "" >> "$OUTPUT_SCRIPT"

        echo "bamCoverage -b \"$BAM\" -p $NPROC --binSize $BIN_SIZE --ignoreForNormalization KJ947872.2 \\" >> "$OUTPUT_SCRIPT"
        echo "  --normalizeUsing $NORM_METHOD --exactScaling --extendReads --skipNonCoveredRegions \\" >> "$OUTPUT_SCRIPT"
        echo "  --blackListFileName \"$BLACKLIST\" --filterRNAstrand forward -of bigwig \\" >> "$OUTPUT_SCRIPT"
        echo "  -o \"$OUTPUT_DIR/${BASENAME}.fw.bigwig\" &" >> "$OUTPUT_SCRIPT"
    else
        echo "bamCoverage -b \"$BAM\" -p $NPROC --binSize $BIN_SIZE --ignoreForNormalization KJ947872.2 \\" >> "$OUTPUT_SCRIPT"
        echo "  --normalizeUsing $NORM_METHOD --exactScaling --extendReads --skipNonCoveredRegions \\" >> "$OUTPUT_SCRIPT"
        echo "  --blackListFileName \"$BLACKLIST\" -of bigwig \\" >> "$OUTPUT_SCRIPT"
        echo "  -o \"$OUTPUT_DIR/${BASENAME}.bigwig\" &" >> "$OUTPUT_SCRIPT"
    fi

    echo "" >> "$OUTPUT_SCRIPT"
done

echo "wait" >> "$OUTPUT_SCRIPT"
chmod +x "$OUTPUT_SCRIPT"

echo "Command script written to: $OUTPUT_SCRIPT"
