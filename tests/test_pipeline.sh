#!/usr/bin/env bash
set -euo pipefail

Green='\033[0;32m'
Red='\033[0;31m'
Cyan='\033[0;36m'
NC='\033[0m'

# move to repo root
REPO_ROOT="$(git rev-parse --show-toplevel)"
cd "$REPO_ROOT/tests"

# create temporary symlink to run the pipeline
ln -s ../run_cutandrun.sh
ln -s ../build_sample_sheet.sh
ln -s ../map_reads_primary.sh
ln -s ../map_reads_spikein.sh
ln -s ../deduplicate_reads.sh
ln -s ../trim_reads.sh .

# preparing output folders
INDEX_DIR="index"
OUT_DIR="output"

rm -rf "$INDEX_DIR" "$OUT_DIR"

mkdir -p "$INDEX_DIR" "$OUT_DIR"

echo "Building Bowtie2 indexes"

bowtie2-build references/chr2R.fa "$INDEX_DIR/dromel"
bowtie2-build references/chrIV.fa "$INDEX_DIR/scer"

echo "Running pipeline"

bash build_sample_sheet.sh data

bash run_cutandrun.sh \
    sample_sheet.csv \
    "$INDEX_DIR/dromel" \
    "$INDEX_DIR/scer"

echo "Checking results..."

echo "Checking summary after filtering from fastp json report..."

EXP_TRIMMED_REP1=23456
EXP_TRIMMED_REP2=29574
OBS_TRIMMED_REP1=$(jq '.summary.after_filtering.total_reads' output/sample_rep1/trimmed/fastp_report.json)
OBS_TRIMMED_REP2=$(jq '.summary.after_filtering.total_reads' output/sample_rep2/trimmed/fastp_report.json)

if [[ "$EXP_TRIMMED_REP1" -eq "$OBS_TRIMMED_REP1" ]] && [[ "$EXP_TRIMMED_REP2" -eq "$OBS_TRIMMED_REP2" ]]; then
    echo -e "${Cyan}Summary matches expected number of processed reads${NC}"
else
    echo -e "${Red}Number of observed processed reads doesn't match expected result${NC}"
    exit 1
fi


[[ "$EXP_TRIMMED_REP1" -eq "$OBS_TRIMMED_REP1" ]] || exit 1

FILE="output/flagstat_summary.tsv"

if [[ ! -f "$FILE" ]]; then
    echo -e "${Red}ERROR: $FILE not found${NC}"
    exit 1
fi

# Expected values
cat <<EOF > /tmp/expected.tsv
sample_rep1 primary 23456   19438   82.87
sample_rep1 primary_dedup   23376   19358   82.81
sample_rep1 spikein 23456   4000    17.05
sample_rep1 spikein_dedup   23424   3968    16.94
sample_rep2 primary 29574   19552   66.11
sample_rep2 primary_dedup   29486   19464   66.01
sample_rep2 spikein 29574   10000   33.81
sample_rep2 spikein_dedup   29452   9878    33.54
EOF

echo "Checking flagstat summary..."

sed 1d $FILE | awk 'NR>1 {print $1,$2,$3,$4,$5}' | sort > /tmp/observed
awk 'NR>1 {print $1,$2,$3,$4,$5}' /tmp/expected.tsv | sort > /tmp/expected

if diff -u /tmp/expected /tmp/observed; then
    echo -e "${Cyan}flagstat_summary.tsv metrics matches expected results${NC}"
else
    echo -e "${Red}flagstat_summary.tsv differs from expected results${NC}"
    exit 1
fi

# Remove temporary files
rm -f /tmp/expected /tmp/observed
rm sample_sheet.csv run_cutandrun.sh build_sample_sheet.sh map_reads_primary.sh map_reads_spikein.sh deduplicate_reads.sh trim_reads.sh 


echo -e "${Green}Test passed${NC}"
