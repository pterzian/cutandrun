#!/usr/bin/env bash
set -euo pipefail

# Get repository root (parent of tests/)
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT/tests"

INDEX_DIR="tests/index"
OUT_DIR="tests/output"

mkdir -p "$INDEX_DIR"
mkdir -p "$OUT_DIR"

echo "Building Bowtie2 indexes"

bowtie2-build tests/references/chr2R.fa "$INDEX_DIR/dromel"
bowtie2-build tests/references/chrIV.fa "$INDEX_DIR/scer"

echo "Running pipeline"

bash build_sample_sheet.sh tests/data

./run_cutandrun.sh \
    sample_sheet.csv \
    "$INDEX_DIR/dromel" \
    "$INDEX_DIR/scer"
