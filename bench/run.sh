#!/usr/bin/env bash
# Reproducible vcfclick ingest benchmark.
#
# Reproduces the per-configuration numbers in bench/BENCHMARK.md against
# the 1000 Genomes Phase 3 chr17:40M-50M slice (235,768 variants ×
# 3,202 samples, ~114 MB compressed).
#
# Usage:
#     ./bench/run.sh                     # all three configs
#     ./bench/run.sh serial              # one config only
#     ./bench/run.sh parallel-4 parallel-8
#
# Requires: uv, tabix, bgzip on PATH. Downloads ~110 MB on first run
# and caches it under bench/data/. The benchmark itself runs entirely
# in bench/.vcfclick/ — no other vcfclick database is touched.

set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
DATA_DIR="$REPO/bench/data"
SLICE="$DATA_DIR/chr17_40_50M.vcf.gz"
SLICE_URL="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr17.filtered.shapeit2-duohmm-phased.vcf.gz"
SLICE_REGION="chr17:40000000-50000000"

# Isolated home so this benchmark cannot touch the user's real
# ~/.vcfclick/dbs/. Removed between configurations for cold starts.
export VCFCLICK_HOME="$REPO/bench/.vcfclick"

# --- configurations ---
# `case` rather than an associative array — stock macOS bash is 3.2 and
# does not support `declare -A`. Keep the script portable.
CONFIG_ORDER=(serial parallel-4 parallel-8)

flags_for() {
    case "$1" in
        serial)     echo "--serial" ;;
        parallel-4) echo "--workers 4" ;;
        parallel-8) echo "--workers 8" ;;
        *)          echo "" ;;
    esac
}

# Pick which configs to run from argv, default = all. Validate first so
# a typo cannot trigger the 110 MB download below.
if [[ $# -gt 0 ]]; then
    SELECTED=("$@")
    for c in "${SELECTED[@]}"; do
        if [[ -z "$(flags_for "$c")" ]]; then
            echo "unknown config: $c (choose from: ${CONFIG_ORDER[*]})" >&2
            exit 1
        fi
    done
else
    SELECTED=("${CONFIG_ORDER[@]}")
fi

# --- tool check ---
for tool in uv tabix bgzip; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "missing: $tool — install before running this benchmark" >&2
        exit 1
    fi
done

# --- data slice (cached) ---
mkdir -p "$DATA_DIR"
if [[ ! -f "$SLICE" ]]; then
    echo "[bench] fetching $SLICE_REGION from 1000G (~110 MB, one-time)..."
    # tabix on a remote URL writes the downloaded .tbi to the current
    # working directory. cd into DATA_DIR so the side-effect file lands
    # under bench/data/ (which is gitignored) rather than the repo root.
    (
        cd "$DATA_DIR"
        tabix -h "$SLICE_URL" "$SLICE_REGION" | bgzip > "$(basename "$SLICE").tmp"
        mv "$(basename "$SLICE").tmp" "$(basename "$SLICE")"
        tabix -p vcf "$(basename "$SLICE")"
    )
fi
echo "[bench] slice: $SLICE ($(du -h "$SLICE" | cut -f1))"

# --- run ---
run_one() {
    local label="$1"
    local flags
    flags="$(flags_for "$label")"
    echo
    echo "===== $label ====="
    # Cold start: wipe the embedded DB so every config starts from zero.
    rm -rf "$VCFCLICK_HOME"
    uv run vcfclick db create bench >/dev/null
    # The ingester prints its own throughput line:
    #   [ingest] done. N variants in X.Xs (Y/s)
    # We rely on that rather than wrapping in `time` so the reported
    # numbers match what a manual reproducer would see.
    # shellcheck disable=SC2086
    uv run vcfclick db ingest bench "$SLICE" \
        --cohort bench --ingest-id "$label" $flags
    uv run vcfclick db query bench "SELECT count() FROM variants"
}

for c in "${SELECTED[@]}"; do
    run_one "$c"
done

echo
echo "[bench] done. Working dir: $VCFCLICK_HOME ($(du -sh "$VCFCLICK_HOME" 2>/dev/null | cut -f1))"
echo "[bench] safe to delete with: rm -rf $VCFCLICK_HOME"
