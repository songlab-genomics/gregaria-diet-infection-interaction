#!/usr/bin/env bash
set -euo pipefail

# Process only the missing library while reusing the exact immutable reference
# products and STAR indexes from the completed 44-library run.
PROJECT_ROOT="${PROJECT_ROOT:-/scratch/mtecher/gregaria-diet-infection-interaction}"
CODE_DIR="${PROJECT_ROOT}/code"
PREVIOUS_RUN_ID="${PREVIOUS_RUN_ID:-host_pathogen_dual_20260727_000952}"

cd "${CODE_DIR}"

export DUAL_RUN_ID="${DUAL_RUN_ID:-host_pathogen_dual_1044_$(date +%Y%m%d_%H%M%S)}"
export DUAL_SAMPLE_TABLE="config/fatbody_dual_rnaseq_sample_1044.tsv"
export DUAL_INPUT_MODE="raw"
export DUAL_RAW_DIR="${DUAL_RAW_DIR:-/data/songlab/sequencing_data/RNAseq/mehreen}"
export DUAL_REFERENCE_CACHE_DIR="${DUAL_REFERENCE_CACHE_DIR:-${PROJECT_ROOT}/output/runs/${PREVIOUS_RUN_ID}/00-reference}"
export DUAL_RESUME_EXISTING="false"
export DUAL_REQUIRE_R="false"

echo "Checking the uploaded merged 1044 reads and cached references..."
bash pilot/check_host_pathogen_dual_hpc.sh

echo "Submitting one-library add-on run: ${DUAL_RUN_ID}"
sbatch --export=ALL snakemake.host_pathogen_dual.slurm
