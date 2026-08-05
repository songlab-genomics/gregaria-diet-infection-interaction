#!/usr/bin/env bash
#SBATCH --job-name=metatx_env
#SBATCH --output=slurm-metatx-env-%j.out
#SBATCH --error=slurm-metatx-env-%j.err
#SBATCH --partition=htc
#SBATCH --qos=public
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=0-02:00:00
#SBATCH --export=ALL

set -euo pipefail

# Keep the reusable Kraken2/Bracken/minimap2 environment on scratch rather than
# consuming home-directory quota.
ENV_PREFIX="${METATX_CONDA_ENV:-/scratch/mtecher/conda-envs/metatranscriptomics-qc}"
INSTALLER_VERSION="scaffold-origin-v2"

module purge
module load mamba

mkdir -p "$(dirname "${ENV_PREFIX}")"

# Prevent concurrent submissions from changing the same environment.
exec 9>"${ENV_PREFIX}.install.lock"
if ! flock -n 9; then
  echo "ERROR: another installation is using ${ENV_PREFIX}" >&2
  exit 2
fi

echo "Environment: ${ENV_PREFIX}"
echo "Installer version: ${INSTALLER_VERSION}"
echo "Required packages: kraken2 bracken minimap2"
echo "Started: $(date)"

if [ -d "${ENV_PREFIX}/conda-meta" ]; then
  echo "Updating the existing environment without deleting it."
  mamba install --yes \
    --prefix "${ENV_PREFIX}" \
    kraken2 bracken minimap2 \
    --channel conda-forge \
    --channel bioconda \
    --strict-channel-priority
else
  echo "Creating the Kraken2/Bracken/minimap2 environment."
  mamba create --yes \
    --prefix "${ENV_PREFIX}" \
    kraken2 bracken minimap2 \
    --channel conda-forge \
    --channel bioconda \
    --strict-channel-priority
fi

source activate "${ENV_PREFIX}"

command -v kraken2
kraken2 --version
command -v bracken
bracken -v 2>&1 | head -n 1 || true
command -v minimap2
minimap2 --version

# Preserve an exact package manifest for reproducibility.
mamba list --explicit --prefix "${ENV_PREFIX}" \
  > "${ENV_PREFIX}/explicit-spec.txt"

echo "Environment installation complete."
echo "Package manifest: ${ENV_PREFIX}/explicit-spec.txt"
echo "Finished: $(date)"
