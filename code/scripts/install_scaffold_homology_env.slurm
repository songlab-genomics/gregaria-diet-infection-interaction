#!/usr/bin/env bash
#SBATCH --job-name=install_homology
#SBATCH --output=slurm-install-homology-%j.out
#SBATCH --error=slurm-install-homology-%j.err
#SBATCH --partition=public
#SBATCH --qos=public
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=1-00:00:00
#SBATCH --export=ALL

set -euo pipefail

ENV_DIR="${HOMOLOGY_CONDA_ENV:-/scratch/mtecher/conda-envs/scaffold-homology}"

module purge
module load mamba

# BLAST+ provides blastn/update_blastdb.pl; DIAMOND searches nr_cluster_seq.
mamba create -y -p "${ENV_DIR}" \
  -c conda-forge -c bioconda \
  python=3.12 blast diamond=2.1.24

source activate "${ENV_DIR}"
blastn -version
diamond version
update_blastdb.pl --version

echo "Scaffold homology environment ready: ${ENV_DIR}"
