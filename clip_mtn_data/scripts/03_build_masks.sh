#!/bin/bash
#SBATCH -J mtn_masks
#SBATCH -A ICT26_ESP_0          # Leonardo account: change -A and -p for your HPC
#SBATCH -p lrd_all_serial
#SBATCH -c 4
#SBATCH --mem=30G
#SBATCH --time=04:00:00
#SBATCH -o logs/%x_%j.out
cd "${SLURM_SUBMIT_DIR:-.}"   # submit from the project folder: sbatch scripts/<this>.sh
export PYTHONUNBUFFERED=1 OMP_NUM_THREADS=4
"${PY:-python3}" scripts/03_build_masks.py
