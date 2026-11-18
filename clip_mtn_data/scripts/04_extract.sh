#!/bin/bash
# sbatch scripts/04_extract.sh DOMAIN DRIVER EXP "VARS" [extra 04_extract.py args]
#SBATCH -J mtn_extract
#SBATCH -A ICT26_ESP_0          # Leonardo account: change -A and -p for your HPC
#SBATCH -p lrd_all_serial
#SBATCH -c 4
#SBATCH --mem=30G
#SBATCH --time=04:00:00
#SBATCH -o logs/%x_%j.out
cd "${SLURM_SUBMIT_DIR:-.}"   # submit from the project folder: sbatch scripts/<this>.sh
export PYTHONUNBUFFERED=1 OMP_NUM_THREADS=1 HDF5_USE_FILE_LOCKING=FALSE
DOMAIN=$1; DRIVER=$2; EXP=$3; VARS=$4; shift 4
for v in $VARS; do
  "${PY:-python3}" scripts/04_extract.py --domain "$DOMAIN" --driver "$DRIVER" --exp "$EXP" --var "$v" "$@" || exit 1
done
