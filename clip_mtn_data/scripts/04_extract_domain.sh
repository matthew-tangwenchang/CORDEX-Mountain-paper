#!/bin/bash
# Extract daily pr and tas over mountains for every run of one domain, then build its regional series.
# Usage: sbatch scripts/04_extract_domain.sh DOMAIN      (resumable: finished years are skipped)
#SBATCH -J mtn_extract
#SBATCH -A ICT26_ESP_0          # Leonardo account: change -A and -p for your HPC
#SBATCH -p lrd_all_serial
#SBATCH -c 4
#SBATCH --mem=30G
#SBATCH --time=04:00:00
#SBATCH -o logs/%x_%j.out
set -u
cd "${SLURM_SUBMIT_DIR:-.}"   # submit from the project folder: sbatch scripts/<this>.sh
export PYTHONUNBUFFERED=1 OMP_NUM_THREADS=1 HDF5_USE_FILE_LOCKING=FALSE
PY="${PY:-python3}"   # python with the packages in requirements.txt
ROOT="${CORDEX_ROOT:-/leonardo_work/ICT26_ESP/CORDEX-CMIP6/DD}"
DOMAIN=$1
fail=0
for run in "$ROOT/$DOMAIN"/ICTP/*/*/r1i1p1f1/RegCM5-0/v1-r1; do
  exp=$(basename "$(dirname "$(dirname "$(dirname "$run")")")")
  driver=$(basename "$(dirname "$(dirname "$(dirname "$(dirname "$run")")")")")
  case "$driver" in wrong_*) echo "skip withdrawn $driver"; continue;; esac
  for v in pr tas; do
    if ! compgen -G "$run/day/$v/*.nc" >/dev/null && ! compgen -G "$run/day/$v/*/*.nc" >/dev/null; then echo "no $v files: $DOMAIN $driver $exp"; continue; fi
    echo "=== $DOMAIN $driver $exp $v  $(date +%T)"
    "$PY" scripts/04_extract.py --domain "$DOMAIN" --driver "$driver" --exp "$exp" --var "$v" || { echo "FAILED $DOMAIN $driver $exp $v"; fail=1; }
  done
done
echo "=== regional series $(date +%T)"
"$PY" scripts/06_regional_series.py --domain "$DOMAIN" || fail=1
echo "exit status $fail"
exit $fail
