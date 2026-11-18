#!/bin/bash
# Rebuild the masks, then extract every 12 km domain once the masks are done.
#   bash scripts/submit_all.sh           masks + all domains
#   bash scripts/submit_all.sh --no-mask all domains only (masks already up to date)
# Serial QoS allows 2 running / 10 queued jobs per user, so one job per domain.
cd "$(dirname "$0")/.."
dep=""
if [ "${1:-}" != "--no-mask" ]; then
  mid=$(sbatch --parsable scripts/03_build_masks.sh) || exit 1
  echo "masks: job $mid"
  dep="--dependency=afterok:$mid"
fi
for d in EUR-12 NAM-12 SAM-12 AUS-12 EAS-12 AFR-12; do
  jid=$(sbatch --parsable $dep -J "mtn_${d}" scripts/04_extract_domain.sh "$d") || exit 1
  echo "$d: job $jid"
done
