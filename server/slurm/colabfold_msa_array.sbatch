#!/usr/bin/env bash
# Submit with: sbatch --array=1-N%2 server/slurm/colabfold_msa_array.sbatch
# Adapt the CPU partition, memory, time, and log paths to your cluster policy.
#SBATCH --job-name=pocketstate-msa
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=server/logs/colabfold_msa_%A_%a.out
#SBATCH --error=server/logs/colabfold_msa_%A_%a.err

set -euo pipefail

: "${PROJECT_ROOT:?source server/config.private.env first}"
: "${CONTAINER_IMAGE:?source server/config.private.env first}"
: "${BATCH_MANIFEST:?source server/config.private.env first}"
: "${OUTPUT_ROOT:?source server/config.private.env first}"
: "${RUN_LABEL:=screen}"

command -v apptainer >/dev/null || { echo "Apptainer is required" >&2; exit 1; }
[[ -f "$CONTAINER_IMAGE" ]] || { echo "Missing container: $CONTAINER_IMAGE" >&2; exit 1; }
[[ -f "$BATCH_MANIFEST" ]] || { echo "Missing manifest: $BATCH_MANIFEST" >&2; exit 1; }

line_number=$((SLURM_ARRAY_TASK_ID + 1))
input_name=$(awk -F '\t' -v row="$line_number" 'NR == row {print $2}' "$BATCH_MANIFEST")
[[ -n "$input_name" ]] || { echo "No manifest row for task $SLURM_ARRAY_TASK_ID" >&2; exit 1; }

input_fasta="$(dirname "$BATCH_MANIFEST")/$input_name"
output_dir="$OUTPUT_ROOT/$RUN_LABEL/${input_name%.fasta}"
mkdir -p "$output_dir" "$PROJECT_ROOT/server/logs"

if [[ -f "$output_dir/.msa_completed" ]]; then
  echo "Already completed: $output_dir (MSA)"
  exit 0
fi

apptainer exec \
  --bind "$PROJECT_ROOT:$PROJECT_ROOT" \
  "$CONTAINER_IMAGE" \
  colabfold_batch --msa-only "$input_fasta" "$output_dir"

touch "$output_dir/.msa_completed"
echo "Completed: $output_dir (MSA)"
