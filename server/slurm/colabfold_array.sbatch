#!/usr/bin/env bash
# Submit with: sbatch --array=1-N%2 server/slurm/colabfold_array.sbatch
# Values below are conservative placeholders; adapt partition, GPU, CPU, memory,
# wall time, and module name to the server's documented scheduler policy.
#SBATCH --job-name=pocketstate-cf
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=server/logs/colabfold_%A_%a.out
#SBATCH --error=server/logs/colabfold_%A_%a.err

set -euo pipefail

: "${PROJECT_ROOT:?source server/config.private.env first}"
: "${CONTAINER_IMAGE:?source server/config.private.env first}"
: "${BATCH_MANIFEST:?source server/config.private.env first}"
: "${OUTPUT_ROOT:?source server/config.private.env first}"
: "${SCRATCH_ROOT:?source server/config.private.env first}"
: "${RUN_LABEL:=screen}"
: "${NUM_MODELS:=1}"
: "${NUM_RECYCLES:=3}"

command -v apptainer >/dev/null || { echo "Apptainer is required" >&2; exit 1; }
[[ -f "$CONTAINER_IMAGE" ]] || { echo "Missing container: $CONTAINER_IMAGE" >&2; exit 1; }
[[ -f "$BATCH_MANIFEST" ]] || { echo "Missing manifest: $BATCH_MANIFEST" >&2; exit 1; }

line_number=$((SLURM_ARRAY_TASK_ID + 1))
input_name=$(awk -F '\t' -v row="$line_number" 'NR == row {print $2}' "$BATCH_MANIFEST")
[[ -n "$input_name" ]] || { echo "No manifest row for task $SLURM_ARRAY_TASK_ID" >&2; exit 1; }

input_fasta="$(dirname "$BATCH_MANIFEST")/$input_name"
output_dir="$OUTPUT_ROOT/$RUN_LABEL/${input_name%.fasta}"
scratch_dir="$SCRATCH_ROOT/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$output_dir" "$scratch_dir" "$PROJECT_ROOT/server/logs"

completion_marker="$output_dir/.prediction_completed"

if [[ -f "$completion_marker" ]]; then
  echo "Already completed: $output_dir (predict)"
  exit 0
fi

trap 'rm -rf "$scratch_dir"' EXIT
[[ -f "$output_dir/.msa_completed" ]] || { echo "MSA phase has not completed: $output_dir" >&2; exit 1; }
apptainer exec --nv \
  --bind "$PROJECT_ROOT:$PROJECT_ROOT" \
  --bind "$scratch_dir:/scratch" \
  "$CONTAINER_IMAGE" \
  colabfold_batch \
    --num-models "$NUM_MODELS" \
    --num-recycle "$NUM_RECYCLES" \
    --model-type alphafold2_ptm \
    --sort-queries-by length \
    --compile-mode tuned \
    "$input_fasta" "$output_dir"

touch "$completion_marker"
echo "Completed: $output_dir (predict)"
