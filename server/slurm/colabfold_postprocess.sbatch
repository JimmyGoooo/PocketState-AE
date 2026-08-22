#!/usr/bin/env bash
# Submit after the complete prediction array succeeds, for example:
# sbatch --dependency=afterok:<prediction_job_id> server/slurm/colabfold_postprocess.sbatch
# Adapt the CPU partition, memory, time, and log paths to your cluster policy.
#SBATCH --job-name=pocketstate-cf-post
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=server/logs/colabfold_postprocess_%j.out
#SBATCH --error=server/logs/colabfold_postprocess_%j.err

set -euo pipefail

: "${PROJECT_ROOT:?source server/config.private.env first}"
: "${BATCH_MANIFEST:?source server/config.private.env first}"
: "${OUTPUT_ROOT:?source server/config.private.env first}"
: "${RUN_LABEL:?source server/config.private.env first}"
: "${POSTPROCESS_MODE:=sampling}"
: "${PYTHON_BIN:=python}"

[[ -d "$PROJECT_ROOT" ]] || { echo "Missing project root: $PROJECT_ROOT" >&2; exit 1; }
[[ -f "$BATCH_MANIFEST" ]] || { echo "Missing manifest: $BATCH_MANIFEST" >&2; exit 1; }
[[ "$RUN_LABEL" =~ ^[A-Za-z0-9._-]+$ ]] || { echo "RUN_LABEL contains unsafe characters: $RUN_LABEL" >&2; exit 1; }
command -v "$PYTHON_BIN" >/dev/null || { echo "Python executable not found: $PYTHON_BIN" >&2; exit 1; }

validator="$PROJECT_ROOT/scripts/validate_colabfold_chunk.py"
ingester="$PROJECT_ROOT/scripts/ingest_colabfold_ensemble.py"
feature_extractor="$PROJECT_ROOT/scripts/extract_pocket_features.py"
pocket_definition="${POCKET_DEFINITION:-$PROJECT_ROOT/data/processed/pocket_definition/pocket_definition.json}"
raw_run_dir="$OUTPUT_ROOT/$RUN_LABEL"
conformer_output_dir="$PROJECT_ROOT/data/processed/conformers/$RUN_LABEL"
feature_output_dir="$PROJECT_ROOT/data/processed/features/$RUN_LABEL"

for required_file in "$validator" "$ingester" "$feature_extractor" "$pocket_definition"; do
  [[ -f "$required_file" ]] || { echo "Missing required file: $required_file" >&2; exit 1; }
done
[[ -d "$raw_run_dir" ]] || { echo "Missing ColabFold run directory: $raw_run_dir" >&2; exit 1; }

mkdir -p "$PROJECT_ROOT/server/logs" "$conformer_output_dir" "$feature_output_dir"

validated_chunks=0
{
  IFS= read -r _header || { echo "Manifest is empty: $BATCH_MANIFEST" >&2; exit 1; }
  while IFS=$'\t' read -r array_index input_name _sequence_count _minimum_length _maximum_length; do
    [[ -n "$array_index" && -n "$input_name" ]] || continue
    input_fasta="$(dirname "$BATCH_MANIFEST")/$input_name"
    chunk_output_dir="$raw_run_dir/${input_name%.fasta}"
    [[ -f "$input_fasta" ]] || { echo "Missing chunk FASTA: $input_fasta" >&2; exit 1; }
    [[ -d "$chunk_output_dir" ]] || { echo "Missing chunk output: $chunk_output_dir" >&2; exit 1; }
    "$PYTHON_BIN" "$validator" \
      --fasta "$input_fasta" \
      --output-dir "$chunk_output_dir" \
      --phase prediction
    validated_chunks=$((validated_chunks + 1))
  done
} < "$BATCH_MANIFEST"

(( validated_chunks > 0 )) || { echo "Manifest contains no chunk rows: $BATCH_MANIFEST" >&2; exit 1; }

"$PYTHON_BIN" "$ingester" \
  "$raw_run_dir" \
  --run-label "$RUN_LABEL" \
  --batch-manifest "$BATCH_MANIFEST" \
  --pocket-definition "$pocket_definition" \
  --output-dir "$conformer_output_dir" \
  --mode "$POSTPROCESS_MODE"

conformer_manifest="$conformer_output_dir/conformer_manifest.csv"
[[ -f "$conformer_manifest" ]] || { echo "Ingester did not create: $conformer_manifest" >&2; exit 1; }

"$PYTHON_BIN" "$feature_extractor" \
  --manifest "$conformer_manifest" \
  --pocket-definition "$pocket_definition" \
  --output-dir "$feature_output_dir"

touch "$conformer_output_dir/.postprocess_completed"
echo "Validated $validated_chunks chunk(s)."
echo "Conformer manifest: $conformer_manifest"
echo "Pocket features: $feature_output_dir"
