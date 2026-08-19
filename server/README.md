# ColabFold server workflow

This package is for Linux servers with NVIDIA GPUs and Slurm. It is designed
for 1,000–2,000 single-chain PTGER4 candidate sequences, not for predicting
ligand-bound complexes. Keep the experimental 5YWY/7D7M structures as the
source of binding-pocket evidence.

## Recommended two-stage design

1. **MSA phase:** query/cache MSAs at low concurrency. This protects the public
   MSA server and avoids GPU jobs waiting on network work.
2. **Screen all candidates:** one model, three recycles, no Amber relaxation.
   This produces a comparable first-pass confidence and structure set while
   controlling GPU demand.
3. **Confirm selected candidates:** choose a pre-specified shortlist using the
   summary metrics and rerun only that shortlist with five models. Do not select
   solely on mean pLDDT; also inspect confidence at the 6 A pocket residues and
   structural displacement after aligning transmembrane residues.

For 2,000 candidates this is 2,000 screening models. A 10% confirmation
shortlist adds 1,000 models, rather than immediately computing 10,000 models.

## Server requirements

- Linux + Slurm + NVIDIA driver compatible with the cluster CUDA runtime.
- Apptainer/Singularity with `--nv` enabled.
- One writable project directory and one scratch directory per job.
- Internet access to the ColabFold MSA server **or** a locally managed MSA
  database. Limit MSA jobs to two concurrent jobs unless the administrator
  approves a higher limit.

## Initial setup

From the project root on the server:

```bash
mkdir -p server/containers data/derived/server_runs
apptainer pull server/containers/colabfold_1.6.2_cuda12.sif \
  docker://ghcr.io/sokrypton/colabfold:1.6.2-cuda12
```

Copy `server/config.example.env` to a private server-specific configuration
file and fill the absolute paths. Do not commit that private configuration.
Submit jobs from the project root so the Slurm log paths resolve correctly.

## Prepare candidates

Put one or more FASTA records into, for example,
`data/derived/colabfold_inputs/large_candidate_panel.fasta`, with a unique
identifier per candidate. Then split the panel into chunks of 25 sequences:

```bash
python server/scripts/split_fasta_batches.py \
  data/derived/colabfold_inputs/large_candidate_panel.fasta \
  data/derived/server_runs/input_chunks --chunk-size 25
```

This writes `manifest.tsv`. Inspect its count before submission. Chunks must
not be edited after the manifest has been created.

## Submit MSA generation, then the screen

```bash
source server/config.private.env
export RUN_LABEL=screen
export NUM_MODELS=1
export NUM_RECYCLES=3
export MAX_CONCURRENT_MSA=2
N=$(awk 'END {print NR - 1}' "$BATCH_MANIFEST")
sbatch --array=1-${N}%${MAX_CONCURRENT_MSA} server/slurm/colabfold_msa_array.sbatch
```

After all MSA tasks complete, submit prediction jobs. Set the concurrency to the
number of GPUs that the administrator allocates to this project:

```bash
sbatch --array=1-${N}%${MAX_CONCURRENT_GPU} server/slurm/colabfold_array.sbatch
```

The Slurm array task writes one output directory per chunk. It creates
`.msa_completed` and `.prediction_completed` only after each respective command
exits successfully; incomplete chunks can be resubmitted safely.

## Confirmation stage

Create a FASTA containing only the pre-specified shortlist, split it again,
then rerun with:

```bash
export RUN_LABEL=confirm
export NUM_MODELS=5
export NUM_RECYCLES=3
sbatch --array=1-N%${MAX_CONCURRENT_MSA} server/slurm/colabfold_msa_array.sbatch
# After MSA jobs complete:
sbatch --array=1-N%${MAX_CONCURRENT_GPU} server/slurm/colabfold_array.sbatch
```

Use a separate output root for confirmation. Retain the exact input FASTA,
manifest, container tag, command-line parameters, and Slurm logs for every run.

## Safety notes

- Do not launch 1,000–2,000 parallel requests against the public MSA server.
- Use `NUM_MODELS=1` only for a coarse computational screen; it is not a final
  confidence assessment.
- This workflow does not establish ligand affinity, signaling activity, or an
  active/inactive state. Those claims require the experimental structures and,
  where needed, independent functional evidence.
