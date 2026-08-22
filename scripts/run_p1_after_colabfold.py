#!/usr/bin/env python3
"""One-command PocketState-AE P1 processing after an existing ColabFold run.

This driver does not change or rerun the existing server workflow.  It validates
the completed raw output, freezes the experimental pocket definition, ingests
QC-passing WT conformers, and builds the 196-column geometry feature matrix.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path

from finalize_pocket_definition import finalize_pocket_definition
from validate_colabfold_chunk import validate_chunk


PROJECT = Path(__file__).resolve().parents[1]


def manifest_chunks(manifest: Path, input_root: Path) -> list[tuple[Path, Path]]:
    with manifest.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or "input_fasta" not in rows[0]:
        raise ValueError(f"{manifest}: manifest is empty or lacks input_fasta")
    chunks = []
    for row in rows:
        fasta = manifest.parent / row["input_fasta"]
        output = input_root / Path(row["input_fasta"]).stem
        chunks.append((fasta, output))
    return chunks


def run_command(arguments: list[str]) -> None:
    print("+ " + " ".join(arguments))
    subprocess.run(arguments, check=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run PocketState-AE P1 steps after an existing ColabFold server run")
    parser.add_argument("--input-root", type=Path, required=True, help="Existing raw run: OUTPUT_ROOT/RUN_LABEL")
    parser.add_argument("--batch-manifest", type=Path, required=True, help="Manifest used by the existing server run")
    parser.add_argument("--run-label", required=True, help="Stable name for this completed run")
    parser.add_argument("--mode", choices=("sampling", "variant_screen", "confirm"), default="sampling")
    parser.add_argument("--output-root", type=Path, default=PROJECT / "data" / "processed")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    input_root = args.input_root.resolve()
    manifest = args.batch_manifest.resolve()
    if not input_root.is_dir():
        raise FileNotFoundError(f"ColabFold run directory does not exist: {input_root}")
    if not manifest.is_file():
        raise FileNotFoundError(f"Batch manifest does not exist: {manifest}")

    print("[1/4] Freeze the experimental 6 A pocket definition")
    pocket_dir = args.output_root / "pocket_definition"
    _, pocket_json = finalize_pocket_definition(
        PROJECT / "data" / "derived" / "c3_ligand_pockets" / "ep4_ligand_pocket_residues_within_6A.csv",
        pocket_dir,
    )

    print("[2/4] Validate every FASTA ID in every completed ColabFold chunk")
    failures: list[str] = []
    for fasta, chunk_output in manifest_chunks(manifest, input_root):
        results = validate_chunk(fasta, chunk_output, "prediction")
        for result in results:
            if result.ok:
                print(f"[OK] {fasta.name}: {result.identifier}")
            else:
                failures.append(f"{fasta.name}/{result.identifier}: " + "; ".join(result.reasons))
    if failures:
        print("Validation failed; raw files were not modified and ingestion did not start.", file=sys.stderr)
        for failure in failures:
            print("[FAIL] " + failure, file=sys.stderr)
        return 2

    conformer_dir = args.output_root / "conformers" / args.run_label
    print("[3/4] Ingest structures and apply pocket-centric QC")
    run_command([
        sys.executable,
        str(PROJECT / "scripts" / "ingest_colabfold_ensemble.py"),
        str(input_root),
        "--run-label", args.run_label,
        "--batch-manifest", str(manifest),
        "--pocket-definition", str(pocket_json),
        "--output-dir", str(conformer_dir),
        "--mode", args.mode,
    ])

    conformer_manifest = conformer_dir / "conformer_manifest.csv"
    if not conformer_manifest.is_file() or len(conformer_manifest.read_text(encoding="utf-8-sig").splitlines()) < 2:
        print("No WT conformer passed QC. See colabfold_models.csv and query_summary.csv; feature extraction was not run.", file=sys.stderr)
        return 3

    feature_dir = args.output_root / "features" / args.run_label
    print("[4/4] Build the fixed 196-dimensional pocket feature matrix")
    run_command([
        sys.executable,
        str(PROJECT / "scripts" / "extract_pocket_features.py"),
        "--manifest", str(conformer_manifest),
        "--pocket-definition", str(pocket_json),
        "--output-dir", str(feature_dir),
    ])
    print("P1 completed without modifying the raw ColabFold run.")
    print(f"QC: {conformer_dir}")
    print(f"Features: {feature_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
