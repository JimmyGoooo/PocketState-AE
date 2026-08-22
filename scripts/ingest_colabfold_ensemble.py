#!/usr/bin/env python3
"""Ingest ColabFold outputs into an auditable PocketState-AE conformer set.

The raw server output is read-only.  All models are inventoried; only WT models
that pass pocket-centric QC are copied into the AE conformer manifest.  Variant
models remain available in the full audit table but are never silently mixed
into the state-learning set.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import sys
from pathlib import Path
from typing import Any

from p1_utils import (
    atom_lookup,
    distance,
    load_scores,
    mean,
    parse_pdb_atoms,
    parse_score_name,
    pocket_positions,
    project_path,
    read_fasta,
    safe_identifier,
    sha256_file,
)


PROJECT = Path(__file__).resolve().parents[1]
DEFAULT_POCKET = PROJECT / "data" / "processed" / "pocket_definition" / "pocket_definition.json"
DEFAULT_WT = PROJECT / "data" / "derived" / "colabfold_inputs" / "PTGER4_WT.fasta"
STATUS_ORDER = {"PASS": 0, "WARN": 1, "FAIL": 2}


def manifest_sequences(manifest: Path) -> dict[str, tuple[str, Path]]:
    sequences: dict[str, tuple[str, Path]] = {}
    with manifest.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or "input_fasta" not in rows[0]:
        raise ValueError(f"{manifest}: missing input_fasta column")
    for row in rows:
        fasta = manifest.parent / row["input_fasta"]
        if not fasta.is_file():
            raise FileNotFoundError(f"Manifest input does not exist: {fasta}")
        for identifier, sequence in read_fasta(fasta):
            if identifier in sequences and sequences[identifier][0] != sequence:
                raise ValueError(f"Conflicting sequences for {identifier}")
            sequences[identifier] = (sequence, fasta)
    return sequences


def match_expected_query(query: str, expected: dict[str, tuple[str, Path]]) -> str | None:
    if query in expected:
        return query
    candidates = [identifier for identifier in expected if safe_identifier(identifier) == query]
    return candidates[0] if len(candidates) == 1 else None


def find_structure(score_file: Path, query: str, rank: int, tag: str) -> Path | None:
    suffix = f"rank_{rank:03d}_{tag}"
    candidates = [
        path for path in score_file.parent.iterdir()
        if path.is_file()
        and path.suffix.lower() in {".pdb", ".cif", ".mmcif"}
        and path.name.startswith(query + "_")
        and suffix.lower() in path.stem.lower()
    ]
    if not candidates:
        return None
    return sorted(candidates, key=lambda p: (0 if "relaxed" in p.name.lower() and "unrelaxed" not in p.name.lower() else 1, p.name))[0]


def scalar(data: dict[str, Any], key: str) -> float | None:
    value = data.get(key)
    return float(value) if isinstance(value, (int, float)) and math.isfinite(float(value)) else None


def pocket_pae(data: dict[str, Any], positions: list[int]) -> float | None:
    matrix = data.get("pae") or data.get("predicted_aligned_error")
    if not isinstance(matrix, list):
        return None
    values = [float(matrix[first - 1][second - 1]) for first in positions for second in positions if first != second]
    return mean(values)


def ca_clashes(atoms: list[Any], cutoff: float = 2.5) -> int:
    ca = sorted((atom for atom in atoms if atom.atom_name == "CA"), key=lambda atom: atom.residue_number)
    count = 0
    for index, first in enumerate(ca):
        for second in ca[index + 1:]:
            if abs(first.residue_number - second.residue_number) <= 1:
                continue
            if distance(first.xyz, second.xyz) < cutoff:
                count += 1
    return count


def qc_model(
    score_data: dict[str, Any],
    structure: Path | None,
    sequence: str,
    positions: list[int],
    args: argparse.Namespace,
) -> dict[str, Any]:
    plddt = score_data["plddt"]
    reasons: list[str] = []
    warnings: list[str] = []
    if len(plddt) != len(sequence):
        reasons.append(f"plddt_length_{len(plddt)}_expected_{len(sequence)}")
    usable_positions = [position for position in positions if position <= len(plddt)]
    pocket_values = [plddt[position - 1] for position in usable_positions]
    core_values = plddt[max(0, args.core_start - 1): min(len(plddt), args.core_end)]
    missing = [position for position in positions if position > len(plddt)]
    clashes = ""
    structure_residues: set[int] = set()
    if structure is None:
        reasons.append("missing_structure")
    elif structure.suffix.lower() != ".pdb":
        warnings.append("cif_structure_not_geometry_checked")
    else:
        atoms = parse_pdb_atoms(structure)
        structure_residues = {atom.residue_number for atom in atoms if atom.atom_name == "CA"}
        missing.extend(position for position in positions if position not in structure_residues)
        clashes = ca_clashes(atoms)
        if clashes > args.max_ca_clashes:
            reasons.append(f"ca_clashes_{clashes}")
    missing = sorted(set(missing))
    if missing:
        reasons.append("missing_pocket_positions=" + ";".join(map(str, missing)))
    pocket_mean = mean(pocket_values)
    core_mean = mean(core_values)
    fraction_ge70 = (sum(value >= 70 for value in pocket_values) / len(pocket_values)) if pocket_values else None
    pae = pocket_pae(score_data, usable_positions)
    if pocket_mean is None or pocket_mean < args.min_pocket_plddt:
        reasons.append("low_pocket_plddt")
    if fraction_ge70 is None or fraction_ge70 < args.min_pocket_fraction_ge70:
        reasons.append("low_pocket_confident_fraction")
    if core_mean is None or core_mean < args.min_core_plddt:
        warnings.append("low_7tm_core_plddt")
    if pae is None:
        warnings.append("missing_pocket_pae")
    elif pae > args.max_pocket_pae:
        warnings.append("high_pocket_pae")
    status = "FAIL" if reasons else ("WARN" if warnings else "PASS")
    return {
        "qc_status": status,
        "qc_fail_reasons": ";".join(reasons),
        "qc_warnings": ";".join(warnings),
        "sequence_length": len(sequence),
        "mean_plddt": mean(plddt),
        "core_mean_plddt": core_mean,
        "pocket_mean_plddt": pocket_mean,
        "pocket_min_plddt": min(pocket_values) if pocket_values else None,
        "pocket_fraction_ge70": fraction_ge70,
        "pocket_mean_pae": pae,
        "missing_pocket_positions": ";".join(map(str, missing)),
        "ca_clash_count": clashes,
    }


def format_row(row: dict[str, Any]) -> dict[str, Any]:
    result = dict(row)
    for key, value in result.items():
        if isinstance(value, float):
            result[key] = f"{value:.4f}"
        elif value is None:
            result[key] = ""
    return result


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = fields or (list(rows[0]) if rows else [])
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        if not fields:
            handle.write("")
            return
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(format_row(row) for row in rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Pocket-centric ingestion and QC of ColabFold ensembles")
    parser.add_argument("input_root", type=Path, help="One ColabFold run root, usually OUTPUT_ROOT/RUN_LABEL")
    parser.add_argument("--run-label", required=True)
    parser.add_argument("--batch-manifest", type=Path, required=True)
    parser.add_argument("--pocket-definition", type=Path, default=DEFAULT_POCKET)
    parser.add_argument("--wt-fasta", type=Path, default=DEFAULT_WT)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--mode", choices=("sampling", "variant_screen", "confirm"), default="sampling")
    parser.add_argument("--core-start", type=int, default=19)
    parser.add_argument("--core-end", type=int, default=345)
    parser.add_argument("--min-core-plddt", type=float, default=70.0)
    parser.add_argument("--min-pocket-plddt", type=float, default=70.0)
    parser.add_argument("--min-pocket-fraction-ge70", type=float, default=0.80)
    parser.add_argument("--max-pocket-pae", type=float, default=10.0)
    parser.add_argument("--max-ca-clashes", type=int, default=0)
    parser.add_argument("--include-warn", action=argparse.BooleanOptionalAction, default=True)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    input_root = args.input_root.resolve()
    output = args.output_dir.resolve()
    expected = manifest_sequences(args.batch_manifest.resolve())
    wt_records = read_fasta(args.wt_fasta.resolve())
    if len(wt_records) != 1:
        raise ValueError("WT FASTA must contain exactly one sequence")
    wt_sequence = wt_records[0][1]
    definition = json.loads(args.pocket_definition.read_text(encoding="utf-8"))
    positions = pocket_positions(definition)
    output.mkdir(parents=True, exist_ok=True)
    coordinates = output / "coordinates"
    coordinates.mkdir(exist_ok=True)

    model_rows: list[dict[str, Any]] = []
    seen_models: set[tuple[str, int]] = set()
    query_errors: dict[str, list[str]] = {identifier: [] for identifier in expected}
    score_files = sorted(input_root.rglob("*_scores_rank_*.json"))
    for score_file in score_files:
        parsed = parse_score_name(score_file)
        if not parsed:
            continue
        raw_query, rank, tag, model_number, seed = parsed
        query = match_expected_query(raw_query, expected)
        if query is None:
            continue
        key = (query, rank)
        if key in seen_models:
            query_errors[query].append(f"duplicate_rank_{rank:03d}")
            continue
        seen_models.add(key)
        sequence, input_fasta = expected[query]
        try:
            scores = load_scores(score_file)
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            query_errors[query].append(f"rank_{rank:03d}_score_error:{exc}")
            continue
        structure = find_structure(score_file, raw_query, rank, tag)
        qc = qc_model(scores, structure, sequence, positions, args)
        is_wt = sequence == wt_sequence
        include = is_wt and qc["qc_status"] == "PASS"
        if is_wt and args.include_warn and qc["qc_status"] == "WARN":
            include = True
        model_id = f"{safe_identifier(query)}__rank_{rank:03d}__model_{model_number or 'NA'}__seed_{seed or 'NA'}"
        row = {
            "model_id": model_id,
            "run_label": args.run_label,
            "query_id": query,
            "rank": rank,
            "model_number": model_number,
            "seed": seed,
            "source_type": "colabfold_wt" if is_wt else "colabfold_variant",
            "is_wild_type": "YES" if is_wt else "NO",
            "include_for_pocketstate_ae": "YES" if include else "NO",
            "input_fasta": project_path(input_fasta, PROJECT),
            "score_json": project_path(score_file, PROJECT),
            "score_sha256": sha256_file(score_file),
            "structure_path": project_path(structure, PROJECT) if structure else "",
            "structure_sha256": sha256_file(structure) if structure else "",
            "ptm": scalar(scores, "ptm"),
            **qc,
        }
        if include and structure:
            extension = structure.suffix.lower()
            copied = coordinates / f"{model_id}{extension}"
            shutil.copy2(structure, copied)
            row["ae_structure_path"] = project_path(copied, PROJECT)
        else:
            row["ae_structure_path"] = ""
        model_rows.append(row)

    query_rows: list[dict[str, Any]] = []
    for query, (sequence, input_fasta) in expected.items():
        candidates = [row for row in model_rows if row["query_id"] == query]
        if not candidates:
            query_errors[query].append("no_valid_score_models")
        best = sorted(
            candidates,
            key=lambda row: (
                STATUS_ORDER[row["qc_status"]],
                -(row["pocket_mean_plddt"] or -1),
                row["pocket_mean_pae"] if row["pocket_mean_pae"] is not None else 9999,
                row["rank"],
            ),
        )[0] if candidates else None
        query_rows.append({
            "query_id": query,
            "sequence_length": len(sequence),
            "is_wild_type": "YES" if sequence == wt_sequence else "NO",
            "models_found": len(candidates),
            "included_models": sum(row["include_for_pocketstate_ae"] == "YES" for row in candidates),
            "best_model_id": best["model_id"] if best else "",
            "best_qc_status": best["qc_status"] if best else "FAIL",
            "best_pocket_mean_plddt": best["pocket_mean_plddt"] if best else None,
            "errors": ";".join(query_errors[query]),
            "input_fasta": project_path(input_fasta, PROJECT),
        })

    conformers = [row for row in model_rows if row["include_for_pocketstate_ae"] == "YES"]
    conformer_fields = [
        "model_id", "run_label", "query_id", "rank", "model_number", "seed", "source_type",
        "is_wild_type", "qc_status", "ae_structure_path", "structure_sha256", "score_json",
        "score_sha256", "mean_plddt", "core_mean_plddt", "pocket_mean_plddt",
        "pocket_min_plddt", "pocket_fraction_ge70", "pocket_mean_pae", "qc_warnings",
    ]
    write_csv(output / "colabfold_models.csv", model_rows)
    write_csv(output / "query_summary.csv", query_rows)
    write_csv(output / "conformer_manifest.csv", conformers, conformer_fields)
    provenance = {
        "schema_version": "1.0",
        "run_label": args.run_label,
        "mode": args.mode,
        "input_root": str(input_root),
        "batch_manifest": project_path(args.batch_manifest, PROJECT),
        "batch_manifest_sha256": sha256_file(args.batch_manifest),
        "pocket_definition": project_path(args.pocket_definition, PROJECT),
        "pocket_definition_sha256": sha256_file(args.pocket_definition),
        "wt_fasta_sha256": sha256_file(args.wt_fasta),
        "thresholds": {
            "core_range": [args.core_start, args.core_end],
            "min_core_plddt": args.min_core_plddt,
            "min_pocket_plddt": args.min_pocket_plddt,
            "min_pocket_fraction_ge70": args.min_pocket_fraction_ge70,
            "max_pocket_pae": args.max_pocket_pae,
            "max_ca_clashes": args.max_ca_clashes,
            "include_warn": args.include_warn,
        },
        "counts": {
            "expected_queries": len(expected),
            "score_files_seen": len(score_files),
            "models_audited": len(model_rows),
            "ae_conformers": len(conformers),
            "variant_models_excluded": sum(row["source_type"] == "colabfold_variant" for row in model_rows),
        },
    }
    (output / "provenance.json").write_text(json.dumps(provenance, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"Audited {len(model_rows)} models for {len(expected)} queries")
    print(f"PocketState-AE conformers retained: {len(conformers)} (WT only)")
    print(f"Output: {output}")
    if not model_rows:
        print("ERROR: no matching ColabFold score files were ingested", file=sys.stderr)
        return 3
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
