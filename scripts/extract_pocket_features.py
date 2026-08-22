#!/usr/bin/env python3
"""Build the fixed-width, geometry-only feature matrix for PocketState-AE.

The model input deliberately excludes pLDDT and source labels so the latent
space cannot trivially cluster by prediction quality or provenance.  Those
values remain in the separate feature_qc.csv metadata table.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from itertools import combinations
from pathlib import Path
from typing import Any

try:
    import numpy as np
except ImportError as exc:  # pragma: no cover - clear server guidance
    raise SystemExit("NumPy is required: python -m pip install -r requirements-p1.txt") from exc

from p1_utils import (
    atom_lookup,
    parse_mmcif_atoms,
    parse_pdb_atoms,
    pocket_positions,
    project_path,
    resolve_project_path,
    sha256_file,
)


PROJECT = Path(__file__).resolve().parents[1]
DEFAULT_POCKET = PROJECT / "data" / "processed" / "pocket_definition" / "pocket_definition.json"
DEFAULT_REFERENCES = {
    "EXP_5YWY_inactive": PROJECT / "data" / "raw" / "structures" / "5YWY.cif",
    "EXP_7D7M_active": PROJECT / "data" / "raw" / "structures" / "7D7M.cif",
}


def kabsch_rmsd(first: np.ndarray, second: np.ndarray) -> float:
    if first.shape != second.shape or first.ndim != 2 or first.shape[1] != 3:
        raise ValueError("Kabsch inputs must be matching N x 3 matrices")
    centered_first = first - first.mean(axis=0)
    centered_second = second - second.mean(axis=0)
    u, _, vt = np.linalg.svd(centered_first.T @ centered_second)
    rotation = u @ vt
    if np.linalg.det(rotation) < 0:
        u[:, -1] *= -1
        rotation = u @ vt
    aligned = centered_first @ rotation
    return float(np.sqrt(np.mean(np.sum((aligned - centered_second) ** 2, axis=1))))


def load_coordinate_maps(path: Path) -> tuple[dict[int, np.ndarray], dict[int, np.ndarray]]:
    if path.suffix.lower() == ".pdb":
        atoms = parse_pdb_atoms(path)
    elif path.suffix.lower() in {".cif", ".mmcif"}:
        atoms = parse_mmcif_atoms(path, label_chain="A")
    else:
        raise ValueError(f"Unsupported structure format: {path}")
    lookup = atom_lookup(atoms)
    ca: dict[int, np.ndarray] = {}
    sidechain: dict[int, np.ndarray] = {}
    residue_numbers = {atom.residue_number for atom in atoms}
    for position in residue_numbers:
        ca_atom = lookup.get((position, "CA"))
        if ca_atom is None:
            continue
        ca[position] = np.asarray(ca_atom.xyz, dtype=float)
        representative = lookup.get((position, "CB"), ca_atom)
        sidechain[position] = np.asarray(representative.xyz, dtype=float)
    return ca, sidechain


def geometry_features(
    ca: dict[int, np.ndarray],
    sidechain: dict[int, np.ndarray],
    positions: list[int],
) -> tuple[dict[str, float], np.ndarray]:
    missing = [position for position in positions if position not in ca or position not in sidechain]
    if missing:
        raise ValueError("missing feature-panel positions: " + ",".join(map(str, missing)))
    ca_matrix = np.vstack([ca[position] for position in positions])
    sidechain_matrix = np.vstack([sidechain[position] for position in positions])
    features: dict[str, float] = {}
    pair_distances: list[float] = []
    for first_index, second_index in combinations(range(len(positions)), 2):
        value = float(np.linalg.norm(ca_matrix[first_index] - ca_matrix[second_index]))
        pair_distances.append(value)
        features[f"ca_dist_{positions[first_index]}_{positions[second_index]}"] = value
    centroid = ca_matrix.mean(axis=0)
    for index, position in enumerate(positions):
        features[f"sidechain_radius_{position}"] = float(np.linalg.norm(sidechain_matrix[index] - centroid))
    pair_array = np.asarray(pair_distances, dtype=float)
    ca_centered = ca_matrix - ca_matrix.mean(axis=0)
    side_centered = sidechain_matrix - sidechain_matrix.mean(axis=0)
    features.update({
        "summary_ca_radius_gyration": float(np.sqrt(np.mean(np.sum(ca_centered**2, axis=1)))),
        "summary_sidechain_radius_gyration": float(np.sqrt(np.mean(np.sum(side_centered**2, axis=1)))),
        "summary_mean_ca_distance": float(pair_array.mean()),
        "summary_contact_fraction_8A": float(np.mean(pair_array <= 8.0)),
        "summary_min_ca_distance": float(pair_array.min()),
        "summary_max_ca_distance": float(pair_array.max()),
    })
    return features, ca_matrix


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    required = {"model_id", "ae_structure_path", "source_type", "qc_status"}
    if not rows:
        raise ValueError(f"{path}: conformer manifest is empty")
    missing = required - set(rows[0])
    if missing:
        raise ValueError(f"{path}: missing columns {sorted(missing)}")
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            formatted = {key: (f"{value:.6f}" if isinstance(value, float) else value) for key, value in row.items()}
            writer.writerow(formatted)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Extract deterministic 196-dimensional EP4 pocket geometry features")
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--pocket-definition", type=Path, default=DEFAULT_POCKET)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--include-experimental-anchors", action=argparse.BooleanOptionalAction, default=True)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    definition = json.loads(args.pocket_definition.read_text(encoding="utf-8"))
    positions = pocket_positions(definition, panel_only=True)
    if len(positions) != 19:
        raise ValueError(f"Expected a frozen 19-residue feature panel, found {len(positions)}")
    manifest_rows = read_manifest(args.manifest)
    entries: list[dict[str, str]] = []
    if args.include_experimental_anchors:
        entries.extend([
            {
                "model_id": model_id,
                "ae_structure_path": project_path(path, PROJECT),
                "source_type": "experimental_reference",
                "qc_status": "REFERENCE",
                "run_label": "experimental_anchor",
                "query_id": model_id,
                "rank": "",
                "model_number": "",
                "seed": "",
                "pocket_mean_plddt": "",
                "pocket_mean_pae": "",
            }
            for model_id, path in DEFAULT_REFERENCES.items()
        ])
    entries.extend(manifest_rows)

    coordinate_cache: dict[str, tuple[dict[int, np.ndarray], dict[int, np.ndarray]]] = {}
    anchor_ca: dict[str, np.ndarray] = {}
    for model_id, path in DEFAULT_REFERENCES.items():
        ca, side = load_coordinate_maps(path)
        _, matrix = geometry_features(ca, side, positions)
        anchor_ca[model_id] = matrix

    feature_rows: list[dict[str, Any]] = []
    qc_rows: list[dict[str, Any]] = []
    feature_names: list[str] | None = None
    for entry in entries:
        structure_value = entry.get("ae_structure_path", "")
        if not structure_value:
            qc_rows.append({"model_id": entry.get("model_id", ""), "feature_status": "FAIL", "reason": "empty_structure_path"})
            continue
        structure = resolve_project_path(structure_value, PROJECT).resolve()
        if not structure.is_file():
            qc_rows.append({"model_id": entry.get("model_id", ""), "feature_status": "FAIL", "reason": f"missing_structure:{structure}"})
            continue
        try:
            ca, sidechain = load_coordinate_maps(structure)
            features, ca_matrix = geometry_features(ca, sidechain, positions)
        except (OSError, ValueError) as exc:
            qc_rows.append({"model_id": entry.get("model_id", ""), "feature_status": "FAIL", "reason": str(exc)})
            continue
        if feature_names is None:
            feature_names = list(features)
        elif list(features) != feature_names:
            raise RuntimeError("Feature order changed between conformers")
        rmsd_inactive = kabsch_rmsd(ca_matrix, anchor_ca["EXP_5YWY_inactive"])
        rmsd_active = kabsch_rmsd(ca_matrix, anchor_ca["EXP_7D7M_active"])
        metadata = {
            "model_id": entry["model_id"],
            "source_type": entry.get("source_type", ""),
            "structure_path": project_path(structure, PROJECT),
            "structure_sha256": sha256_file(structure),
        }
        feature_rows.append({**metadata, **features})
        qc_rows.append({
            **metadata,
            "feature_status": "PASS",
            "reason": "",
            "run_label": entry.get("run_label", ""),
            "query_id": entry.get("query_id", ""),
            "rank": entry.get("rank", ""),
            "model_number": entry.get("model_number", ""),
            "seed": entry.get("seed", ""),
            "input_qc_status": entry.get("qc_status", ""),
            "pocket_mean_plddt": entry.get("pocket_mean_plddt", ""),
            "pocket_mean_pae": entry.get("pocket_mean_pae", ""),
            "pocket_rmsd_to_5YWY_A": rmsd_inactive,
            "pocket_rmsd_to_7D7M_A": rmsd_active,
            "nearest_reference": "5YWY_inactive" if rmsd_inactive <= rmsd_active else "7D7M_active",
        })

    if not feature_rows or feature_names is None:
        print("ERROR: no conformers produced a complete feature vector", file=sys.stderr)
        return 3
    metadata_fields = ["model_id", "source_type", "structure_path", "structure_sha256"]
    write_csv(args.output_dir / "pocket_features.csv", feature_rows, metadata_fields + feature_names)
    qc_fields = [
        "model_id", "source_type", "structure_path", "structure_sha256", "feature_status", "reason",
        "run_label", "query_id", "rank", "model_number", "seed", "input_qc_status",
        "pocket_mean_plddt", "pocket_mean_pae", "pocket_rmsd_to_5YWY_A",
        "pocket_rmsd_to_7D7M_A", "nearest_reference",
    ]
    write_csv(args.output_dir / "feature_qc.csv", qc_rows, qc_fields)
    schema = {
        "schema_version": "1.0",
        "target": "PTGER4/EP4",
        "feature_panel_positions": positions,
        "feature_count": len(feature_names),
        "feature_names": feature_names,
        "model_input_columns": feature_names,
        "excluded_from_model_input": metadata_fields + ["pLDDT", "PAE", "reference labels", "source type"],
        "feature_groups": {
            "pairwise_ca_distances": len(list(combinations(positions, 2))),
            "sidechain_radial_distances": len(positions),
            "geometry_summaries": 6,
        },
        "standardization_rule": "Fit mean and standard deviation on the training split only; do not standardize before splitting.",
        "pocket_definition_sha256": sha256_file(args.pocket_definition),
        "manifest_sha256": sha256_file(args.manifest),
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "feature_schema.json").write_text(json.dumps(schema, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"Feature-complete conformers: {len(feature_rows)}")
    print(f"Feature dimension: {len(feature_names)}")
    print(f"Output: {args.output_dir.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
