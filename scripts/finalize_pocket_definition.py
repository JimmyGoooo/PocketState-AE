#!/usr/bin/env python3
"""Freeze the experimentally anchored PTGER4 pocket definition.

The input is the combined 6 A contact table produced by
``extract_c3_ligand_pockets.py``.  This command validates that table, records
its SHA-256 digest, and writes a deterministic residue table plus a JSON
definition for downstream PocketState-AE stages.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = (
    PROJECT_ROOT
    / "data"
    / "derived"
    / "c3_ligand_pockets"
    / "ep4_ligand_pocket_residues_within_6A.csv"
)
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "data" / "processed" / "pocket_definition"

SCHEMA_VERSION = "1.0"
EXPECTED_CUTOFF_ANGSTROM = 6.0
EXPECTED_POSITION_COUNT = 38
EXPECTED_FEATURE_PANEL_COUNT = 19

ANCHORS = (
    ("5YWY_7UR", "5YWY", "7UR"),
    ("7D7M_P2E", "7D7M", "P2E"),
)
ANCHOR_BY_PAIR = {(pdb_id, ligand_id): anchor_id for anchor_id, pdb_id, ligand_id in ANCHORS}
REQUIRED_COLUMNS = {
    "pdb_id",
    "ligand_id",
    "cutoff_angstrom",
    "receptor_label_chain",
    "auth_seq_id",
    "residue_name_3letter",
    "minimum_heavy_atom_distance_angstrom",
}


@dataclass(frozen=True)
class PocketPosition:
    """One UniProt position and its observations in the two anchors."""

    uniprot_position: int
    residue_name_3letter: str
    distances: dict[str, float]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _parse_positive_int(value: str, *, row_number: int, field: str) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Row {row_number}: {field} must be an integer") from exc
    if parsed < 1:
        raise ValueError(f"Row {row_number}: {field} must be positive")
    return parsed


def _parse_finite_float(value: str, *, row_number: int, field: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Row {row_number}: {field} must be numeric") from exc
    if not math.isfinite(parsed):
        raise ValueError(f"Row {row_number}: {field} must be finite")
    return parsed


def load_pocket_positions(input_csv: Path) -> list[PocketPosition]:
    """Load and strictly validate the combined anchor-contact table."""

    if not input_csv.is_file():
        raise FileNotFoundError(f"Pocket contact table does not exist: {input_csv}")

    residues: dict[int, str] = {}
    observations: dict[tuple[int, str], float] = {}
    seen_cutoffs: set[float] = set()

    with input_csv.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError("Pocket contact table has no header")
        missing = sorted(REQUIRED_COLUMNS - set(reader.fieldnames))
        if missing:
            raise ValueError(f"Pocket contact table is missing columns: {', '.join(missing)}")

        for row_number, row in enumerate(reader, start=2):
            pdb_id = (row.get("pdb_id") or "").strip().upper()
            ligand_id = (row.get("ligand_id") or "").strip().upper()
            anchor_id = ANCHOR_BY_PAIR.get((pdb_id, ligand_id))
            if anchor_id is None:
                raise ValueError(
                    f"Row {row_number}: unexpected anchor pair {pdb_id or '<blank>'}/"
                    f"{ligand_id or '<blank>'}"
                )
            chain = (row.get("receptor_label_chain") or "").strip()
            if chain != "A":
                raise ValueError(f"Row {row_number}: expected receptor label chain A, found {chain!r}")

            cutoff = _parse_finite_float(
                row.get("cutoff_angstrom", ""), row_number=row_number, field="cutoff_angstrom"
            )
            seen_cutoffs.add(cutoff)
            if not math.isclose(cutoff, EXPECTED_CUTOFF_ANGSTROM, rel_tol=0.0, abs_tol=1e-9):
                raise ValueError(
                    f"Row {row_number}: expected a {EXPECTED_CUTOFF_ANGSTROM:g} A cutoff, found {cutoff:g}"
                )

            position = _parse_positive_int(
                row.get("auth_seq_id", ""), row_number=row_number, field="auth_seq_id"
            )
            residue_name = (row.get("residue_name_3letter") or "").strip().upper()
            if len(residue_name) != 3 or not residue_name.isalpha():
                raise ValueError(
                    f"Row {row_number}: residue_name_3letter must contain three letters"
                )
            prior_residue = residues.get(position)
            if prior_residue is not None and prior_residue != residue_name:
                raise ValueError(
                    f"Row {row_number}: conflicting residue names at position {position}: "
                    f"{prior_residue} versus {residue_name}"
                )
            residues[position] = residue_name

            distance = _parse_finite_float(
                row.get("minimum_heavy_atom_distance_angstrom", ""),
                row_number=row_number,
                field="minimum_heavy_atom_distance_angstrom",
            )
            if distance < 0 or distance > cutoff + 1e-9:
                raise ValueError(
                    f"Row {row_number}: distance {distance:g} lies outside the 0-{cutoff:g} A cutoff"
                )
            key = (position, anchor_id)
            if key in observations:
                prior_distance = observations[key]
                if not math.isclose(prior_distance, distance, rel_tol=0.0, abs_tol=1e-12):
                    raise ValueError(
                        f"Row {row_number}: conflicting {anchor_id} distances at position "
                        f"{position}: {prior_distance:g} versus {distance:g}"
                    )
                raise ValueError(
                    f"Row {row_number}: duplicate {anchor_id} observation at position {position}"
                )
            observations[key] = distance

    if not observations:
        raise ValueError("Pocket contact table contains no observations")
    if len(seen_cutoffs) != 1:
        raise ValueError("Pocket contact table contains conflicting cutoffs")

    positions = [
        PocketPosition(
            uniprot_position=position,
            residue_name_3letter=residues[position],
            distances={
                anchor_id: observations[(position, anchor_id)]
                for anchor_id, _, _ in ANCHORS
                if (position, anchor_id) in observations
            },
        )
        for position in sorted(residues)
    ]
    if len(positions) != EXPECTED_POSITION_COUNT:
        raise ValueError(
            f"Expected {EXPECTED_POSITION_COUNT} unique pocket positions, found {len(positions)}"
        )
    return positions


def select_feature_panel(positions: Sequence[PocketPosition]) -> dict[int, str]:
    """Select shared sites plus the two nearest sites specific to each anchor."""

    anchor_ids = [anchor_id for anchor_id, _, _ in ANCHORS]
    shared = [position for position in positions if all(a in position.distances for a in anchor_ids)]
    selected: dict[int, str] = {
        position.uniprot_position: "shared_between_anchors" for position in shared
    }

    for anchor_id in anchor_ids:
        other_anchor = next(candidate for candidate in anchor_ids if candidate != anchor_id)
        specific = [
            position
            for position in positions
            if anchor_id in position.distances and other_anchor not in position.distances
        ]
        specific.sort(key=lambda position: (position.distances[anchor_id], position.uniprot_position))
        if len(specific) < 2:
            raise ValueError(f"Anchor {anchor_id} has fewer than two anchor-specific positions")
        for position in specific[:2]:
            selected[position.uniprot_position] = f"nearest_specific_to_{anchor_id}"

    if len(selected) != EXPECTED_FEATURE_PANEL_COUNT:
        raise ValueError(
            f"Expected a {EXPECTED_FEATURE_PANEL_COUNT}-position feature panel, found {len(selected)}"
        )
    return selected


def _format_distance(value: float | None) -> str:
    return "" if value is None else format(value, ".15g")


def _source_path(path: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(PROJECT_ROOT.resolve()).as_posix()
    except ValueError:
        return resolved.as_posix()


def build_outputs(input_csv: Path) -> tuple[str, str]:
    """Return deterministic CSV and JSON documents without writing files."""

    positions = load_pocket_positions(input_csv)
    feature_panel = select_feature_panel(positions)
    input_sha256 = sha256_file(input_csv)
    anchor_ids = [anchor_id for anchor_id, _, _ in ANCHORS]

    csv_fields = [
        "order",
        "uniprot_position",
        "residue_name_3letter",
        "in_5YWY_7UR_pocket",
        "distance_to_7UR_angstrom",
        "in_7D7M_P2E_pocket",
        "distance_to_P2E_angstrom",
        "in_feature_panel",
        "feature_panel_reason",
    ]
    csv_buffer = io.StringIO(newline="")
    writer = csv.DictWriter(csv_buffer, fieldnames=csv_fields, lineterminator="\n")
    writer.writeheader()

    json_positions = []
    for order, position in enumerate(positions, start=1):
        anchor_values = {
            anchor_id: {
                "present": anchor_id in position.distances,
                "minimum_heavy_atom_distance_angstrom": position.distances.get(anchor_id),
            }
            for anchor_id in anchor_ids
        }
        reason = feature_panel.get(position.uniprot_position)
        writer.writerow(
            {
                "order": order,
                "uniprot_position": position.uniprot_position,
                "residue_name_3letter": position.residue_name_3letter,
                "in_5YWY_7UR_pocket": "YES" if "5YWY_7UR" in position.distances else "NO",
                "distance_to_7UR_angstrom": _format_distance(
                    position.distances.get("5YWY_7UR")
                ),
                "in_7D7M_P2E_pocket": "YES" if "7D7M_P2E" in position.distances else "NO",
                "distance_to_P2E_angstrom": _format_distance(
                    position.distances.get("7D7M_P2E")
                ),
                "in_feature_panel": "YES" if reason is not None else "NO",
                "feature_panel_reason": reason or "",
            }
        )
        json_positions.append(
            {
                "order": order,
                "uniprot_position": position.uniprot_position,
                "residue_name_3letter": position.residue_name_3letter,
                "anchors": anchor_values,
                "in_feature_panel": reason is not None,
                "feature_panel_reason": reason,
            }
        )

    document = {
        "schema_version": SCHEMA_VERSION,
        "target": {
            "gene": "PTGER4",
            "protein": "prostaglandin E2 receptor EP4 subtype",
            "uniprot_id": "P35408",
        },
        "cutoff_angstrom": EXPECTED_CUTOFF_ANGSTROM,
        "sources": {
            "input_csv": {
                "path": _source_path(input_csv),
                "sha256": input_sha256,
            },
            "anchors": [
                {
                    "anchor_id": anchor_id,
                    "pdb_id": pdb_id,
                    "ligand_id": ligand_id,
                    "receptor_label_chain": "A",
                }
                for anchor_id, pdb_id, ligand_id in ANCHORS
            ],
        },
        "position_count": len(json_positions),
        "positions": json_positions,
        "feature_panel": {
            "selection_rule": (
                "all positions present in both anchors plus the two anchor-specific "
                "positions with the shortest ligand heavy-atom distance for each anchor; "
                "ties are resolved by ascending UniProt position"
            ),
            "position_count": len(feature_panel),
            "positions": sorted(feature_panel),
        },
    }
    json_text = json.dumps(document, ensure_ascii=False, indent=2) + "\n"
    return csv_buffer.getvalue(), json_text


def finalize_pocket_definition(input_csv: Path, output_dir: Path) -> tuple[Path, Path]:
    """Validate, freeze, and atomically write the pocket definition."""

    csv_text, json_text = build_outputs(input_csv)
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / "pocket_residues.csv"
    json_path = output_dir / "pocket_definition.json"
    csv_tmp = output_dir / ".pocket_residues.csv.tmp"
    json_tmp = output_dir / ".pocket_definition.json.tmp"
    try:
        csv_tmp.write_text(csv_text, encoding="utf-8", newline="")
        json_tmp.write_text(json_text, encoding="utf-8", newline="")
        csv_tmp.replace(csv_path)
        json_tmp.replace(json_path)
    finally:
        if csv_tmp.exists():
            csv_tmp.unlink()
        if json_tmp.exists():
            json_tmp.unlink()
    return csv_path, json_path


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Freeze the validated PTGER4 6 A pocket definition for PocketState-AE."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help=f"Combined 6 A pocket CSV (default: {DEFAULT_INPUT})",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Destination directory (default: {DEFAULT_OUTPUT_DIR})",
    )
    args = parser.parse_args(argv)
    csv_path, json_path = finalize_pocket_definition(args.input, args.output_dir)
    print(f"Wrote {csv_path}")
    print(f"Wrote {json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
