#!/usr/bin/env python3
"""Extract EP4 chain-A residues within a ligand-centered distance cutoff.

This C3 helper reads the raw mmCIF files without modifying them and writes
derived, reproducible pocket tables.  It reports the closest heavy-atom
distance from each receptor residue to its cognate ligand.
"""

from __future__ import annotations

import csv
import math
import shlex
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
RAW = PROJECT / "data" / "raw" / "structures"
OUT = PROJECT / "data" / "derived" / "c3_ligand_pockets"
CUTOFF_ANGSTROM = 6.0

TARGETS = {
    "5YWY": {"ligand": "7UR", "ligand_label_chain": "D"},
    "7D7M": {"ligand": "P2E", "ligand_label_chain": "F"},
}


def atom_site_rows(cif_path: Path):
    """Yield dictionaries from the _atom_site loop of a standard mmCIF."""
    lines = cif_path.read_text(encoding="utf-8").splitlines()
    i = 0
    while i < len(lines):
        if lines[i].strip() != "loop_":
            i += 1
            continue
        i += 1
        fields = []
        while i < len(lines) and lines[i].startswith("_atom_site."):
            fields.append(lines[i].strip().split(".", 1)[1])
            i += 1
        if not fields:
            continue
        while i < len(lines):
            line = lines[i].strip()
            if not line or line.startswith("#") or line == "loop_" or line.startswith("_"):
                return
            values = shlex.split(line, posix=True)
            if len(values) != len(fields):
                raise ValueError(f"Unexpected atom_site row in {cif_path.name}: {line[:80]}")
            yield dict(zip(fields, values))
            i += 1
        return
    raise ValueError(f"No _atom_site loop found in {cif_path.name}")


def is_hydrogen(row: dict[str, str]) -> bool:
    return row["type_symbol"].upper() in {"H", "D"}


def coordinates(row: dict[str, str]) -> tuple[float, float, float]:
    return (float(row["Cartn_x"]), float(row["Cartn_y"]), float(row["Cartn_z"]))


def extract_target(pdb_id: str, target: dict[str, str]) -> list[dict[str, object]]:
    ligand_atoms = []
    receptor_residues: dict[tuple[str, str, str], list[dict[str, str]]] = {}
    for row in atom_site_rows(RAW / f"{pdb_id}.cif"):
        if row.get("pdbx_PDB_model_num") != "1" or is_hydrogen(row):
            continue
        if row["label_comp_id"] == target["ligand"] and row["label_asym_id"] == target["ligand_label_chain"]:
            ligand_atoms.append(row)
        if row["group_PDB"] == "ATOM" and row["label_asym_id"] == "A":
            key = (row["auth_seq_id"], row["label_seq_id"], row["auth_comp_id"])
            receptor_residues.setdefault(key, []).append(row)
    if not ligand_atoms:
        raise ValueError(f"{pdb_id}: no ligand atoms found for {target['ligand']}")

    cutoff_squared = CUTOFF_ANGSTROM**2
    results = []
    for (auth_seq_id, label_seq_id, residue_name), atoms in receptor_residues.items():
        best = None
        for receptor_atom in atoms:
            rx, ry, rz = coordinates(receptor_atom)
            for ligand_atom in ligand_atoms:
                lx, ly, lz = coordinates(ligand_atom)
                d2 = (rx - lx) ** 2 + (ry - ly) ** 2 + (rz - lz) ** 2
                if best is None or d2 < best[0]:
                    best = (d2, receptor_atom["auth_atom_id"], ligand_atom["label_atom_id"])
        assert best is not None
        if best[0] <= cutoff_squared:
            results.append({
                "pdb_id": pdb_id,
                "ligand_id": target["ligand"],
                "cutoff_angstrom": CUTOFF_ANGSTROM,
                "receptor_label_chain": "A",
                "auth_seq_id": int(auth_seq_id),
                "label_seq_id": int(label_seq_id),
                "residue_name_3letter": residue_name,
                "minimum_heavy_atom_distance_angstrom": round(math.sqrt(best[0]), 3),
                "closest_receptor_atom": best[1],
                "closest_ligand_atom": best[2],
            })
    return sorted(results, key=lambda item: item["auth_seq_id"])


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    all_rows = []
    for pdb_id, target in TARGETS.items():
        rows = extract_target(pdb_id, target)
        all_rows.extend(rows)
        destination = OUT / f"{pdb_id}_{target['ligand']}_ep4_chain_a_within_{CUTOFF_ANGSTROM:g}A.csv"
        with destination.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        print(f"{pdb_id}: {len(rows)} EP4 residues within {CUTOFF_ANGSTROM:g} A of {target['ligand']}")
    combined = OUT / f"ep4_ligand_pocket_residues_within_{CUTOFF_ANGSTROM:g}A.csv"
    with combined.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(all_rows[0]))
        writer.writeheader()
        writer.writerows(all_rows)
    print(f"Wrote {combined.relative_to(PROJECT)}")


if __name__ == "__main__":
    main()
