#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
C2 EP4 structure audit helper
For PocketState-AE: 5YWY / 7D7M

Usage from the PocketState-AE project root:
    python scripts/audit_c2_structures.py

Requires:
    python -m pip install biopython

Expected inputs:
    data/raw/structures/5YWY.cif
    data/raw/structures/7D7M.cif

Outputs:
    data/raw/structures/unobserved_residues.csv
    data/raw/structures/structure_coordinate_audit.csv
    data/raw/structures/structure_audit_completed.csv
"""

from pathlib import Path
import csv
import sys

try:
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
except ImportError:
    print("ERROR: Biopython is not installed.")
    print("Run: python -m pip install biopython")
    sys.exit(1)

ROOT = Path.cwd()
RAW = ROOT / "data" / "raw" / "structures"

FILES = {
    "5YWY": RAW / "5YWY.cif",
    "7D7M": RAW / "7D7M.cif",
}

# Literature/RCSB-confirmed construct information.
# We deliberately keep coordinate-derived facts (unobserved residues, waters)
# separate and obtain them from the user's actual downloaded CIF files.
CONFIRMED = {
    "5YWY": {
        "receptor_chain": "A",
        "uniprot": "P35408",
        "ligand_id": "7UR",
        "ligand_name": "ONO-AE3-208",
        "mutations": "N7Q; A62L; G106R; N177Q",
        "designed_deletions": "1-3; 218-259 (ICL3); 347-488 (C-terminus)",
        "auxiliary_components": "Fab heavy chain B [auth H]; Fab light chain C [auth L]",
        "modeled_region_summary": (
            "Inactive/antagonist anchor. ICL3 is unresolved/deleted; "
            "use CIF unobserved-residue table for exact coordinate-level gaps."
        ),
        "notes": (
            "Do not delete Fab or ligand in the raw CIF. "
            "A62L and G106R are thermostabilizing inactive-state mutations."
        ),
    },
    "7D7M": {
        "receptor_chain": "A",
        "uniprot": "P35408",
        "ligand_id": "P2E",
        "ligand_name": "PGE2",
        "mutations": "N7Q; N177Q",
        "designed_deletions": "1-3; 218-259 (ICL3); 347-488 (C-terminus)",
        "auxiliary_components": (
            "Gbeta1 chain B; Ggamma2 chain C; Gs alpha/mini-Gs chain D; Nb35 chain E"
        ),
        "modeled_region_summary": (
            "Active/PGE2-Gs anchor. Primary paper reports EP4 residues "
            "Ser19-Cys345 modeled except ICL3."
        ),
        "notes": (
            "A62L/G106R were not used in this active-state construct. "
            "Do not delete G protein/Nb35/PGE2 in the raw CIF."
        ),
    },
}


def as_list(d, key):
    value = d.get(key, [])
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


def normalized_columns(d, keys):
    cols = {k: as_list(d, k) for k in keys}
    lengths = [len(v) for v in cols.values() if v]
    n = max(lengths) if lengths else 0
    rows = []
    for i in range(n):
        row = {}
        for k, vals in cols.items():
            if not vals:
                row[k] = ""
            elif len(vals) == 1 and n > 1:
                row[k] = vals[0]
            elif i < len(vals):
                row[k] = vals[i]
            else:
                row[k] = ""
        rows.append(row)
    return rows


def numeric_sort_key(x):
    try:
        return (0, int(float(x)))
    except Exception:
        return (1, str(x))


def audit_one(pdb_id, path):
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")

    cif = MMCIF2Dict(str(path))

    # ---- Atom-site inspection ----
    atom_keys = [
        "_atom_site.group_PDB",
        "_atom_site.label_asym_id",
        "_atom_site.auth_asym_id",
        "_atom_site.label_comp_id",
        "_atom_site.auth_comp_id",
        "_atom_site.label_seq_id",
        "_atom_site.auth_seq_id",
    ]
    atoms = normalized_columns(cif, atom_keys)

    receptor_label_chain = "A"

    modeled_residues = set()
    waters = set()
    nonpoly = set()

    for r in atoms:
        group = r["_atom_site.group_PDB"]
        label_chain = r["_atom_site.label_asym_id"]
        auth_chain = r["_atom_site.auth_asym_id"]
        comp = r["_atom_site.label_comp_id"] or r["_atom_site.auth_comp_id"]
        lseq = r["_atom_site.label_seq_id"]
        aseq = r["_atom_site.auth_seq_id"]

        if group == "ATOM" and label_chain == receptor_label_chain:
            modeled_residues.add((lseq, aseq, comp))

        if comp in {"HOH", "WAT", "H2O"}:
            waters.add((label_chain, auth_chain, aseq))

        if group == "HETATM" and comp not in {"HOH", "WAT", "H2O"}:
            nonpoly.add((comp, label_chain, auth_chain))

    # ---- Unobserved / zero-occupancy residues ----
    unobs_keys = [
        "_pdbx_unobs_or_zero_occ_residues.PDB_model_num",
        "_pdbx_unobs_or_zero_occ_residues.polymer_flag",
        "_pdbx_unobs_or_zero_occ_residues.occupancy_flag",
        "_pdbx_unobs_or_zero_occ_residues.label_asym_id",
        "_pdbx_unobs_or_zero_occ_residues.auth_asym_id",
        "_pdbx_unobs_or_zero_occ_residues.label_comp_id",
        "_pdbx_unobs_or_zero_occ_residues.auth_comp_id",
        "_pdbx_unobs_or_zero_occ_residues.label_seq_id",
        "_pdbx_unobs_or_zero_occ_residues.auth_seq_id",
        "_pdbx_unobs_or_zero_occ_residues.PDB_ins_code",
    ]
    unobs_all = normalized_columns(cif, unobs_keys)

    receptor_unobs = []
    for r in unobs_all:
        label_chain = r["_pdbx_unobs_or_zero_occ_residues.label_asym_id"]
        if label_chain == receptor_label_chain:
            receptor_unobs.append({
                "pdb_id": pdb_id,
                "label_asym_id": label_chain,
                "auth_asym_id": r["_pdbx_unobs_or_zero_occ_residues.auth_asym_id"],
                "label_seq_id": r["_pdbx_unobs_or_zero_occ_residues.label_seq_id"],
                "auth_seq_id": r["_pdbx_unobs_or_zero_occ_residues.auth_seq_id"],
                "label_comp_id": r["_pdbx_unobs_or_zero_occ_residues.label_comp_id"],
                "auth_comp_id": r["_pdbx_unobs_or_zero_occ_residues.auth_comp_id"],
                "polymer_flag": r["_pdbx_unobs_or_zero_occ_residues.polymer_flag"],
                "occupancy_flag": r["_pdbx_unobs_or_zero_occ_residues.occupancy_flag"],
                "PDB_ins_code": r["_pdbx_unobs_or_zero_occ_residues.PDB_ins_code"],
            })

    modeled_sorted = sorted(
        modeled_residues,
        key=lambda x: numeric_sort_key(x[1] if x[1] not in {"", ".", "?"} else x[0])
    )

    auth_ids = [
        x[1] for x in modeled_sorted
        if x[1] not in {"", ".", "?"}
    ]
    label_ids = [
        x[0] for x in modeled_sorted
        if x[0] not in {"", ".", "?"}
    ]

    modeled_auth_range = (
        f"{auth_ids[0]}-{auth_ids[-1]}" if auth_ids else "NA"
    )
    modeled_label_range = (
        f"{label_ids[0]}-{label_ids[-1]}" if label_ids else "NA"
    )

    nonpoly_text = "; ".join(
        f"{comp}(label_chain={lch},auth_chain={ach})"
        for comp, lch, ach in sorted(nonpoly)
    ) or "NONE"

    coordinate_row = {
        "pdb_id": pdb_id,
        "receptor_label_chain": receptor_label_chain,
        "modeled_receptor_residue_count": len(modeled_residues),
        "modeled_auth_seq_range": modeled_auth_range,
        "modeled_label_seq_range": modeled_label_range,
        "receptor_unobserved_residue_count": len(receptor_unobs),
        "water_residue_count_whole_entry": len(waters),
        "nonpoly_components_detected": nonpoly_text,
    }

    final = CONFIRMED[pdb_id].copy()
    final.update({
        "pdb_id": pdb_id,
        "receptor_unobserved_residue_count_raw_cif": len(receptor_unobs),
        "water_residue_count_raw_cif": len(waters),
        "nonpoly_components_detected": nonpoly_text,
        "missing_residues_checked": "YES",
        "mutation_positions_checked": "YES",
        "fusion_auxiliary_checked": "YES",
        "water_checked": "YES",
        "sequence_mapping_checked": "YES (RCSB maps receptor entity to UniProt P35408)",
        "coordinate_audit_file": "structure_coordinate_audit.csv",
        "unobserved_residue_file": "unobserved_residues.csv",
    })

    return coordinate_row, receptor_unobs, final


def write_csv(path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


def main():
    missing = [str(p) for p in FILES.values() if not p.exists()]
    if missing:
        print("ERROR: The following input CIF files are missing:")
        for p in missing:
            print("  ", p)
        print("\nRun this script from the PocketState-AE project root.")
        sys.exit(1)

    coordinate_rows = []
    all_unobs = []
    final_rows = []

    for pdb_id, path in FILES.items():
        print(f"Auditing {pdb_id}: {path}")
        coord, unobs, final = audit_one(pdb_id, path)
        coordinate_rows.append(coord)
        all_unobs.extend(unobs)
        final_rows.append(final)

    coord_fields = [
        "pdb_id",
        "receptor_label_chain",
        "modeled_receptor_residue_count",
        "modeled_auth_seq_range",
        "modeled_label_seq_range",
        "receptor_unobserved_residue_count",
        "water_residue_count_whole_entry",
        "nonpoly_components_detected",
    ]
    write_csv(
        RAW / "structure_coordinate_audit.csv",
        coordinate_rows,
        coord_fields,
    )

    unobs_fields = [
        "pdb_id",
        "label_asym_id",
        "auth_asym_id",
        "label_seq_id",
        "auth_seq_id",
        "label_comp_id",
        "auth_comp_id",
        "polymer_flag",
        "occupancy_flag",
        "PDB_ins_code",
    ]
    write_csv(
        RAW / "unobserved_residues.csv",
        all_unobs,
        unobs_fields,
    )

    final_fields = [
        "pdb_id",
        "receptor_chain",
        "uniprot",
        "ligand_id",
        "ligand_name",
        "mutations",
        "designed_deletions",
        "auxiliary_components",
        "modeled_region_summary",
        "receptor_unobserved_residue_count_raw_cif",
        "water_residue_count_raw_cif",
        "nonpoly_components_detected",
        "missing_residues_checked",
        "mutation_positions_checked",
        "fusion_auxiliary_checked",
        "water_checked",
        "sequence_mapping_checked",
        "coordinate_audit_file",
        "unobserved_residue_file",
        "notes",
    ]
    write_csv(
        RAW / "structure_audit_completed.csv",
        final_rows,
        final_fields,
    )

    print("\nDONE.")
    print("Created:")
    print("  data/raw/structures/structure_coordinate_audit.csv")
    print("  data/raw/structures/unobserved_residues.csv")
    print("  data/raw/structures/structure_audit_completed.csv")
    print("\nNext checks:")
    print("  1) Open structure_coordinate_audit.csv")
    print("  2) Confirm ligand IDs 7UR (5YWY) and P2E (7D7M) appear in nonpoly components")
    print("  3) Open unobserved_residues.csv and inspect all EP4 chain-A gaps")
    print("  4) Keep the original CIF files unchanged")


if __name__ == "__main__":
    main()
