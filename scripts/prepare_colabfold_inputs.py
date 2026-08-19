#!/usr/bin/env python3
"""Prepare validated, provenance-preserving ColabFold FASTA inputs for PTGER4.

The script always writes the full-length wild-type sequence.  It also creates
a candidate-site table from the experimental 6 Å pockets.  Variant sequences
are written only after the user explicitly fills selected_variants.csv.
"""

from __future__ import annotations

import csv
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
RAW_FASTA = PROJECT / "data" / "raw" / "structures" / "PTGER4.fasta"
POCKET_TABLE = PROJECT / "data" / "derived" / "c3_ligand_pockets" / "ep4_ligand_pocket_residues_within_6A.csv"
OUT = PROJECT / "data" / "derived" / "colabfold_inputs"
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


def read_single_fasta(path: Path) -> tuple[str, str]:
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(lines) < 2 or not lines[0].startswith(">"):
        raise ValueError(f"Expected one FASTA record in {path.name}")
    sequence = "".join(lines[1:]).upper()
    if not sequence or set(sequence) - VALID_AA:
        raise ValueError("Reference FASTA contains invalid amino-acid characters")
    return lines[0][1:], sequence


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def write_candidate_sites() -> None:
    by_site: dict[tuple[int, str], dict[str, str]] = {}
    with POCKET_TABLE.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            key = (int(row["auth_seq_id"]), row["residue_name_3letter"])
            current = by_site.setdefault(key, {
                "uniprot_position": str(key[0]),
                "residue_name_3letter": key[1],
                "in_5YWY_7UR_pocket": "NO",
                "distance_to_7UR_angstrom": "",
                "in_7D7M_P2E_pocket": "NO",
                "distance_to_P2E_angstrom": "",
            })
            if row["pdb_id"] == "5YWY":
                current["in_5YWY_7UR_pocket"] = "YES"
                current["distance_to_7UR_angstrom"] = row["minimum_heavy_atom_distance_angstrom"]
            elif row["pdb_id"] == "7D7M":
                current["in_7D7M_P2E_pocket"] = "YES"
                current["distance_to_P2E_angstrom"] = row["minimum_heavy_atom_distance_angstrom"]
    destination = OUT / "experimental_pocket_candidate_sites.csv"
    with destination.open("w", encoding="utf-8", newline="") as handle:
        fields = list(next(iter(by_site.values())).keys())
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(by_site[key] for key in sorted(by_site))


def create_variant_template(path: Path) -> None:
    if path.exists():
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            "variant_id", "wild_type_aa", "position", "mutant_aa", "rationale", "include_in_colabfold"
        ])
        writer.writeheader()


def read_selected_variants(path: Path, wild_type_sequence: str) -> list[tuple[str, str]]:
    variants = []
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            if row.get("include_in_colabfold", "").strip().upper() != "YES":
                continue
            variant_id = row["variant_id"].strip()
            wt = row["wild_type_aa"].strip().upper()
            mutant = row["mutant_aa"].strip().upper()
            position = int(row["position"])
            if not variant_id or len(wt) != 1 or len(mutant) != 1 or wt not in VALID_AA or mutant not in VALID_AA:
                raise ValueError(f"Invalid variant definition for {variant_id or 'unnamed row'}")
            if not 1 <= position <= len(wild_type_sequence):
                raise ValueError(f"{variant_id}: position {position} lies outside PTGER4")
            if wild_type_sequence[position - 1] != wt:
                raise ValueError(f"{variant_id}: expected {wt}{position}, but PTGER4 has {wild_type_sequence[position - 1]}{position}")
            sequence = wild_type_sequence[:position - 1] + mutant + wild_type_sequence[position:]
            variants.append((variant_id, sequence))
    return variants


def write_readme(reference_id: str, length: int) -> None:
    (OUT / "README.md").write_text(
        "# ColabFold input package\n\n"
        f"- Reference: `{reference_id}`\n"
        f"- Wild-type sequence length: {length} aa\n"
        "- Model objective: compare full-length PTGER4 WT with explicitly selected single-site variants.\n"
        "- Experimental pocket evidence: `experimental_pocket_candidate_sites.csv`, derived from 6 Å heavy-atom contacts to 7UR (5YWY) and P2E (7D7M).\n\n"
        "## Files\n\n"
        "- `PTGER4_WT.fasta`: validated WT input for ColabFold.\n"
        "- `PTGER4_WT_and_variant_panel.fasta`: WT plus explicitly selected variants, when present.\n"
        "- `experimental_pocket_candidate_sites.csv`: candidate locations, not recommended mutations.\n"
        "- `selected_variants.csv`: user decision sheet. Set `include_in_colabfold` to `YES` only after recording a rationale.\n\n"
        "## Generate a variant panel\n\n"
        "1. Fill one or more rows in `selected_variants.csv`.\n"
        "2. Run `python scripts/prepare_colabfold_inputs.py` from the project root.\n"
        "3. Submit either `PTGER4_WT.fasta` or `PTGER4_variant_panel.fasta` to ColabFold.\n\n"
        "ColabFold predictions are sequence-based; preserve 5YWY and 7D7M as the evidence for ligand-bound inactive/active pockets.\n",
        encoding="utf-8",
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    reference_id, wild_type_sequence = read_single_fasta(RAW_FASTA)
    write_fasta(OUT / "PTGER4_WT.fasta", [("PTGER4_WT_P35408_full_length", wild_type_sequence)])
    write_candidate_sites()
    selection_sheet = OUT / "selected_variants.csv"
    create_variant_template(selection_sheet)
    variants = read_selected_variants(selection_sheet, wild_type_sequence)
    if variants:
        write_fasta(OUT / "PTGER4_variant_panel.fasta", variants)
        write_fasta(
            OUT / "PTGER4_WT_and_variant_panel.fasta",
            [("PTGER4_WT_P35408_full_length", wild_type_sequence), *variants],
        )
        print(f"Wrote {len(variants)} selected variant sequences")
    else:
        print("No selected variants: wrote WT input and decision template only")
    write_readme(reference_id, len(wild_type_sequence))
    print(f"Prepared ColabFold inputs in {OUT.relative_to(PROJECT)}")


if __name__ == "__main__":
    main()
