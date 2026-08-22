#!/usr/bin/env python3
"""Shared, dependency-light helpers for the PocketState-AE P1 pipeline."""

from __future__ import annotations

import hashlib
import json
import math
import re
import shlex
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Iterator


VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")
SCORE_RE = re.compile(r"^(?P<query>.+)_scores_rank_(?P<rank>\d+)_(?P<tag>.+)\.json$", re.IGNORECASE)
MODEL_SEED_RE = re.compile(r"(?:^|_)model_(?P<model>\d+)(?:_|$).*?(?:^|_)seed_(?P<seed>\d+)(?:_|$)", re.IGNORECASE)


@dataclass(frozen=True)
class Atom:
    residue_number: int
    residue_name: str
    atom_name: str
    xyz: tuple[float, float, float]
    bfactor: float | None


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    identifier: str | None = None
    sequence: list[str] = []
    for line_number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if identifier is not None:
                records.append((identifier, "".join(sequence)))
            identifier = line[1:].split()[0]
            sequence = []
            if not identifier:
                raise ValueError(f"{path}: line {line_number} has an empty FASTA identifier")
        elif identifier is None:
            raise ValueError(f"{path}: line {line_number} contains sequence before an identifier")
        else:
            sequence.append(line.upper())
    if identifier is not None:
        records.append((identifier, "".join(sequence)))
    if not records:
        raise ValueError(f"{path}: no FASTA records")
    seen: set[str] = set()
    for name, seq in records:
        if name in seen:
            raise ValueError(f"{path}: duplicate FASTA identifier {name}")
        if not seq or set(seq) - VALID_AA:
            raise ValueError(f"{path}: {name} has an empty or non-standard sequence")
        seen.add(name)
    return records


def safe_identifier(identifier: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", identifier).strip("._")


def parse_score_name(path: Path) -> tuple[str, int, str, str, str] | None:
    match = SCORE_RE.match(path.name)
    if not match:
        return None
    query = match.group("query")
    rank = int(match.group("rank"))
    tag = match.group("tag")
    model = ""
    seed = ""
    model_seed = MODEL_SEED_RE.search(tag)
    if model_seed:
        model = model_seed.group("model")
        seed = model_seed.group("seed")
    return query, rank, tag, model, seed


def load_scores(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"{path}: score JSON is not an object")
    plddt = data.get("plddt")
    if not isinstance(plddt, list) or not plddt:
        raise ValueError(f"{path}: missing plddt array")
    values = []
    for index, value in enumerate(plddt, start=1):
        if not isinstance(value, (int, float)) or not math.isfinite(float(value)) or not 0 <= float(value) <= 100:
            raise ValueError(f"{path}: invalid pLDDT at position {index}")
        values.append(float(value))
    data["plddt"] = values
    pae = data.get("pae") or data.get("predicted_aligned_error")
    if pae is not None:
        if not isinstance(pae, list) or len(pae) != len(values):
            raise ValueError(f"{path}: PAE dimension does not match pLDDT")
        for row in pae:
            if not isinstance(row, list) or len(row) != len(values):
                raise ValueError(f"{path}: PAE is not a square residue matrix")
    return data


def parse_pdb_atoms(path: Path, chain: str | None = None) -> list[Atom]:
    atoms: list[Atom] = []
    seen: set[tuple[int, str]] = set()
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("ATOM  ") or len(line) < 54:
                continue
            atom_name = line[12:16].strip()
            altloc = line[16:17]
            chain_id = line[21:22].strip() or "A"
            if chain is not None and chain_id != chain:
                continue
            if altloc not in (" ", "A"):
                continue
            try:
                residue_number = int(line[22:26])
                xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue
            key = (residue_number, atom_name)
            if key in seen:
                continue
            seen.add(key)
            try:
                bfactor = float(line[60:66])
            except (ValueError, IndexError):
                bfactor = None
            atoms.append(Atom(residue_number, line[17:20].strip(), atom_name, xyz, bfactor))
    return atoms


def mmcif_atom_site_rows(path: Path) -> Iterator[dict[str, str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    index = 0
    while index < len(lines):
        if lines[index].strip() != "loop_":
            index += 1
            continue
        index += 1
        fields: list[str] = []
        while index < len(lines) and lines[index].startswith("_atom_site."):
            fields.append(lines[index].strip().split(".", 1)[1])
            index += 1
        if not fields:
            continue
        while index < len(lines):
            line = lines[index].strip()
            if not line or line.startswith("#") or line == "loop_" or line.startswith("_"):
                return
            values = shlex.split(line, posix=True)
            if len(values) != len(fields):
                raise ValueError(f"{path}: unsupported wrapped _atom_site row near line {index + 1}")
            yield dict(zip(fields, values))
            index += 1
        return
    raise ValueError(f"{path}: no _atom_site loop found")


def parse_mmcif_atoms(path: Path, label_chain: str = "A") -> list[Atom]:
    atoms: list[Atom] = []
    seen: set[tuple[int, str]] = set()
    for row in mmcif_atom_site_rows(path):
        if row.get("pdbx_PDB_model_num", "1") != "1" or row.get("group_PDB") != "ATOM":
            continue
        if row.get("label_asym_id") != label_chain:
            continue
        altloc = row.get("label_alt_id", ".")
        if altloc not in (".", "?", "A"):
            continue
        try:
            residue_number = int(row["auth_seq_id"])
            xyz = (float(row["Cartn_x"]), float(row["Cartn_y"]), float(row["Cartn_z"]))
        except (KeyError, ValueError):
            continue
        atom_name = row.get("auth_atom_id") or row.get("label_atom_id") or ""
        key = (residue_number, atom_name)
        if key in seen:
            continue
        seen.add(key)
        try:
            bfactor = float(row.get("B_iso_or_equiv", ""))
        except ValueError:
            bfactor = None
        atoms.append(Atom(residue_number, row.get("auth_comp_id") or row.get("label_comp_id", ""), atom_name, xyz, bfactor))
    return atoms


def atom_lookup(atoms: Iterable[Atom]) -> dict[tuple[int, str], Atom]:
    return {(atom.residue_number, atom.atom_name): atom for atom in atoms}


def distance(first: tuple[float, float, float], second: tuple[float, float, float]) -> float:
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(first, second)))


def mean(values: Iterable[float]) -> float | None:
    items = list(values)
    return sum(items) / len(items) if items else None


def project_path(path: Path, project: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(project.resolve()).as_posix()
    except ValueError:
        return str(resolved)


def resolve_project_path(value: str, project: Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else project / path


def pocket_positions(definition: dict[str, Any], panel_only: bool = False) -> list[int]:
    source: Any = definition.get("feature_panel") if panel_only else definition.get("positions")
    if isinstance(source, dict):
        source = source.get("positions") or source.get("selected_positions") or source.get("residues")
    if not isinstance(source, list):
        raise ValueError("Pocket definition has no usable positions list")
    positions: list[int] = []
    for item in source:
        if isinstance(item, int):
            positions.append(item)
        elif isinstance(item, str) and item.isdigit():
            positions.append(int(item))
        elif isinstance(item, dict):
            value = item.get("uniprot_position") or item.get("position")
            positions.append(int(value))
        else:
            raise ValueError(f"Invalid pocket position entry: {item!r}")
    if len(positions) != len(set(positions)):
        raise ValueError("Pocket definition contains duplicate positions")
    return positions
