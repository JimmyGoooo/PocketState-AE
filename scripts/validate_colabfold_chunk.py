#!/usr/bin/env python3
"""Validate per-query artifacts from one ColabFold batch chunk.

The module intentionally depends only on the Python standard library so it can
run inside a lightweight Slurm post-processing job.  It also exposes
``read_fasta`` and ``validate_chunk`` for unit tests and downstream callers.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence


SCORE_RANK_ONE_RE = re.compile(
    r"(?:^|_)scores?_rank_001(?:_|\.|$)", re.IGNORECASE
)
STRUCTURE_RANK_ONE_RE = re.compile(
    r"(?:^|_)rank_001(?:_|\.|$)", re.IGNORECASE
)
PICKLE_SUFFIXES = {".pickle", ".pkl"}
STRUCTURE_SUFFIXES = {".pdb", ".cif", ".mmcif"}


@dataclass(frozen=True)
class FastaRecord:
    """One FASTA query, using the first whitespace-delimited header token."""

    identifier: str
    sequence: str


@dataclass
class ValidationResult:
    """Validation outcome for one FASTA identifier."""

    identifier: str
    reasons: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.reasons


def read_fasta(path: Path | str) -> list[FastaRecord]:
    """Read and minimally validate a multi-record protein FASTA file."""

    fasta_path = Path(path)
    records: list[FastaRecord] = []
    identifier: str | None = None
    sequence_parts: list[str] = []

    def finish_record() -> None:
        nonlocal identifier, sequence_parts
        if identifier is None:
            return
        sequence = "".join(sequence_parts).upper()
        if not sequence:
            raise ValueError(f"{identifier}: FASTA sequence is empty")
        if not re.fullmatch(r"[A-Z]+", sequence):
            raise ValueError(
                f"{identifier}: FASTA sequence contains non-letter characters"
            )
        records.append(FastaRecord(identifier, sequence))

    try:
        lines = fasta_path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise ValueError(f"cannot read FASTA {fasta_path}: {exc}") from exc

    for line_number, raw_line in enumerate(lines, start=1):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            finish_record()
            header = line[1:].strip()
            identifier = header.split(maxsplit=1)[0] if header else ""
            sequence_parts = []
            if not identifier:
                raise ValueError(f"line {line_number}: FASTA identifier is empty")
        elif identifier is None:
            raise ValueError(
                f"line {line_number}: sequence appears before a FASTA identifier"
            )
        else:
            sequence_parts.append(line)
    finish_record()

    if not records:
        raise ValueError("FASTA contains no records")
    identifiers = [record.identifier for record in records]
    duplicates = sorted(
        identifier for identifier in set(identifiers) if identifiers.count(identifier) > 1
    )
    if duplicates:
        raise ValueError(f"duplicate FASTA identifiers: {', '.join(duplicates)}")
    return records


def safe_identifier_aliases(identifier: str) -> tuple[str, ...]:
    """Return common ColabFold filename forms for a FASTA identifier."""

    aliases = [identifier]
    # Current ColabFold replaces filename-unsafe characters with underscores.
    aliases.append(re.sub(r"[^\w_.-]", "_", identifier, flags=re.UNICODE))
    # Some wrappers use an explicitly ASCII-safe equivalent.
    aliases.append(re.sub(r"[^A-Za-z0-9_.-]", "_", identifier))
    return tuple(dict.fromkeys(alias for alias in aliases if alias))


def _matching_alias_length(filename: str, aliases: Iterable[str]) -> int:
    """Return the longest alias forming the query prefix of ``filename``."""

    folded = filename.casefold()
    best = 0
    for alias in aliases:
        alias_folded = alias.casefold()
        if folded == alias_folded or folded.startswith(alias_folded + ".") or folded.startswith(alias_folded + "_"):
            best = max(best, len(alias_folded))
    return best


def _assign_files(
    records: Sequence[FastaRecord], files: Sequence[Path]
) -> dict[str, list[Path]]:
    """Assign each artifact to the most specific matching query prefix.

    Longest-prefix assignment prevents an identifier such as ``sample`` from
    claiming files that belong to ``sample_variant``.
    """

    aliases = {
        record.identifier: safe_identifier_aliases(record.identifier)
        for record in records
    }
    assigned = {record.identifier: [] for record in records}
    for path in files:
        lengths = {
            identifier: _matching_alias_length(path.name, query_aliases)
            for identifier, query_aliases in aliases.items()
        }
        longest = max(lengths.values(), default=0)
        if longest == 0:
            continue
        winners = [identifier for identifier, length in lengths.items() if length == longest]
        # Sanitization collisions are diagnosed separately by validate_chunk;
        # leaving the file unassigned avoids silently selecting either query.
        if len(winners) == 1:
            assigned[winners[0]].append(path)
    return assigned


def _sanitization_collisions(
    records: Sequence[FastaRecord],
) -> dict[str, set[str]]:
    by_alias: dict[str, set[str]] = {}
    for record in records:
        for alias in safe_identifier_aliases(record.identifier):
            by_alias.setdefault(alias.casefold(), set()).add(record.identifier)
    collisions: dict[str, set[str]] = {record.identifier: set() for record in records}
    for identifiers in by_alias.values():
        if len(identifiers) > 1:
            for identifier in identifiers:
                collisions[identifier].update(identifiers - {identifier})
    return collisions


def _load_plddt(score_path: Path) -> tuple[list[object] | None, str | None]:
    try:
        payload = json.loads(score_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        return None, f"cannot read scores JSON {score_path.name}: {exc}"
    if not isinstance(payload, dict):
        return None, f"scores JSON {score_path.name} is not an object"
    plddt = payload.get("plddt")
    if not isinstance(plddt, list):
        return None, f"scores JSON {score_path.name} has no plddt list"
    return plddt, None


def _validate_plddt(
    score_path: Path, plddt: Sequence[object], expected_length: int
) -> list[str]:
    reasons: list[str] = []
    if len(plddt) != expected_length:
        reasons.append(
            f"{score_path.name}: plddt length {len(plddt)} != sequence length {expected_length}"
        )
    invalid_positions: list[int] = []
    for index, value in enumerate(plddt, start=1):
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(float(value))
            or not 0.0 <= float(value) <= 100.0
        ):
            invalid_positions.append(index)
    if invalid_positions:
        preview = ",".join(str(index) for index in invalid_positions[:10])
        if len(invalid_positions) > 10:
            preview += ",..."
        reasons.append(
            f"{score_path.name}: plddt contains invalid value(s) at position(s) {preview}; expected finite numbers in [0, 100]"
        )
    return reasons


def validate_chunk(
    fasta_path: Path | str,
    output_dir: Path | str,
    phase: str,
) -> list[ValidationResult]:
    """Validate all FASTA records for the requested ColabFold phase."""

    if phase not in {"msa", "prediction"}:
        raise ValueError("phase must be 'msa' or 'prediction'")
    records = read_fasta(fasta_path)
    output = Path(output_dir)
    if not output.is_dir():
        return [
            ValidationResult(
                record.identifier, [f"output directory does not exist: {output}"]
            )
            for record in records
        ]

    files = sorted(path for path in output.rglob("*") if path.is_file())
    assigned = _assign_files(records, files)
    collisions = _sanitization_collisions(records)
    results: list[ValidationResult] = []

    for record in records:
        reasons: list[str] = []
        if collisions[record.identifier]:
            reasons.append(
                "filename sanitization collides with FASTA ID(s): "
                + ", ".join(sorted(collisions[record.identifier]))
            )
        query_files = assigned[record.identifier]

        if phase == "msa":
            pickles = [
                path for path in query_files if path.suffix.casefold() in PICKLE_SUFFIXES
            ]
            if not pickles:
                reasons.append("missing matching MSA pickle (.pickle or .pkl)")
        else:
            aliases = safe_identifier_aliases(record.identifier)
            done_files = [
                path
                for path in query_files
                if any(path.name.casefold() == f"{alias}.done.txt".casefold() for alias in aliases)
            ]
            if not done_files:
                reasons.append("missing matching .done.txt completion marker")

            score_files = [
                path
                for path in query_files
                if path.suffix.casefold() == ".json"
                and SCORE_RANK_ONE_RE.search(path.name)
            ]
            if not score_files:
                reasons.append("missing rank_001 scores JSON")
            elif len(score_files) > 1:
                reasons.append(
                    "expected exactly one rank_001 scores JSON, found "
                    + str(len(score_files))
                    + ": "
                    + ", ".join(path.name for path in score_files)
                )
            else:
                plddt, load_error = _load_plddt(score_files[0])
                if load_error:
                    reasons.append(load_error)
                else:
                    assert plddt is not None
                    reasons.extend(
                        _validate_plddt(score_files[0], plddt, len(record.sequence))
                    )

            structures = [
                path
                for path in query_files
                if path.suffix.casefold() in STRUCTURE_SUFFIXES
                and STRUCTURE_RANK_ONE_RE.search(path.name)
            ]
            if not structures:
                reasons.append("missing rank_001 PDB/CIF structure")

        results.append(ValidationResult(record.identifier, reasons))
    return results


def print_results(results: Sequence[ValidationResult]) -> None:
    for result in results:
        if result.ok:
            print(f"[OK] {result.identifier}")
        else:
            print(f"[FAIL] {result.identifier}: " + "; ".join(result.reasons))
    passed = sum(result.ok for result in results)
    print(f"Summary: {passed}/{len(results)} IDs passed")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Validate one ColabFold FASTA chunk and its output artifacts"
    )
    parser.add_argument("--fasta", required=True, type=Path, help="input chunk FASTA")
    parser.add_argument(
        "--output-dir", required=True, type=Path, help="ColabFold output directory"
    )
    parser.add_argument("--phase", required=True, choices=("msa", "prediction"))
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        results = validate_chunk(args.fasta, args.output_dir, args.phase)
    except ValueError as exc:
        print(f"[FAIL] input: {exc}", file=sys.stderr)
        return 2
    print_results(results)
    return 0 if all(result.ok for result in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
