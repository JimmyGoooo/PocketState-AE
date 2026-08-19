#!/usr/bin/env python3
"""Validate a FASTA panel and split it into deterministic ColabFold batches."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    title: str | None = None
    sequence: list[str] = []
    for line_number, raw_line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if title is not None:
                records.append((title, "".join(sequence)))
            title, sequence = line[1:].strip(), []
            if not title or any(char.isspace() for char in title):
                raise ValueError(f"Line {line_number}: FASTA identifier must be non-empty and whitespace-free")
        elif title is None:
            raise ValueError(f"Line {line_number}: sequence found before FASTA identifier")
        else:
            sequence.append(line.upper())
    if title is not None:
        records.append((title, "".join(sequence)))
    if not records:
        raise ValueError("No FASTA records found")
    identifiers = [identifier for identifier, _ in records]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("FASTA identifiers must be unique")
    for identifier, seq in records:
        if not seq or set(seq) - VALID_AA:
            raise ValueError(f"{identifier}: sequence is empty or contains non-standard amino-acid characters")
    return records


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for identifier, sequence in records:
            handle.write(f">{identifier}\n")
            for offset in range(0, len(sequence), 80):
                handle.write(sequence[offset:offset + 80] + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_fasta", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--chunk-size", type=int, default=25)
    args = parser.parse_args()
    if args.chunk_size < 1:
        raise ValueError("--chunk-size must be positive")
    records = read_fasta(args.input_fasta)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows = []
    for number, start in enumerate(range(0, len(records), args.chunk_size), start=1):
        chunk = records[start:start + args.chunk_size]
        filename = f"chunk_{number:04d}.fasta"
        write_fasta(args.output_dir / filename, chunk)
        manifest_rows.append({
            "array_index": number,
            "input_fasta": filename,
            "sequence_count": len(chunk),
            "minimum_length": min(len(sequence) for _, sequence in chunk),
            "maximum_length": max(len(sequence) for _, sequence in chunk),
        })
    with (args.output_dir / "manifest.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)
    print(f"Validated {len(records)} sequences and wrote {len(manifest_rows)} chunks to {args.output_dir}")


if __name__ == "__main__":
    main()
