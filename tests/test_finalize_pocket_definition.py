from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_ROOT / "scripts" / "finalize_pocket_definition.py"
SPEC = importlib.util.spec_from_file_location("finalize_pocket_definition", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class FinalizePocketDefinitionTests(unittest.TestCase):
    def test_real_input_freezes_38_positions_and_19_feature_sites(self) -> None:
        expected_panel = [
            23, 24, 27, 69, 72, 73, 76, 80, 99, 166,
            168, 169, 312, 315, 316, 317, 318, 319, 320,
        ]
        with tempfile.TemporaryDirectory() as temporary_directory:
            output_dir = Path(temporary_directory) / "nested" / "definition"
            csv_path, json_path = MODULE.finalize_pocket_definition(
                MODULE.DEFAULT_INPUT, output_dir
            )

            self.assertEqual(csv_path, output_dir / "pocket_residues.csv")
            self.assertEqual(json_path, output_dir / "pocket_definition.json")
            document = json.loads(json_path.read_text(encoding="utf-8"))
            self.assertEqual(document["schema_version"], "1.0")
            self.assertEqual(document["target"]["uniprot_id"], "P35408")
            self.assertEqual(document["cutoff_angstrom"], 6.0)
            self.assertEqual(document["position_count"], 38)
            self.assertEqual(document["feature_panel"]["position_count"], 19)
            self.assertEqual(document["feature_panel"]["positions"], expected_panel)

            positions = document["positions"]
            self.assertEqual(
                [item["uniprot_position"] for item in positions],
                sorted(item["uniprot_position"] for item in positions),
            )
            self.assertEqual(len({item["uniprot_position"] for item in positions}), 38)
            by_position = {item["uniprot_position"]: item for item in positions}
            self.assertFalse(by_position[20]["anchors"]["5YWY_7UR"]["present"])
            self.assertIsNone(
                by_position[20]["anchors"]["5YWY_7UR"][
                    "minimum_heavy_atom_distance_angstrom"
                ]
            )
            self.assertTrue(by_position[20]["anchors"]["7D7M_P2E"]["present"])
            self.assertTrue(by_position[23]["anchors"]["5YWY_7UR"]["present"])
            self.assertTrue(by_position[23]["anchors"]["7D7M_P2E"]["present"])
            self.assertEqual(
                by_position[320]["feature_panel_reason"],
                "nearest_specific_to_5YWY_7UR",
            )
            self.assertEqual(
                by_position[69]["feature_panel_reason"],
                "nearest_specific_to_7D7M_P2E",
            )

            expected_hash = hashlib.sha256(MODULE.DEFAULT_INPUT.read_bytes()).hexdigest()
            self.assertEqual(document["sources"]["input_csv"]["sha256"], expected_hash)
            self.assertEqual(
                [anchor["anchor_id"] for anchor in document["sources"]["anchors"]],
                ["5YWY_7UR", "7D7M_P2E"],
            )

            with csv_path.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 38)
            self.assertEqual([int(row["order"]) for row in rows], list(range(1, 39)))
            self.assertEqual(sum(row["in_feature_panel"] == "YES" for row in rows), 19)

    def test_cli_overrides_input_and_output_and_is_deterministic(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            temporary = Path(temporary_directory)
            copied_input = temporary / "contacts.csv"
            copied_input.write_bytes(MODULE.DEFAULT_INPUT.read_bytes())
            output_dir = temporary / "output"

            self.assertEqual(
                MODULE.main(
                    ["--input", str(copied_input), "--output-dir", str(output_dir)]
                ),
                0,
            )
            first_csv = (output_dir / "pocket_residues.csv").read_bytes()
            first_json = (output_dir / "pocket_definition.json").read_bytes()
            MODULE.main(["--input", str(copied_input), "--output-dir", str(output_dir)])
            self.assertEqual((output_dir / "pocket_residues.csv").read_bytes(), first_csv)
            self.assertEqual((output_dir / "pocket_definition.json").read_bytes(), first_json)

    def test_conflicting_residue_names_are_rejected(self) -> None:
        rows = self._real_rows()
        duplicate = dict(rows[0])
        duplicate["pdb_id"] = "7D7M"
        duplicate["ligand_id"] = "P2E"
        duplicate["residue_name_3letter"] = "ZZZ"
        duplicate["minimum_heavy_atom_distance_angstrom"] = "1.0"
        with tempfile.TemporaryDirectory() as temporary_directory:
            input_csv = Path(temporary_directory) / "conflict.csv"
            self._write_rows(input_csv, [*rows, duplicate])
            with self.assertRaisesRegex(ValueError, "conflicting residue names"):
                MODULE.load_pocket_positions(input_csv)

    def test_conflicting_duplicate_anchor_distance_is_rejected(self) -> None:
        rows = self._real_rows()
        duplicate = dict(rows[0])
        duplicate["minimum_heavy_atom_distance_angstrom"] = "0.5"
        with tempfile.TemporaryDirectory() as temporary_directory:
            input_csv = Path(temporary_directory) / "conflict.csv"
            self._write_rows(input_csv, [*rows, duplicate])
            with self.assertRaisesRegex(ValueError, "conflicting .* distances"):
                MODULE.load_pocket_positions(input_csv)

    def test_unknown_anchor_is_rejected(self) -> None:
        rows = self._real_rows()
        rows[0]["pdb_id"] = "XXXX"
        with tempfile.TemporaryDirectory() as temporary_directory:
            input_csv = Path(temporary_directory) / "unknown.csv"
            self._write_rows(input_csv, rows)
            with self.assertRaisesRegex(ValueError, "unexpected anchor pair"):
                MODULE.load_pocket_positions(input_csv)

    @staticmethod
    def _real_rows() -> list[dict[str, str]]:
        with MODULE.DEFAULT_INPUT.open(encoding="utf-8", newline="") as handle:
            return list(csv.DictReader(handle))

    @staticmethod
    def _write_rows(path: Path, rows: list[dict[str, str]]) -> None:
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)


if __name__ == "__main__":
    unittest.main()
