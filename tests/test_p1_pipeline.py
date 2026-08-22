from __future__ import annotations

import csv
import importlib
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np


PROJECT = Path(__file__).resolve().parents[1]
SCRIPTS = PROJECT / "scripts"
sys.path.insert(0, str(SCRIPTS))

utils = importlib.import_module("p1_utils")
ingest = importlib.import_module("ingest_colabfold_ensemble")
features = importlib.import_module("extract_pocket_features")


class P1PipelineTests(unittest.TestCase):
    def test_frozen_panel_produces_196_geometry_features_on_both_anchors(self) -> None:
        definition = json.loads(
            (PROJECT / "data" / "processed" / "pocket_definition" / "pocket_definition.json").read_text(encoding="utf-8")
        )
        positions = utils.pocket_positions(definition, panel_only=True)
        self.assertEqual(len(positions), 19)
        expected_names = None
        for path in (
            PROJECT / "data" / "raw" / "structures" / "5YWY.cif",
            PROJECT / "data" / "raw" / "structures" / "7D7M.cif",
        ):
            ca, sidechain = features.load_coordinate_maps(path)
            values, matrix = features.geometry_features(ca, sidechain, positions)
            self.assertEqual(len(values), 196)
            self.assertEqual(matrix.shape, (19, 3))
            expected_names = expected_names or list(values)
            self.assertEqual(list(values), expected_names)

    def test_kabsch_rmsd_is_invariant_to_rigid_transform(self) -> None:
        original = np.asarray([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]])
        rotation = np.asarray([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
        transformed = original @ rotation + np.asarray([10.0, -4.0, 2.0])
        self.assertAlmostEqual(features.kabsch_rmsd(original, transformed), 0.0, places=10)

    def test_synthetic_colabfold_output_ingests_and_reaches_feature_matrix(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            chunks = root / "chunks"
            raw = root / "raw" / "chunk_0001"
            output = root / "processed"
            feature_output = root / "features"
            chunks.mkdir(parents=True)
            raw.mkdir(parents=True)
            wt_sequence = "A" * 488
            fasta = chunks / "chunk_0001.fasta"
            fasta.write_text(f">PTGER4_WT\n{wt_sequence}\n", encoding="utf-8")
            wt_fasta = root / "wt.fasta"
            wt_fasta.write_text(f">PTGER4_WT\n{wt_sequence}\n", encoding="utf-8")
            manifest = chunks / "manifest.tsv"
            manifest.write_text(
                "array_index\tinput_fasta\tsequence_count\tminimum_length\tmaximum_length\n"
                "1\tchunk_0001.fasta\t1\t488\t488\n",
                encoding="utf-8",
            )
            tag = "alphafold2_ptm_model_1_seed_000"
            (raw / f"PTGER4_WT_scores_rank_001_{tag}.json").write_text(
                json.dumps({"plddt": [90.0] * 488, "ptm": 0.8}), encoding="utf-8"
            )
            pdb_lines = []
            for residue in range(1, 489):
                x = residue * 3.8
                pdb_lines.append(
                    f"ATOM  {residue:5d}  CA  ALA A{residue:4d}    {x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00 90.00           C"
                )
            pdb_lines.append("END")
            (raw / f"PTGER4_WT_unrelaxed_rank_001_{tag}.pdb").write_text("\n".join(pdb_lines) + "\n", encoding="utf-8")
            pocket_definition = PROJECT / "data" / "processed" / "pocket_definition" / "pocket_definition.json"
            argv = [
                "ingest_colabfold_ensemble.py", str(root / "raw"), "--run-label", "synthetic",
                "--batch-manifest", str(manifest), "--pocket-definition", str(pocket_definition),
                "--wt-fasta", str(wt_fasta), "--output-dir", str(output), "--mode", "sampling",
            ]
            with mock.patch.object(sys, "argv", argv):
                self.assertEqual(ingest.main(), 0)
            with (output / "conformer_manifest.csv").open(encoding="utf-8-sig", newline="") as handle:
                retained = list(csv.DictReader(handle))
            self.assertEqual(len(retained), 1)
            self.assertEqual(retained[0]["source_type"], "colabfold_wt")
            self.assertEqual(retained[0]["qc_status"], "WARN")
            argv = [
                "extract_pocket_features.py", "--manifest", str(output / "conformer_manifest.csv"),
                "--pocket-definition", str(pocket_definition), "--output-dir", str(feature_output),
            ]
            with mock.patch.object(sys, "argv", argv):
                self.assertEqual(features.main(), 0)
            schema = json.loads((feature_output / "feature_schema.json").read_text(encoding="utf-8"))
            self.assertEqual(schema["feature_count"], 196)
            with (feature_output / "pocket_features.csv").open(encoding="utf-8-sig", newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 3)  # two experimental anchors + one predicted conformer


if __name__ == "__main__":
    unittest.main()
