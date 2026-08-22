from __future__ import annotations

import contextlib
import importlib.util
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "validate_colabfold_chunk.py"
SPEC = importlib.util.spec_from_file_location("validate_colabfold_chunk", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
validator = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = validator
SPEC.loader.exec_module(validator)


class ValidateColabFoldChunkTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        self.fasta = self.root / "chunk.fasta"
        self.output = self.root / "output"
        self.output.mkdir()

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def write_fasta(self, text: str) -> None:
        self.fasta.write_text(text, encoding="utf-8")

    def write_prediction(
        self,
        identifier: str,
        plddt: list[object],
        structure_suffix: str = ".pdb",
    ) -> None:
        (self.output / f"{identifier}.done.txt").write_text("done\n", encoding="utf-8")
        (self.output / f"{identifier}_scores_rank_001_model_1_seed_000.json").write_text(
            json.dumps({"plddt": plddt}), encoding="utf-8"
        )
        (self.output / f"{identifier}_unrelaxed_rank_001_model_1_seed_000{structure_suffix}").write_text(
            "MODEL\nEND\n", encoding="utf-8"
        )

    def test_msa_accepts_common_sanitized_identifier_and_reports_missing_id(self) -> None:
        self.write_fasta(">sp|P35408|EP4 description\nACDE\n>missing\nAAAA\n")
        (self.output / "sp_P35408_EP4.pickle").write_bytes(b"placeholder")

        results = validator.validate_chunk(self.fasta, self.output, "msa")

        self.assertTrue(results[0].ok)
        self.assertFalse(results[1].ok)
        self.assertIn("missing matching MSA pickle", results[1].reasons[0])

    def test_prediction_accepts_valid_rank_one_artifacts(self) -> None:
        self.write_fasta(">PTGER4_WT\nACDE\n")
        self.write_prediction("PTGER4_WT", [91.0, 82, 73.5, 64], ".cif")

        results = validator.validate_chunk(self.fasta, self.output, "prediction")

        self.assertEqual(len(results), 1)
        self.assertTrue(results[0].ok, results[0].reasons)

    def test_prediction_reports_all_missing_artifacts_per_identifier(self) -> None:
        self.write_fasta(">query_a\nAAA\n>query_b\nCCCC\n")

        results = validator.validate_chunk(self.fasta, self.output, "prediction")

        self.assertEqual([result.identifier for result in results], ["query_a", "query_b"])
        for result in results:
            combined = "; ".join(result.reasons)
            self.assertIn(".done.txt", combined)
            self.assertIn("scores JSON", combined)
            self.assertIn("PDB/CIF", combined)

    def test_prediction_rejects_duplicate_scores_and_bad_plddt(self) -> None:
        self.write_fasta(">query\nAAAA\n")
        self.write_prediction("query", [90, 80, 70, 60])
        (self.output / "query_scores_rank_001_model_2_seed_000.json").write_text(
            json.dumps({"plddt": [90, 80, 70, 101]}), encoding="utf-8"
        )

        duplicate = validator.validate_chunk(self.fasta, self.output, "prediction")[0]
        self.assertFalse(duplicate.ok)
        self.assertIn("exactly one", "; ".join(duplicate.reasons))

        (self.output / "query_scores_rank_001_model_2_seed_000.json").unlink()
        score = self.output / "query_scores_rank_001_model_1_seed_000.json"
        score.write_text(json.dumps({"plddt": [90, 80, 101]}), encoding="utf-8")
        invalid = validator.validate_chunk(self.fasta, self.output, "prediction")[0]
        combined = "; ".join(invalid.reasons)
        self.assertIn("plddt length 3 != sequence length 4", combined)
        self.assertIn("invalid value", combined)

    def test_longest_prefix_prevents_cross_assignment(self) -> None:
        self.write_fasta(">sample\nAAA\n>sample_variant\nAAA\n")
        (self.output / "sample.pickle").write_bytes(b"one")
        (self.output / "sample_variant.pickle").write_bytes(b"two")

        results = validator.validate_chunk(self.fasta, self.output, "msa")

        self.assertTrue(all(result.ok for result in results))

    def test_main_returns_nonzero_and_prints_each_id_reason(self) -> None:
        self.write_fasta(">first\nAA\n>second\nAA\n")
        stream = io.StringIO()
        with contextlib.redirect_stdout(stream):
            exit_code = validator.main(
                [
                    "--fasta",
                    str(self.fasta),
                    "--output-dir",
                    str(self.output),
                    "--phase",
                    "msa",
                ]
            )

        output = stream.getvalue()
        self.assertEqual(exit_code, 1)
        self.assertIn("[FAIL] first:", output)
        self.assertIn("[FAIL] second:", output)
        self.assertIn("Summary: 0/2 IDs passed", output)


if __name__ == "__main__":
    unittest.main()
