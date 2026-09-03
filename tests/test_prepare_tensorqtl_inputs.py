from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "prepare_tensorqtl_inputs.py"
EXAMPLES = REPO_ROOT / "examples"


class PrepareTensorQTLInputsTests(unittest.TestCase):
    def test_prepare_filters_and_orders(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            tmp_path = Path(directory)
            output_bed = tmp_path / "phenotype.bed"
            output_covariates = tmp_path / "covariates.tsv"
            summary_json = tmp_path / "summary.json"

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--expression",
                    str(EXAMPLES / "expression.tsv"),
                    "--gene-positions",
                    str(EXAMPLES / "gene_positions.hg38.tsv"),
                    "--covariates",
                    str(EXAMPLES / "covariates.tsv"),
                    "--output-bed",
                    str(output_bed),
                    "--output-covariates",
                    str(output_covariates),
                    "--summary-json",
                    str(summary_json),
                    "--max-missing-fraction",
                    "0.5",
                ],
                check=True,
                capture_output=True,
                text=True,
            )

            self.assertEqual(result.returncode, 0)
            bed = pd.read_csv(output_bed, sep="\t")
            self.assertEqual(
                bed.columns.tolist(),
                [
                    "#chr",
                    "start",
                    "end",
                    "phenotype_id",
                    "sample_01",
                    "sample_02",
                    "sample_03",
                ],
            )
            self.assertEqual(bed["phenotype_id"].tolist(), ["GENE1"])
            self.assertEqual(bed.loc[0, "start"], 100000)
            self.assertEqual(bed.loc[0, "end"], 100001)

            covariates = pd.read_csv(output_covariates, sep="\t")
            self.assertEqual(covariates["covariate_id"].tolist(), ["PC1", "sex"])

            summary = json.loads(summary_json.read_text(encoding="utf-8"))
            self.assertEqual(summary["phenotypes_output"], 1)
            self.assertEqual(summary["matched_samples_output"], 3)
            self.assertEqual(summary["phenotypes_removed_for_zero_variance"], 1)

    def test_inverse_normal_option(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            tmp_path = Path(directory)
            output_bed = tmp_path / "phenotype.bed"
            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--expression",
                    str(EXAMPLES / "expression.tsv"),
                    "--gene-positions",
                    str(EXAMPLES / "gene_positions.hg38.tsv"),
                    "--covariates",
                    str(EXAMPLES / "covariates.tsv"),
                    "--output-bed",
                    str(output_bed),
                    "--output-covariates",
                    str(tmp_path / "covariates.tsv"),
                    "--summary-json",
                    str(tmp_path / "summary.json"),
                    "--inverse-normal",
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0)
            bed = pd.read_csv(output_bed, sep="\t")
            mean = float(bed[["sample_01", "sample_02", "sample_03"]].mean(axis=1).iloc[0])
            self.assertLess(abs(mean), 1e-8)


if __name__ == "__main__":
    unittest.main()
