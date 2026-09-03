#!/usr/bin/env python3
"""Prepare matched phenotype and covariate files for TensorQTL."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Convert a genes-by-samples expression matrix, hg38 gene positions, "
            "and covariates into matched TensorQTL inputs."
        )
    )
    parser.add_argument("--expression", required=True, type=Path)
    parser.add_argument("--gene-positions", required=True, type=Path)
    parser.add_argument("--covariates", required=True, type=Path)
    parser.add_argument("--fam", type=Path, help="Optional PLINK .fam file used to order samples")
    parser.add_argument("--output-bed", required=True, type=Path)
    parser.add_argument("--output-covariates", required=True, type=Path)
    parser.add_argument("--summary-json", type=Path)
    parser.add_argument("--expression-gene-column", default="gene_id")
    parser.add_argument("--position-gene-column", default="gene_id")
    parser.add_argument("--chrom-column", default="chrom")
    parser.add_argument("--tss-column", default="tss")
    parser.add_argument("--covariate-id-column", default="covariate_id")
    parser.add_argument("--max-missing-fraction", type=float, default=0.10)
    parser.add_argument("--inverse-normal", action="store_true")
    parser.add_argument("--include-sex-chromosomes", action="store_true")
    parser.add_argument("--require-complete-sample-match", action="store_true")
    return parser.parse_args()


def fail(message: str) -> None:
    raise ValueError(message)


def read_table(path: Path) -> pd.DataFrame:
    if not path.is_file():
        fail(f"Input file not found: {path}")
    return pd.read_csv(path, sep="\t", compression="infer")


def ensure_unique(values: pd.Index | pd.Series, label: str) -> None:
    duplicated = pd.Series(values, dtype="string").duplicated(keep=False)
    if duplicated.any():
        examples = pd.Series(values, dtype="string")[duplicated].drop_duplicates().head(5).tolist()
        fail(f"Duplicate {label}: {examples}")


def read_fam_ids(path: Path) -> list[str]:
    if not path.is_file():
        fail(f"FAM file not found: {path}")
    fam = pd.read_csv(path, sep=r"\s+", header=None, dtype=str)
    if fam.shape[1] < 2:
        fail(f"FAM file has fewer than two columns: {path}")
    ids = fam.iloc[:, 1].astype(str).tolist()
    ensure_unique(pd.Index(ids), "IID values in FAM")
    return ids


def normalize_chrom(value: object) -> str | None:
    text = str(value).strip()
    if text.lower().startswith("chr"):
        text = text[3:]
    text = text.upper()
    if text in {str(i) for i in range(1, 23)} | {"X", "Y"}:
        return f"chr{text}"
    return None


def inverse_normal_rows(values: pd.DataFrame) -> pd.DataFrame:
    transformed = np.empty(values.shape, dtype=float)
    for row_index, row in enumerate(values.to_numpy(dtype=float)):
        ranks = rankdata(row, method="average")
        transformed[row_index, :] = norm.ppf((ranks - 0.5) / len(row))
    return pd.DataFrame(transformed, index=values.index, columns=values.columns)


def write_bed(frame: pd.DataFrame, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.name.endswith(".bed.gz"):
        bgzip = shutil.which("bgzip")
        tabix = shutil.which("tabix")
        if not bgzip or not tabix:
            fail("bgzip and tabix are required when --output-bed ends in .bed.gz")
        plain = Path(str(output)[:-3])
        frame.to_csv(plain, sep="\t", index=False, float_format="%.8g")
        subprocess.run([bgzip, "-f", str(plain)], check=True)
        subprocess.run([tabix, "-f", "-p", "bed", str(output)], check=True)
    else:
        frame.to_csv(output, sep="\t", index=False, float_format="%.8g")


def main() -> int:
    args = parse_args()
    if not 0 <= args.max_missing_fraction <= 1:
        fail("--max-missing-fraction must be between 0 and 1")

    expression = read_table(args.expression)
    positions = read_table(args.gene_positions)
    covariates = read_table(args.covariates)

    for column, frame, label in [
        (args.expression_gene_column, expression, "expression"),
        (args.position_gene_column, positions, "gene positions"),
        (args.chrom_column, positions, "gene positions"),
        (args.tss_column, positions, "gene positions"),
        (args.covariate_id_column, covariates, "covariates"),
    ]:
        if column not in frame.columns:
            fail(f"Missing column '{column}' in {label}")

    expression[args.expression_gene_column] = expression[args.expression_gene_column].astype(str)
    positions[args.position_gene_column] = positions[args.position_gene_column].astype(str)
    covariates[args.covariate_id_column] = covariates[args.covariate_id_column].astype(str)
    ensure_unique(expression[args.expression_gene_column], "expression phenotype IDs")
    ensure_unique(positions[args.position_gene_column], "gene-position phenotype IDs")
    ensure_unique(covariates[args.covariate_id_column], "covariate IDs")

    expression_samples = [str(c) for c in expression.columns if c != args.expression_gene_column]
    covariate_samples = [str(c) for c in covariates.columns if c != args.covariate_id_column]
    ensure_unique(pd.Index(expression_samples), "expression sample IDs")
    ensure_unique(pd.Index(covariate_samples), "covariate sample IDs")

    expression.columns = [args.expression_gene_column] + expression_samples
    covariates.columns = [args.covariate_id_column] + covariate_samples
    expression_sample_set = set(expression_samples)
    covariate_sample_set = set(covariate_samples)

    if args.fam:
        genotype_samples = read_fam_ids(args.fam)
        sample_order = [
            sample
            for sample in genotype_samples
            if sample in expression_sample_set and sample in covariate_sample_set
        ]
        genotype_sample_set = set(genotype_samples)
    else:
        genotype_samples = []
        genotype_sample_set = expression_sample_set
        sample_order = [sample for sample in expression_samples if sample in covariate_sample_set]

    if len(sample_order) < 3:
        fail(f"Fewer than three matched samples were found: {len(sample_order)}")

    if args.require_complete_sample_match:
        common = set(sample_order)
        if common != expression_sample_set or common != covariate_sample_set:
            fail("Expression and covariate sample sets are not identical")
        if args.fam and not common.issubset(genotype_sample_set):
            fail("Some expression/covariate samples are absent from the FAM file")

    expression_values = expression.set_index(args.expression_gene_column)[sample_order]
    expression_values = expression_values.apply(pd.to_numeric, errors="coerce")
    missing_fraction = expression_values.isna().mean(axis=1)
    keep_missing = missing_fraction <= args.max_missing_fraction
    expression_values = expression_values.loc[keep_missing].copy()
    missing_values = int(expression_values.isna().sum().sum())
    expression_values = expression_values.T.fillna(expression_values.mean(axis=1)).T

    variance = expression_values.var(axis=1, ddof=1)
    keep_variable = variance.notna() & (variance > 0)
    expression_values = expression_values.loc[keep_variable].copy()
    if args.inverse_normal:
        expression_values = inverse_normal_rows(expression_values)

    positions = positions[
        [args.position_gene_column, args.chrom_column, args.tss_column]
    ].rename(columns={args.position_gene_column: "phenotype_id"})
    positions["#chr"] = positions[args.chrom_column].map(normalize_chrom)
    positions["tss"] = pd.to_numeric(positions[args.tss_column], errors="coerce")
    positions = positions.dropna(subset=["#chr", "tss"]).copy()
    positions["tss"] = positions["tss"].astype(int)
    positions = positions.loc[positions["tss"] >= 1].copy()
    if not args.include_sex_chromosomes:
        positions = positions.loc[positions["#chr"].isin([f"chr{i}" for i in range(1, 23)])]

    values = expression_values.copy()
    values.index.name = "phenotype_id"
    merged = positions.merge(values.reset_index(), on="phenotype_id", how="inner", validate="one_to_one")
    if merged.empty:
        fail("No phenotype IDs overlap between expression and gene-position files")

    merged["start"] = merged["tss"] - 1
    merged["end"] = merged["tss"]
    merged["chrom_order"] = merged["#chr"].str.removeprefix("chr").map(
        {**{str(i): i for i in range(1, 23)}, "X": 23, "Y": 24}
    )
    merged = merged.sort_values(["chrom_order", "start", "end", "phenotype_id"])
    bed = merged[["#chr", "start", "end", "phenotype_id"] + sample_order]

    covariate_values = covariates.set_index(args.covariate_id_column)[sample_order]
    covariate_values = covariate_values.apply(pd.to_numeric, errors="coerce")
    if covariate_values.isna().any().any():
        bad = covariate_values.index[covariate_values.isna().any(axis=1)].tolist()[:5]
        fail(f"Covariates contain missing or non-numeric values; examples: {bad}")
    covariate_variance = covariate_values.var(axis=1, ddof=1)
    constant_covariates = covariate_values.index[covariate_variance <= 0].tolist()
    covariate_values = covariate_values.loc[covariate_variance > 0].copy()
    covariate_values.index.name = "covariate_id"

    write_bed(bed, args.output_bed)
    args.output_covariates.parent.mkdir(parents=True, exist_ok=True)
    covariate_values.to_csv(args.output_covariates, sep="\t", float_format="%.8g")

    summary = {
        "expression_phenotypes_input": int(expression.shape[0]),
        "phenotypes_removed_for_missingness": int((~keep_missing).sum()),
        "phenotypes_removed_for_zero_variance": int((~keep_variable).sum()),
        "phenotypes_output": int(bed.shape[0]),
        "expression_samples_input": len(expression_samples),
        "covariate_samples_input": len(covariate_samples),
        "genotype_samples_input": len(genotype_samples) if args.fam else None,
        "matched_samples_output": len(sample_order),
        "missing_expression_values_mean_imputed": missing_values,
        "covariates_output": int(covariate_values.shape[0]),
        "constant_covariates_removed": constant_covariates,
        "inverse_normal_transform": bool(args.inverse_normal),
    }
    summary_path = args.summary_json or Path(f"{args.output_bed}.summary.json")
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, subprocess.CalledProcessError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise SystemExit(2)
