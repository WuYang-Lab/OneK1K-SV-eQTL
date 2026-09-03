#!/usr/bin/env python3
"""Validate sample alignment and coordinate fields for TensorQTL inputs."""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from io import StringIO
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--genotype-prefix", required=True, type=Path)
    parser.add_argument("--phenotype-bed", required=True, type=Path)
    parser.add_argument("--covariates", required=True, type=Path)
    parser.add_argument("--output-json", type=Path)
    return parser.parse_args()


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open("rt", encoding="utf-8")


def duplicate_examples(values: list[str]) -> list[str]:
    series = pd.Series(values, dtype="string")
    return series[series.duplicated(keep=False)].drop_duplicates().head(5).tolist()


def read_genotype_metadata(prefix: Path) -> tuple[list[str], set[str], str]:
    fam = Path(f"{prefix}.fam")
    bim = Path(f"{prefix}.bim")
    psam = Path(f"{prefix}.psam")
    pvar = Path(f"{prefix}.pvar")
    if fam.is_file() and bim.is_file():
        fam_df = pd.read_csv(fam, sep=r"\s+", header=None, dtype=str)
        samples = fam_df.iloc[:, 1].astype(str).tolist()
        bim_df = pd.read_csv(bim, sep=r"\s+", header=None, usecols=[0], dtype=str)
        chromosomes = set(bim_df.iloc[:, 0].str.removeprefix("chr"))
        return samples, chromosomes, "bed"
    if psam.is_file() and pvar.is_file():
        psam_df = pd.read_csv(psam, sep=r"\s+", dtype=str)
        iid_column = "#IID" if "#IID" in psam_df.columns else "IID"
        if iid_column not in psam_df.columns:
            raise ValueError(f"No IID column found in {psam}")
        samples = psam_df[iid_column].astype(str).tolist()
        with pvar.open("rt", encoding="utf-8") as handle:
            pvar_text = "".join(line for line in handle if not line.startswith("##"))
        pvar_df = pd.read_csv(StringIO(pvar_text), sep="\t", dtype=str)
        chrom_column = "#CHROM" if "#CHROM" in pvar_df.columns else "CHROM"
        chromosomes = set(pvar_df[chrom_column].str.removeprefix("chr"))
        return samples, chromosomes, "pgen"
    raise ValueError(f"No complete bed/bim/fam or pgen/pvar/psam set found for {prefix}")


def main() -> int:
    args = parse_args()
    genotype_samples, genotype_chromosomes, genotype_format = read_genotype_metadata(
        args.genotype_prefix
    )
    genotype_duplicates = duplicate_examples(genotype_samples)
    if genotype_duplicates:
        raise ValueError(f"Duplicate genotype sample IDs: {genotype_duplicates}")

    if not args.phenotype_bed.is_file():
        raise ValueError(f"Phenotype BED not found: {args.phenotype_bed}")
    with open_text(args.phenotype_bed) as handle:
        header = handle.readline().rstrip("\n").split("\t")
    if header[:4] != ["#chr", "start", "end", "phenotype_id"]:
        raise ValueError("Phenotype BED must start with #chr, start, end, phenotype_id")
    phenotype_samples = header[4:]
    phenotype_duplicates = duplicate_examples(phenotype_samples)
    if phenotype_duplicates:
        raise ValueError(f"Duplicate phenotype sample IDs: {phenotype_duplicates}")

    if not args.covariates.is_file():
        raise ValueError(f"Covariate file not found: {args.covariates}")
    covariates = pd.read_csv(args.covariates, sep="\t", index_col=0)
    covariate_samples = covariates.columns.astype(str).tolist()
    if phenotype_samples != covariate_samples:
        raise ValueError("Phenotype and covariate samples or their order differ")

    missing_from_genotype = sorted(set(phenotype_samples) - set(genotype_samples))
    if missing_from_genotype:
        raise ValueError(
            f"Phenotype samples absent from genotype data: {missing_from_genotype[:5]}"
        )

    phenotype = pd.read_csv(args.phenotype_bed, sep="\t", compression="infer")
    if phenotype["phenotype_id"].duplicated().any():
        raise ValueError("Duplicate phenotype IDs in phenotype BED")
    starts = pd.to_numeric(phenotype["start"], errors="coerce")
    ends = pd.to_numeric(phenotype["end"], errors="coerce")
    if starts.isna().any() or ends.isna().any() or (starts < 0).any() or (ends <= starts).any():
        raise ValueError("Invalid phenotype BED coordinates")

    bed_chromosomes = phenotype["#chr"].astype(str).str.removeprefix("chr")
    missing_chromosomes = sorted(set(bed_chromosomes) - genotype_chromosomes)
    if missing_chromosomes:
        raise ValueError(f"Phenotype chromosomes absent from genotype set: {missing_chromosomes}")
    chrom_order = bed_chromosomes.map({str(i): i for i in range(1, 23)} | {"X": 23, "Y": 24})
    if chrom_order.isna().any():
        raise ValueError("Unsupported chromosome labels in phenotype BED")
    order_frame = pd.DataFrame({"chrom": chrom_order, "start": starts, "end": ends})
    if not order_frame.equals(order_frame.sort_values(["chrom", "start", "end"])):
        raise ValueError("Phenotype BED is not coordinate-sorted")

    numeric_covariates = covariates.apply(pd.to_numeric, errors="coerce")
    if numeric_covariates.isna().any().any():
        raise ValueError("Covariates contain missing or non-numeric values")

    summary = {
        "status": "ok",
        "genotype_format": genotype_format,
        "genotype_samples": len(genotype_samples),
        "phenotype_samples": len(phenotype_samples),
        "phenotypes": int(phenotype.shape[0]),
        "covariates": int(covariates.shape[0]),
        "extra_genotype_samples": len(set(genotype_samples) - set(phenotype_samples)),
        "chromosomes": sorted(set(phenotype["#chr"].astype(str))),
    }
    if args.output_json:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ValueError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise SystemExit(2)
