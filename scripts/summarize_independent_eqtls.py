#!/usr/bin/env python3
"""Summarize conditionally independent TensorQTL signals across cell types."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-dir", required=True, type=Path)
    parser.add_argument("--cell-types", required=True, type=Path)
    parser.add_argument("--genotype-label", default="sv")
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def read_cell_types(path: Path) -> list[str]:
    return [
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]


def main() -> int:
    args = parse_args()
    rows: list[dict[str, object]] = []
    for cell_type in read_cell_types(args.cell_types):
        path = (
            args.base_dir
            / cell_type
            / f"{cell_type}.{args.genotype_label}.cis_independent_qtl.txt.gz"
        )
        if not path.is_file():
            rows.append({"cell_type": cell_type, "status": "missing"})
            continue
        frame = pd.read_csv(path, sep="\t")
        required = {"phenotype_id", "variant_id"}
        missing = required - set(frame.columns)
        if missing:
            rows.append(
                {
                    "cell_type": cell_type,
                    "status": f"missing_columns:{','.join(sorted(missing))}",
                }
            )
            continue
        per_gene = frame.groupby("phenotype_id", observed=True).size()
        rows.append(
            {
                "cell_type": cell_type,
                "status": "ok",
                "independent_signals": int(frame.shape[0]),
                "independent_eGenes": int(frame["phenotype_id"].nunique()),
                "unique_variants": int(frame["variant_id"].nunique()),
                "genes_with_multiple_signals": int((per_gene > 1).sum()),
                "max_rank": int(frame["rank"].max()) if "rank" in frame and not frame.empty else 0,
                "file": str(path),
            }
        )

    output = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, sep="\t", index=False)
    print(output.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
