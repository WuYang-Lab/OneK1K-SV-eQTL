# OneK1K SV–sc-eQTL pipeline

A modular pipeline for structural-variant (SV) and VNTR imputation followed by single-cell/pseudobulk eQTL mapping.

The repository supports:

- PLINK genotype coordinate conversion from hg19/GRCh37 to hg38/GRCh38;
- target VCF preparation and reference-allele checks;
- chromosome-wise SV or VNTR imputation with Beagle and the ImputeSV reference panels;
- post-imputation filtering by imputation quality, allele frequency, missingness, and Hardy–Weinberg equilibrium;
- chromosome merging to PLINK 2 and PLINK 1 formats;
- preparation and validation of TensorQTL phenotype/covariate inputs;
- cis-nominal, permutation-based cis-eQTL, conditionally independent cis-eQTL, and optional trans-eQTL mapping;
- optional joint SNP+SV conditional mapping.


## Workflow

1. Lift SNP-array/WGS genotypes from hg19 to hg38.
2. Export a normalized, indexed target VCF with `chr` chromosome names.
3. Impute SVs or VNTRs chromosome by chromosome using an indexed ImputeSV reference panel.
4. Filter imputed variants and merge chromosomes.
5. Prepare pseudobulk expression and covariates in TensorQTL format.
6. Run TensorQTL cis-nominal, cis-permutation, and independent-signal mapping.
7. Optionally merge SNPs and SVs and repeat mapping for joint conditional analyses.

## Repository layout

```text
config/       Example project configuration and cell-type list
docs/         Input formats, QC guidance, and implementation notes
envs/         Conda environment
examples/     Small synthetic input tables
hpc/          SLURM array templates
scripts/      Executable pipeline modules
tests/        Lightweight tests
```

## Requirements

- Linux/bash
- PLINK 1.9 and PLINK 2
- UCSC `liftOver` and the `hg19ToHg38.over.chain` file
- bcftools, bgzip, and tabix
- Java and Beagle 5.x
- Python 3.10+, TensorQTL, PyTorch, pandas, scipy, and pyarrow
- R with `qvalue` and `data.table` (only needed when a TensorQTL cis result lacks a `qval` column)
- bgzip-indexed SV/VNTR reference panels from [ImputeSV](https://yanglab.westlake.edu.cn/ImputeSV/resource)
- an indexed GRCh38 FASTA matching the reference panel

Create the software environment:

```bash
conda env create -f envs/environment.yml
conda activate onek1k-sv-eqtl
```

Download Beagle separately and set `BEAGLE_JAR` in the configuration. Reference-panel access and use remain subject to the provider's terms.

## Configuration

Copy and edit the example configuration:

```bash
cp config/config.example.env config/config.env
```

All paths in `config/config.env` should be absolute. No cohort data, sample identifiers, or server-specific paths are stored in this repository.

Important parameters include:

| Variable | Purpose | Example default |
| --- | --- | ---: |
| `DR2_MIN` | Minimum Beagle dosage r² | `0.3` |
| `MAF_MIN` | Minimum minor-allele frequency | `0.01` |
| `MAF_MAX` | Maximum alternate-allele frequency | `0.99` |
| `GENO_MAX` | Maximum missing-call rate | `0.05` |
| `HWE_MIN` | Minimum HWE exact-test P value | `1e-6` |
| `CIS_WINDOW` | TensorQTL cis window in bp | `1000000` |
| `TENSORQTL_PERMUTATIONS` | Permutations per phenotype | `10000` |
| `TENSORQTL_FDR` | Gene-level FDR | `0.05` |

## Quick start

### 1. Convert hg19 PLINK coordinates to hg38

```bash
bash scripts/01_liftover_plink_hg19_to_hg38.sh config/config.env
```

This step retains successfully lifted autosomal variants, updates chromosome and position fields, reports unmapped and duplicate coordinates, and writes a new PLINK bed/bim/fam set.

Coordinate liftOver does not by itself prove allele-strand correctness. The next step resets/checks REF alleles against the configured GRCh38 FASTA and stops on unresolved mismatches. See [QC notes](docs/qc-and-harmonization.md).

### 2. Prepare the target VCF

```bash
bash scripts/02_prepare_target_vcf.sh config/config.env
```

The output is a sorted, bgzip-compressed, tabix-indexed, biallelic SNP VCF with `chr1`–`chr22` naming and unique `CHROM:POS:REF:ALT` IDs.

### 3. Impute one chromosome first

```bash
bash scripts/03_impute_chromosome.sh config/config.env SV 22
bash scripts/04_filter_imputed_chromosome.sh config/config.env SV 22
```

Inspect the logs and counts before launching all chromosomes. To process VNTRs, replace `SV` with `VNTR` and configure `VNTR_REFERENCE_VCF`.

### 4. Run all autosomes

Run locally/sequentially:

```bash
for chr in $(awk '!/^($|#)/ {print $1}' config/chromosomes.txt); do
  bash scripts/03_impute_chromosome.sh config/config.env SV "$chr"
  bash scripts/04_filter_imputed_chromosome.sh config/config.env SV "$chr"
done

bash scripts/05_merge_imputed_chromosomes.sh config/config.env SV
```

For a cluster, use the SLURM array examples in `hpc/`. The individual scripts are scheduler-independent and can also be wrapped by SGE/PBS or a site-specific submission helper.

### 5. Prepare TensorQTL inputs

Input schemas are described in [docs/input-formats.md](docs/input-formats.md). For each cell type:

```bash
python scripts/prepare_tensorqtl_inputs.py \
  --expression /path/to/B_IN.expression.tsv.gz \
  --gene-positions /path/to/gene_positions.hg38.tsv.gz \
  --covariates /path/to/B_IN.covariates.tsv \
  --fam /path/to/results/genotypes/SV.merged.bedfmt.fam \
  --output-bed /path/to/phenotypes/B_IN.bed.gz \
  --output-covariates /path/to/phenotypes/B_IN.covariates.tsv
```

The script intersects and orders samples, removes non-autosomal or zero-variance phenotypes, mean-imputes sparse phenotype missingness, and creates a tabix-indexed BED. It does **not** normalize expression unless `--inverse-normal` is specified; normalization should follow the study design.

Validate before mapping:

```bash
python scripts/validate_tensorqtl_inputs.py \
  --genotype-prefix /path/to/results/genotypes/SV.merged.bedfmt \
  --phenotype-bed /path/to/phenotypes/B_IN.bed.gz \
  --covariates /path/to/phenotypes/B_IN.covariates.tsv
```

### 6. Run TensorQTL

The configured file templates use `{cell_type}` as a placeholder:

```bash
bash scripts/06_run_tensorqtl.sh config/config.env B_IN sv
```

By default the script runs:

- `cis_nominal` for all variant–gene summary statistics;
- `cis` with permutations for gene-level FDR;
- `cis_independent` for stepwise conditional signals.

Run only selected modes by adding a fourth argument:

```bash
bash scripts/06_run_tensorqtl.sh config/config.env B_IN sv cis,cis_independent
bash scripts/06_run_tensorqtl.sh config/config.env B_IN sv trans
```

Run all cell types sequentially (or use the SLURM array template for parallel execution):

```bash
bash scripts/run_all_cell_types.sh config/config.env sv
```

Summarize independent signals across cell types:

```bash
python scripts/summarize_independent_eqtls.py \
  --base-dir /path/to/results/tensorqtl/sv \
  --cell-types config/cell_types.example.txt \
  --genotype-label sv \
  --output /path/to/results/tensorqtl/sv.independent_summary.tsv
```

### 7. Optional joint SNP+SV mapping

Create a combined PLINK bed set:

```bash
bash scripts/07_merge_snp_sv_genotypes.sh config/config.env
```

Then run the same TensorQTL workflow with the `joint` genotype label:

```bash
bash scripts/06_run_tensorqtl.sh config/config.env B_IN joint
```

This joint analysis is useful for asking whether an SV remains an independent signal after conditioning within the combined SNP+SV genotype set. It is not equivalent to proving that the SV is causal; LD-aware fine-mapping and functional evidence are still required.

## Reproducibility and QC

Before full-scale analysis, verify:

- genome build and chromosome naming match between target and reference panel;
- sample IDs are unique and identical across genotype, phenotype, and covariate files;
- reference alleles match GRCh38 after liftOver;
- allele frequencies are plausible relative to the reference panel;
- the Beagle output contains `DS`, `AF`, and `DR2` fields;
- variant and sample counts are recorded at every stage;
- a chr22 pilot completes successfully;
- tensorQTL phenotype coordinates are sorted and tabix-indexed.

Detailed guidance is in [docs/qc-and-harmonization.md](docs/qc-and-harmonization.md).

## Testing

The included tests exercise phenotype/covariate preparation on synthetic data:

```bash
python -m unittest discover -s tests -v
bash -n scripts/*.sh hpc/*.sh
```

Large genomic tools and real imputation are not run in continuous tests. Production validation should include one chromosome and one cell type on the target compute environment.

## Methodological notes

- The ImputeSV GRCh38 service/reference convention uses `chr`-prefixed chromosome names.
- TensorQTL `cis` mode supplies phenotype-level permutation statistics and, when the R/qvalue dependency is available, a `qval` column. `scripts/add_qvalues.R` is only a compatibility fallback.
- `cis_independent` uses the significant phenotypes from the permutation result and performs stepwise conditional mapping.
- Trans-eQTL scans have a much larger multiple-testing burden. They are optional and disabled unless explicitly requested.
- Default filters reflect the original OneK1K implementation; sensitivity analyses with stricter DR2 thresholds are recommended.

## Data availability

This repository contains code and synthetic examples only. It does not distribute OneK1K genotypes, expression matrices, individual-level results, or the ImputeSV reference panels.

## Citation

If this workflow contributes to a publication, cite TensorQTL, Beagle, PLINK/PLINK2, UCSC liftOver, and the ImputeSV resource/publication, together with the final OneK1K SV–sc-eQTL study when available. A repository citation record is provided in `CITATION.cff` and should be updated with the publication DOI before release.

## License

BSD 3-Clause License. See [LICENSE](LICENSE).
