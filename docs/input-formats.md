# Input formats

## Genotypes

The starting genotype data are expected in PLINK 1 binary format:

```text
cohort.bed
cohort.bim
cohort.fam
```

`BFILE_HG19` is the common prefix without an extension. Variant IDs in the BIM file must be non-missing and unique. The default workflow processes autosomes only.

## Expression matrix

Tab-delimited, genes × samples. The first column is the phenotype/gene ID and all remaining columns are sample IDs.

```text
gene_id  sample_01  sample_02  sample_03
GENE1    -0.84      0.12       0.72
GENE2     0.31     -1.04       0.55
```

The values should already reflect the study's intended pseudobulk construction and expression normalization. `prepare_tensorqtl_inputs.py` can optionally inverse-normal transform each phenotype with `--inverse-normal`, but this is disabled by default.

## Gene positions

Tab-delimited with one row per phenotype/gene. Coordinates must be hg38/GRCh38 and TSS is 1-based.

```text
gene_id  chrom  tss
GENE1    chr1   100001
GENE2    chr2   200001
```

For strand-aware gene models, supply the transcription start site appropriate for the annotated strand. Alternative column names can be set with command-line arguments.

## Covariates

TensorQTL/FastQTL orientation: covariates × samples. The first column is a covariate identifier.

```text
covariate_id  sample_01  sample_02  sample_03
PC1           0.03       -0.11       0.08
sex           0           1          0
```

All covariates passed to TensorQTL must be numeric. Encode categorical variables before running the pipeline and avoid a redundant full set of dummy variables with an intercept.

## Cell-type list

One cell-type identifier per line. Blank lines and lines beginning with `#` are ignored.

```text
B_IN
B_MEM
CD4_NC
```

The identifier is substituted into `PHENOTYPE_TEMPLATE` and `COVARIATE_TEMPLATE` in the environment configuration.

## TensorQTL phenotype BED

The preparation script writes the required genes × samples BED format:

```text
#chr  start   end     phenotype_id  sample_01  sample_02  sample_03
chr1  100000  100001  GENE1         -0.84      0.12       0.72
```

`start` is 0-based and `end` is 1-based. The file is coordinate-sorted, bgzip-compressed, and tabix-indexed when the output filename ends in `.bed.gz`.
