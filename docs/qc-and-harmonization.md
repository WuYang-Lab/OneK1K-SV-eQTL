# QC and harmonization notes

## Coordinate conversion is not allele harmonization

The original OneK1K workflow used UCSC liftOver to update BIM coordinates. This is appropriate for mapping positions, but it does not by itself resolve strand flips or establish which PLINK allele is the GRCh38 reference allele.

The public pipeline therefore adds a separate target-VCF preparation step that uses PLINK 2 `--ref-from-fa` and a GRCh38 FASTA. It then runs `bcftools norm --check-ref e`, which stops if unresolved REF mismatches remain.

Do not replace this with `bcftools norm --check-ref s` as a strand-fixing shortcut. If mismatches remain, determine whether they arise from the wrong build, strand convention, or ambiguous markers. Resolve them using platform-aware harmonization or exclude the affected sites. Compare allele frequencies against an ancestry-matched reference after harmonization.

## Chromosome naming

All VCFs used together must have identical contig names. The ImputeSV GRCh38 convention uses `chr1`–`chr22`. The pipeline exports the target using this convention and checks that the requested contig exists in both target and reference files.

## Reference-panel compatibility

Before running all chromosomes:

1. Confirm the target and panel use GRCh38.
2. Confirm both are bgzip-compressed and tabix-indexed.
3. Confirm sample IDs are unique.
4. Confirm the target consists of biallelic SNPs with valid REF/ALT alleles.
5. Compare target-panel overlap on chr22.
6. Check Beagle warnings for allele mismatches and excluded markers.

The imputation module writes a per-chromosome count table containing target, reference, and imputed record counts.

## Post-imputation filters

The example defaults reproduce the original analysis:

- `DR2 >= 0.3`
- `MAF >= 0.01`
- `MAF <= 0.99`
- genotype missingness `<= 0.05`
- HWE P value `>= 1e-6`

DR2 0.3 is permissive. A stricter threshold (for example 0.5 or 0.8) and sensitivity analyses may be appropriate for downstream fine-mapping or rare-variant interpretation. Filtering should be prespecified and reported.

## Dosages and PLINK formats

The chromosome QC step retains dosage information in PLINK 2 PGEN format. A PLINK 1 bed/bim/fam copy is also generated for compatibility with tools that require hard calls. Prefer PGEN/dosage-aware analysis when supported and document any hard-call conversion.

## Sample alignment

Genotype, expression, and covariate sample IDs must match exactly. `prepare_tensorqtl_inputs.py` preserves FAM sample order when `--fam` is supplied, and `validate_tensorqtl_inputs.py` fails on duplicates, missing samples, or ordering differences.

## Phenotype preparation

- Use one biological replicate/pseudobulk phenotype per donor and cell type.
- Apply expression filtering before QTL mapping.
- Remove phenotypes with zero variance.
- Treat missing expression values deliberately; the provided script performs row-mean imputation only after reporting their number.
- Use a documented normalization strategy and appropriate technical/biological covariates.
- Avoid including the tested interaction variable twice when performing interaction-QTL mapping.

## TensorQTL significance

Permutation-based `cis` mapping provides phenotype-level empirical/Beta-approximated P values and gene-level q-values. The `cis_independent` step should use the q-value-bearing `cis` output. The compatibility helper adds q-values only when TensorQTL did not generate them; Storey q-values are preferred and Benjamini–Hochberg is the fallback.

Nominal pairwise results should not be declared significant using an arbitrary raw P-value threshold without an analysis-specific multiple-testing procedure.

## Recommended pilot

Run chr22 and one well-powered cell type end to end. Record:

- starting, lifted, and unmapped SNP counts;
- target/reference overlap;
- imputed and retained SV counts;
- sample concordance across all inputs;
- number of tested genes and variants;
- TensorQTL warnings and runtime;
- number of FDR-significant eGenes and independent signals.
