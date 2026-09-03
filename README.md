# Raw genotype sample-QC pipeline

Modular Snakemake implementation of the supplied DAG, from raw PLINK BED/BIM/FAM files to a final sample set, PCA covariates, and a REGENIE Step-1 BED dataset.

## Workflow

1. Standardize chromosome coding and variant IDs, match genotype/metadata IDs, update reported sex, and report duplicate samples, variant IDs, and chromosome-position pairs.
2. Apply an optional external sample exclusion list.
3. Calculate autosomal sample missingness, plot its distribution, and retain samples with missingness at or below `0.10` by default.
4. Run two branches in parallel:
   - common autosomal QC: autosomes, MAF `0.01`, variant missingness `0.10`, HWE `1e-15`, and heterozygosity outliers outside mean ± 3 SD;
   - chromosome-X sex QC after PAR splitting, with configurable F-statistic thresholds.
5. Stop at review gate 1. The analyst creates an approved two-column FID/IID keep file.
6. Create the structure marker set, LD-prune it (`50 kb`, `5`, `r²=0.5`), and retain the unrelated PCA reference selected by KING cutoff `0.0442`.
7. Compute PCA loadings in unrelated samples and project all reviewed samples.
8. Stop at review gate 2 for population-structure review and creation of the final keep file.
9. Recompute final PCA and export the final full-marker BED/BIM/FAM dataset for REGENIE Step 1.

All numerical choices are parameters in `config/config.yaml`; the defaults reproduce the DAG.

## Required inputs

- PLINK 1 BED/BIM/FAM files.
- `resources/sample_metadata.tsv`, tab-separated with `FID`, `IID`, and `SEX` (`1` male, `2` female, `0` unknown).
- `resources/exclusions/external.remove`, a two-column FID/IID file. Leave it comment-only when there are no exclusions.

Edit the paths, genome-build-dependent PAR setting, resources, and thresholds in `config/config.yaml`. The current sex-QC command uses `--split-par b37`; change it to the correct build before running.

## Run

```bash
mamba env create -f environment.yaml
conda activate raw-genotype-qc
snakemake --dry-run
snakemake --cores 8
```

On SLURM:

```bash
snakemake --profile profiles/slurm
```

## Manual review gates

First generate the gate-1 evidence:

```bash
snakemake \
  results/cohort/03_autosomal_qc/heterozygosity_outliers.remove \
  results/cohort/04_sex_qc/candidate_exclusions.remove \
  --cores 8
```

Review both reports and write all approved samples—not only exclusions—to `resources/review/gate1.keep`. Then generate the projected PCA scores:

```bash
snakemake results/cohort/07_pca/projected.sscore --cores 8
```

Review the PCA/population structure and write the final approved samples to `resources/review/gate2.keep`. Finally run the complete workflow.

The gate rules reject empty keep files and IDs absent from the preceding sample set, preventing an accidental silent pass.

## Important adaptation points

- Confirm whether the cohort uses GRCh37 or GRCh38 and adjust `--split-par`.
- Confirm PLINK sex coding and metadata ID columns.
- Add cohort-specific long-range-LD or inversion exclusions if required for PCA.
- The final REGENIE dataset keeps all variants remaining after external exclusions and only filters samples. Add variant-level REGENIE requirements if your analysis protocol calls for them.
- Test on a small dataset before the production SLURM run.
