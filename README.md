# Raw Genotype Quality Control Pipeline
### ALl parameters needs to become from Solène script to compare here results from the output of this pipe
### relatedness > pruning usare plink > PCA (recheck LD parameters in plink) **
## Overview

This pipeline performs quality control (QC) of raw genotype data.

The workflow includes sample and variant quality control, genetic sex verification, population structure assessment, and relatedness analysis to prepare high-quality genotype data for downstream genetic analyses.

## Pipeline Workflow

The pipeline consists of nine main steps:

| Step | Analysis | Description |
|---|---|---|
| 1 | **Input inspection** | Check the integrity and consistency of PLINK input files (BED/BIM/FAM), identifying duplicate sample IDs, duplicate variant IDs, and genomic positions. |
| 2 | **Sample selection** | Select eligible samples based on the target dataset and apply predefined sample exclusions. |
| 3 | **Genotype standardization** | Standardize chromosome naming, sort variants by genomic position, update reported sex information, and generate standardized PLINK files. |
| 4 | **Sample missingness QC** | Calculate the proportion of missing genotypes per individual and exclude samples exceeding the missingness threshold. |
| 5 | **Autosomal variant QC** | Filter autosomal variants based on minor allele frequency (MAF), variant missingness, and Hardy–Weinberg equilibrium (HWE). |
| 6 | **Heterozygosity QC** | Calculate individual heterozygosity rates and identify samples with unusually high or low heterozygosity. |
| 7 | **Genetic sex QC** | Infer genetic sex using X-chromosome markers, compare it with reported sex, and identify discordant or ambiguous samples. |
| 8 | **Population structure assessment** | Use `pcadapt` to identify variants associated with population structure and exclude these variants from the marker set used for relatedness analysis. |
| 9 | **Relatedness** | Estimate pairwise genetic relatedness using KING, identify related individuals, and generate a list of unrelated samples. |

## Quality Control Thresholds

The following thresholds are currently configured in the pipeline:

| QC metric | Threshold |
|---|---|
| Sample missingness | ≤ 10% |
| Variant missingness | ≤ 10% |
| Minor allele frequency (MAF) | ≥ 1% |
| Hardy–Weinberg equilibrium | p ≥ 1 × 10⁻¹⁵ |
| Heterozygosity | Within ±3 SD |
| Female X-chromosome inbreeding coefficient (F) | ≤ 0.49 |
| Male X-chromosome inbreeding coefficient (F) | ≥ 0.51 |
| KING kinship coefficient | 0.0442 |

Thresholds can be adjusted through the pipeline configuration file.

## Pipeline Outputs

The workflow generates:

- **Standardized genotype datasets:** PLINK BED/BIM/FAM files.
- **Sample QC results:** Missingness statistics and exclusion lists.
- **Variant QC results:** Filtered autosomal genotype datasets.
- **Heterozygosity QC results:** Individual heterozygosity statistics, diagnostic plots, and outlier lists.
- **Genetic sex QC results:** Sex inference statistics, discordance reports, and exclusion lists.
- **Population structure results:** `pcadapt` outputs and ancestry-associated variant lists.
- **Relatedness results:** KING pairwise kinship estimates and lists of unrelated individuals.

## Downstream Analyses

The QC outputs support the preparation of genotype datasets for downstream analyses, including:

- Genome-wide association studies (GWAS)
- Protein quantitative trait locus (pQTL) analyses
- Metabolite quantitative trait locus (mQTL) analyses
- Other genetic association studies