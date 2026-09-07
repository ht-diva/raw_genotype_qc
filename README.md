### Input data standardization

The `standardize_input_data.smk` module validates and standardizes the input PLINK BED dataset before downstream quality-control analyses.

The module performs the following steps:

1. **Input inspection:** verifies that the BED, BIM, and FAM files exist and are not empty; checks their basic structure; identifies duplicated sample IDs, duplicated variant IDs, duplicated genomic positions, and missing variant IDs; and generates a chromosome-renaming map.

2. **Reported-sex preparation:** matches genotype sample IDs with the phenotype metadata, converts reported sex to PLINK-compatible codes, and reports genotype or metadata records that could not be matched.

3. **Genotype standardization:** updates the reported sex, standardizes chromosome labels, and sorts variants by chromosome, genomic position, and variant ID. Since PLINK 2 cannot sort variants while directly writing a fixed-width BED file, the dataset is first written as a temporary PGEN fileset and then converted back to BED format.

The final standardized dataset is written to:

```text
results/standardization/genotype.bed
results/standardization/genotype.bim
results/standardization/genotype.fam
```

These standardized files are used as input for the subsequent quality-control steps.
