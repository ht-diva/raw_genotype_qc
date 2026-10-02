Snakemake workflow for quality control of raw/imputed PLINK genotype data before downstream genetic analyses.

The pipeline starts from a PLINK BED/BIM/FAM dataset, validates and standardizes the input, restricts 
the cohort to eligible samples, matches genotype IDs to phenotype metadata, applies sample-level and 
variant-level QC, performs heterozygosity and sex checks, and then evaluates population structure and 
relatedness using `pcadapt` and KING.

