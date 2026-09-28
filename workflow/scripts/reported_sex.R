    )

    ifelse(
        is.na(standardized_values),
        "0",
        ifelse(
            standardized_values %in% male_values,
            "1",
            ifelse(
                standardized_values %in% female_values,
                "2",
                "0"
            )
        )
    )
}

phenotype$PLINK_SEX <- sex_to_plink(
    phenotype[[sex_column]]
)


# ----------------------------------------------------------
# 15. Match genotype samples to phenotype file
# ----------------------------------------------------------

# For every genotype sample, return the position of its MATCH_ID
# in the phenotype. An unmatched sample receives NA.
match_index <- match(
    fam$MATCH_ID,
    phenotype$MATCH_ID
)
# A failed genotype ID extraction must never match a missing phenotype ID.
match_index[is.na(fam$MATCH_ID)] <- NA_integer_

# Use the phenotype sex when a match exists. Assign PLINK sex code
# 0 when the genotype sample is not present in the phenotype.
plink_sex <- ifelse(
    is.na(match_index),
    "0",
    phenotype$PLINK_SEX[match_index]
)


# ----------------------------------------------------------
# 16. Write the PLINK --update-sex file
# ----------------------------------------------------------

# PLINK expects three columns without a header:
# FID IID SEX
write.table(
    data.frame(
        FID = fam$FID,
        IID = fam$IID,
        SEX = plink_sex
    ),
    file = opt$`update-sex`,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
)


# ----------------------------------------------------------
# 17. Write genotype samples not found in the phenotype
# ----------------------------------------------------------

unmatched_genotype_rows <- is.na(match_index)

write.table(
    data.frame(
        FID = fam$FID[unmatched_genotype_rows],
        IID = fam$IID[unmatched_genotype_rows],
        MATCH_ID = fam$MATCH_ID[unmatched_genotype_rows]
    ),
    file = opt$`unmatched-genotypes`,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)


# ----------------------------------------------------------
# 18. Write phenotype IDs not found in the genotype data
# ----------------------------------------------------------

# Identify phenotype IDs used by at least one genotype sample.
used_phenotype_ids <- unique(
    fam$MATCH_ID[!is.na(match_index)]
)

# Phenotype IDs not used by genotype samples are reported.
unmatched_phenotype_rows <- !(
    phenotype$MATCH_ID %in% used_phenotype_ids
)

write.table(
    data.frame(
        PHENOTYPE_ID = phenotype[[phenotype_id_column]][
            unmatched_phenotype_rows
        ],
        MATCH_ID = phenotype$MATCH_ID[
            unmatched_phenotype_rows
        ]
    ),
    file = opt$`unmatched-phenotype`,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)


# ----------------------------------------------------------
# 19. Write the final matching summary
# ----------------------------------------------------------

write_metric(
    c(
        mode = genotype_id_mode,
        id_normalization = opt$`id-normalization`,
        matching_id_class = paste(class(fam$MATCH_ID), collapse = ","),
        genotype_samples = nrow(fam),
        phenotype_rows = nrow(phenotype),
        matched = sum(!is.na(match_index)),
        genotype_only = sum(is.na(match_index)),
        phenotype_only = sum(unmatched_phenotype_rows),
        missing_reported_sex = sum(plink_sex == "0")
