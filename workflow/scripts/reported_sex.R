#!/usr/bin/env Rscript

# Prepare a PLINK-compatible sex update file by matching samples
# from a FAM file with phenotype metadata.
#
# The script:
#   1. Reads the PLINK FAM file.
#   2. Reads phenotype metadata from DTA, CSV, CSV.GZ, TSV, or TXT.
#   3. Normalizes genotype and phenotype sample identifiers.
#   4. Matches genotype samples to phenotype metadata.
#   5. Converts reported sex values to PLINK codes:
#        1 = male
#        2 = female
#        0 = unknown
#   6. Produces a PLINK --update-sex file.
#   7. Reports unmatched genotype and metadata samples.
#   8. Produces a summary report.
#
# The original FAM and metadata files are not modified.


# ----------------------------------------------------------
# 1. Load required R packages
# ----------------------------------------------------------

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
    library(haven)
})


# ----------------------------------------------------------
# 2. Define command-line arguments
# ----------------------------------------------------------

option_list <- list(
    make_option(
        "--fam",
        type = "character",
        help = "Input PLINK FAM file"
    ),
    make_option(
        "--metadata",
        type = "character",
        help = paste(
            "Input phenotype metadata file.",
            "Supported formats: DTA, CSV, CSV.GZ, TSV, and TXT."
        )
    ),
    make_option(
        "--phenotype-id-col",
        type = "character",
        help = "Metadata column containing the sample identifier"
    ),
    make_option(
        "--sex-col",
        type = "character",
        help = "Metadata column containing reported sex"
    ),
    make_option(
        "--genotype-id-mode",
        type = "character",
        default = "direct",
        help = paste(
            "Method used to obtain the genotype matching ID:",
            "fam, direct, or regex [default: %default]"
        )
    ),
    make_option(
        "--genotype-id-regex",
        type = "character",
        default = "NA",
        help = paste(
            "Regular expression used to extract the matching ID",
            "from the FAM IID when genotype-id-mode=regex"
        )
    ),
    make_option(
        "--id-normalization",
        type = "character",
        default = "string",
        help = paste(
            "Method used to normalize IDs:",
            "string or integer_string [default: %default]"
        )
    ),
    make_option(
        "--male-values",
        type = "character",
        default = "1,M,Male",
        help = paste(
            "Comma-separated metadata values interpreted as male",
            "[default: %default]"
        )
    ),
    make_option(
        "--female-values",
        type = "character",
        default = "2,F,Female",
        help = paste(
            "Comma-separated metadata values interpreted as female",
            "[default: %default]"
        )
    ),
    make_option(
        "--update-sex",
        type = "character",
        help = "Output PLINK --update-sex file"
    ),
    make_option(
        "--report",
        type = "character",
        help = "Output summary report"
    ),
    make_option(
        "--unmatched-genotypes",
        type = "character",
        help = "Output file listing genotype samples absent from metadata"
    ),
    make_option(
        "--unmatched-metadata",
        type = "character",
        help = "Output file listing metadata IDs absent from the FAM"
    )
)

opt <- parse_args(
    OptionParser(option_list = option_list)
)


# ----------------------------------------------------------
# 3. Check required command-line arguments
# ----------------------------------------------------------

required_options <- c(
    "fam",
    "metadata",
    "phenotype-id-col",
    "sex-col",
    "update-sex",
    "report",
    "unmatched-genotypes",
    "unmatched-metadata"
)

for (option_name in required_options) {
    if (is.null(opt[[option_name]])) {
        stop(
            paste("Missing required option:", option_name),
            call. = FALSE
        )
    }
}


# ----------------------------------------------------------
# 4. Read and validate the PLINK FAM file
# ----------------------------------------------------------

# Read every FAM column as character. This preserves sample IDs
# containing leading zeros, such as "00123".
fam <- fread(
    opt$fam,
    header = FALSE,
    data.table = FALSE,
    colClasses = "character"
)

# A standard PLINK FAM file must contain at least six fields:
# FID IID PAT MAT SEX PHENOTYPE
if (ncol(fam) < 6) {
    stop(
        "FAM has fewer than 6 fields",
        call. = FALSE
    )
}

# Retain the first six columns and assign standard PLINK names.
fam <- fam[, 1:6]

names(fam) <- c(
    "FID",
    "IID",
    "PAT",
    "MAT",
    "SEX",
    "PHENO"
)


# ----------------------------------------------------------
# 5. Define a function for writing summary metrics
# ----------------------------------------------------------

write_metric <- function(values) {
    report <- data.frame(
        metric = names(values),
        value = unname(values),
        check.names = FALSE
    )

    write.table(
        report,
        file = opt$report,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )
}


# ----------------------------------------------------------
# 6. Optionally use sex values directly from the FAM file
# ----------------------------------------------------------

# External metadata matching is skipped when:
#   - --metadata is set to "NA"; or
#   - --genotype-id-mode is set to "fam".
#
# In this mode, the existing FAM sex values are copied directly
# to the PLINK update-sex file.
if (
    opt$metadata == "NA" ||
    tolower(opt$`genotype-id-mode`) == "fam"
) {
    # Write FID, IID, and existing FAM sex without a header.
    write.table(
        fam[, c("FID", "IID", "SEX")],
        file = opt$`update-sex`,
        sep = "\t",
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE
    )

    # No genotype samples are unmatched because metadata matching
    # is not performed in this mode.
    write.table(
        data.frame(
            FID = character(),
            IID = character()
        ),
        file = opt$`unmatched-genotypes`,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )

    # Create an empty unmatched-metadata report.
    write.table(
        data.frame(
            METADATA_ID = character()
        ),
        file = opt$`unmatched-metadata`,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )

    # Record summary statistics.
    write_metric(
        c(
            mode = "fam",
            genotype_samples = nrow(fam),
            matched = nrow(fam),
            genotype_only = 0,
            metadata_only = 0,
            missing_reported_sex = sum(
                !fam$SEX %in% c("1", "2")
            )
        )
    )

    quit(save = "no")
}


# ----------------------------------------------------------
# 7. Define a function for reading phenotype metadata
# ----------------------------------------------------------

read_metadata <- function(path) {
    lowercase_path <- tolower(path)

    # Read Stata files.
    if (grepl("\\.dta$", lowercase_path)) {
        return(
            as.data.frame(read_dta(path))
        )
    }

    # Read uncompressed or gzip-compressed CSV files.
    if (grepl("\\.csv(\\.gz)?$", lowercase_path)) {
        return(
            fread(
                path,
                data.table = FALSE
            )
        )
    }

    # Automatically detect the delimiter for other text formats,
    # including TSV and TXT files.
    fread(
        path,
        data.table = FALSE,
        sep = "auto"
    )
}


# ----------------------------------------------------------
# 8. Read and validate the phenotype metadata
# ----------------------------------------------------------

metadata <- read_metadata(opt$metadata)

phenotype_id_column <- opt$`phenotype-id-col`
sex_column <- opt$`sex-col`

# Ensure that the requested sample ID and sex columns exist.
if (
    !phenotype_id_column %in% names(metadata) ||
    !sex_column %in% names(metadata)
) {
    stop(
        paste(
            "Phenotype ID or sex column missing from metadata.",
            "Available columns:",
            paste(names(metadata), collapse = ", ")
        ),
        call. = FALSE
    )
}


# ----------------------------------------------------------
# 9. Define a function for normalizing sample IDs
# ----------------------------------------------------------

normalize_id <- function(ids) {
    # Convert IDs to character and remove surrounding whitespace.
    ids <- trimws(as.character(ids))

    # Treat missing and empty identifiers as NA.
    ids[is.na(ids) | ids == ""] <- NA_character_

    # Optionally convert numeric-looking IDs to integer strings.
    #
    # Examples:
    #   "00123" -> "123"
    #   "123.0" -> "123"
    #
    # This should only be used when leading zeros are not meaningful.
    if (opt$`id-normalization` == "integer_string") {
        numeric_ids <- suppressWarnings(
            as.numeric(ids)
        )

        valid_numeric_ids <- !is.na(numeric_ids)

        ids[valid_numeric_ids] <- as.character(
            as.integer(numeric_ids[valid_numeric_ids])
        )
    }

    ids
}


# ----------------------------------------------------------
# 10. Create matching IDs for genotype samples
# ----------------------------------------------------------

# By default, use the FAM IID directly as the matching ID.
fam$MATCH_ID <- normalize_id(fam$IID)

genotype_id_mode <- tolower(
    opt$`genotype-id-mode`
)

if (genotype_id_mode == "regex") {
    genotype_id_regex <- opt$`genotype-id-regex`

    if (genotype_id_regex == "NA") {
        stop(
            paste(
                "genotype-id-regex is required",
                "when genotype-id-mode=regex"
            ),
            call. = FALSE
        )
    }

    # Extract a matching ID from one FAM IID.
    #
    # When the regular expression contains a capture group,
    # the first captured group is returned. Otherwise, the
    # complete regular-expression match is returned.
    extract_one_id <- function(sample_id) {
        match_information <- regexec(
            genotype_id_regex,
            sample_id,
            perl = TRUE
        )

        extracted_values <- regmatches(
            sample_id,
            match_information
        )[[1]]

        # Return NA when the regular expression does not match.
        if (!length(extracted_values)) {
            return(NA_character_)
        }

        # Prefer the first capture group when one is available.
        if (length(extracted_values) >= 2) {
            return(extracted_values[2])
        }

        extracted_values[1]
    }

    fam$MATCH_ID <- normalize_id(
        vapply(
            fam$IID,
            extract_one_id,
            character(1)
        )
    )
} else if (genotype_id_mode != "direct") {
    stop(
        paste(
            "genotype-id-mode must be",
            "fam, direct, or regex"
        ),
        call. = FALSE
    )
}


# ----------------------------------------------------------
# 11. Create normalized matching IDs for the metadata
# ----------------------------------------------------------

metadata$MATCH_ID <- normalize_id(
    metadata[[phenotype_id_column]]
)


# ----------------------------------------------------------
# 12. Check for duplicated phenotype IDs
# ----------------------------------------------------------

# Missing metadata IDs are excluded from the duplicate check.
nonmissing_metadata_ids <- metadata$MATCH_ID[
    !is.na(metadata$MATCH_ID)
]

# Duplicated normalized IDs would make genotype-to-metadata
# matching ambiguous, so the script stops if they are detected.
if (anyDuplicated(nonmissing_metadata_ids)) {
    duplicated_ids <- unique(
        nonmissing_metadata_ids[
            duplicated(nonmissing_metadata_ids) |
            duplicated(
                nonmissing_metadata_ids,
                fromLast = TRUE
            )
        ]
    )

    stop(
        paste(
            "Duplicate phenotype IDs found after ID normalization, e.g.",
            paste(
                head(duplicated_ids, 10),
                collapse = ", "
            )
        ),
        call. = FALSE
    )
}


# ----------------------------------------------------------
# 13. Define the recognized male and female values
# ----------------------------------------------------------

# Convert the comma-separated command-line values into
# lowercase vectors for case-insensitive comparison.
male_values <- tolower(
    trimws(
        strsplit(
            opt$`male-values`,
            ",",
            fixed = TRUE
        )[[1]]
    )
)

female_values <- tolower(
    trimws(
        strsplit(
            opt$`female-values`,
            ",",
            fixed = TRUE
        )[[1]]
    )
)


# ----------------------------------------------------------
# 14. Convert reported sex to PLINK codes
# ----------------------------------------------------------

sex_to_plink <- function(values) {
    standardized_values <- tolower(
        trimws(
            as.character(values)
        )
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

metadata$PLINK_SEX <- sex_to_plink(
    metadata[[sex_column]]
)


# ----------------------------------------------------------
# 15. Match genotype samples to phenotype metadata
# ----------------------------------------------------------

# For every genotype sample, return the position of its MATCH_ID
# in the metadata. An unmatched sample receives NA.
match_index <- match(
    fam$MATCH_ID,
    metadata$MATCH_ID
)

# Use the metadata sex when a match exists. Assign PLINK sex code
# 0 when the genotype sample is not present in the metadata.
plink_sex <- ifelse(
    is.na(match_index),
    "0",
    metadata$PLINK_SEX[match_index]
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
# 17. Write genotype samples not found in the metadata
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
# 18. Write metadata IDs not found in the genotype data
# ----------------------------------------------------------

# Identify metadata IDs used by at least one genotype sample.
used_metadata_ids <- unique(
    fam$MATCH_ID[!is.na(match_index)]
)

# Metadata IDs not used by genotype samples are reported.
unmatched_metadata_rows <- !(
    metadata$MATCH_ID %in% used_metadata_ids
)

write.table(
    data.frame(
        METADATA_ID = metadata[
            [phenotype_id_column]
        ][unmatched_metadata_rows],
        MATCH_ID = metadata$MATCH_ID[
            unmatched_metadata_rows
        ]
    ),
    file = opt$`unmatched-metadata`,
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
        genotype_samples = nrow(fam),
        metadata_rows = nrow(metadata),
        matched = sum(!is.na(match_index)),
        genotype_only = sum(is.na(match_index)),
        metadata_only = sum(unmatched_metadata_rows),
        missing_reported_sex = sum(plink_sex == "0")
    )
)