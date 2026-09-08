#!/usr/bin/env Rscript

# Build the list of samples eligible for downstream genotype QC.
#
# Starting from all samples in a PLINK FAM file, this script:
#   1. Applies the intersection of any external keep files.
#   2. Applies the union of any external remove files.
#   3. Writes a PLINK-compatible FID/IID keep file.
#   4. Reports excluded samples and summary counts.


# ----------------------------------------------------------
# 1. Load packages
# ----------------------------------------------------------

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})


# ----------------------------------------------------------
# 2. Parse command-line options
# ----------------------------------------------------------

option_list <- list(
    make_option(
        "--fam",
        type = "character",
        help = "Input PLINK FAM file"
    ),
    make_option(
        "--keep-output",
        type = "character",
        help = "Output PLINK FID/IID keep file"
    ),
    make_option(
        "--exclusions-output",
        type = "character",
        help = "Output table containing excluded samples"
    ),
    make_option(
        "--report",
        type = "character",
        help = "Output summary report"
    ),
    make_option(
        "--keep-files",
        type = "character",
        default = "",
        help = "Semicolon-delimited list of external keep files"
    ),
    make_option(
        "--remove-files",
        type = "character",
        default = "",
        help = "Semicolon-delimited list of external remove files"
    )
)

options <- parse_args(
    OptionParser(option_list = option_list)
)


# ----------------------------------------------------------
# 3. Check required options
# ----------------------------------------------------------

required_options <- c(
    "fam",
    "keep-output",
    "exclusions-output",
    "report"
)

for (option_name in required_options) {
    if (is.null(options[[option_name]])) {
        stop(
            paste("Missing required option:", option_name),
            call. = FALSE
        )
    }
}


# ----------------------------------------------------------
# 4. Read and validate the FAM file
# ----------------------------------------------------------

fam <- fread(
    options$fam,
    header = FALSE,
    data.table = FALSE,
    colClasses = "character"
)

if (ncol(fam) < 6) {
    stop(
        "Malformed FAM: fewer than six columns",
        call. = FALSE
    )
}

# Retain the six standard PLINK FAM columns.
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
# 5. Parse configured file lists
# ----------------------------------------------------------

# Snakemake passes multiple paths as one semicolon-delimited string.
parse_file_list <- function(value) {
    missing_value <- (
        is.null(value) ||
        value == "" ||
        value == "NA" ||
        value == "None"
    )

    if (missing_value) {
        return(character())
    }

    paths <- strsplit(
        value,
        ";",
        fixed = TRUE
    )[[1]]

    paths <- trimws(paths)
    paths[nzchar(paths)]
}


# ----------------------------------------------------------
# 6. Read sample IDs from an external file
# ----------------------------------------------------------

read_sample_ids <- function(path) {
    if (!file.exists(path)) {
        stop(
            paste("Sample file not found:", path),
            call. = FALSE
        )
    }

    lines <- readLines(
        path,
        warn = FALSE
    )

    lines <- trimws(lines)

    # Ignore empty lines and comment lines.
    lines <- lines[
        nzchar(lines) &
        !grepl("^#", lines)
    ]

    if (!length(lines)) {
        return(character())
    }

    fields <- strsplit(
        lines,
        "[[:space:]]+"
    )

    matched_ids <- character()

    for (fields_in_line in fields) {
        if (!length(fields_in_line)) {
            next
        }

        first_field <- sub(
            "^#",
            "",
            fields_in_line[1]
        )

        # Ignore common FID/IID header lines.
        if (toupper(first_field) %in% c("FID", "IID")) {
            next
        }

        # For two-column PLINK files, prefer the second field (IID).
        if (
            length(fields_in_line) >= 2 &&
            fields_in_line[2] %in% fam$IID
        ) {
            matched_ids <- c(
                matched_ids,
                fields_in_line[2]
            )
        } else if (fields_in_line[1] %in% fam$IID) {
            # Also support one-column IID files.
            matched_ids <- c(
                matched_ids,
                fields_in_line[1]
            )
        }
    }

    unique(matched_ids)
}


# ----------------------------------------------------------
# 7. Apply external keep files
# ----------------------------------------------------------

# Start with every sample in the FAM file.
eligible_ids <- unique(fam$IID)

report_rows <- list()
report_index <- 1

keep_files <- parse_file_list(
    options$`keep-files`
)

# When multiple keep files are provided, a sample must occur
# in every file to remain eligible.
for (path in keep_files) {
    ids_in_file <- read_sample_ids(path)

    eligible_ids <- intersect(
        eligible_ids,
        ids_in_file
    )

    report_rows[[report_index]] <- data.frame(
        file = path,
        operation = "KEEP_INTERSECTION",
        matched_ids = length(ids_in_file)
    )

    report_index <- report_index + 1
}


# ----------------------------------------------------------
# 8. Apply external remove files
# ----------------------------------------------------------

remove_files <- parse_file_list(
    options$`remove-files`
)

excluded_ids <- character()

# When multiple remove files are provided, a sample is excluded
# when it occurs in at least one file.
for (path in remove_files) {
    ids_in_file <- read_sample_ids(path)

    excluded_ids <- union(
        excluded_ids,
        ids_in_file
    )

    report_rows[[report_index]] <- data.frame(
        file = path,
        operation = "REMOVE_UNION",
        matched_ids = length(ids_in_file)
    )

    report_index <- report_index + 1
}


# ----------------------------------------------------------
# 9. Build final sample tables
# ----------------------------------------------------------

final_ids <- setdiff(
    eligible_ids,
    excluded_ids
)

# PLINK-compatible FID/IID keep table.
kept_samples <- fam[
    fam$IID %in% final_ids,
    c("FID", "IID")
]

# Samples not retained for downstream QC.
excluded_samples <- fam[
    !fam$IID %in% final_ids,
    c("FID", "IID")
]

excluded_samples$REASON <- "external_or_eligibility_exclusion"


# ----------------------------------------------------------
# 10. Write sample outputs
# ----------------------------------------------------------

# PLINK --keep expects FID and IID without a header.
write.table(
    kept_samples,
    file = options$`keep-output`,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
)

# The exclusions report includes a header and exclusion reason.
write.table(
    excluded_samples,
    file = options$`exclusions-output`,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)


# ----------------------------------------------------------
# 11. Build and write the summary report
# ----------------------------------------------------------

if (length(report_rows)) {
    report <- do.call(
        rbind,
        report_rows
    )
} else {
    report <- data.frame(
        file = character(),
        operation = character(),
        matched_ids = integer()
    )
}

total_report <- data.frame(
    file = rep("__TOTAL__", 3),
    operation = c(
        "INPUT",
        "KEPT",
        "EXCLUDED"
    ),
    matched_ids = c(
        nrow(fam),
        nrow(kept_samples),
        nrow(excluded_samples)
    )
)

report <- rbind(
    report,
    total_report
)

write.table(
    report,
    file = options$report,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)
