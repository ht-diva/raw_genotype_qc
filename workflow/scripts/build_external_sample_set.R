#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})


# -------------------------------------------------------------------------
# 1. Parse command-line arguments
# -------------------------------------------------------------------------

option_list <- list(
    make_option(
        "--fam",
        type = "character",
        help = "Input PLINK FAM file"
    ),
    make_option(
        "--keep-output",
        type = "character",
        help = "Output PLINK keep file"
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
        help = "Semicolon-separated list of sample keep files"
    ),
    make_option(
        "--remove-files",
        type = "character",
        default = "",
        help = "Semicolon-separated list of sample removal files"
    )
)

opt <- parse_args(
    OptionParser(option_list = option_list)
)


# -------------------------------------------------------------------------
# 2. Check required arguments
# -------------------------------------------------------------------------

required_options <- c(
    "fam",
    "keep-output",
    "exclusions-output",
    "report"
)

for (option_name in required_options) {
    if (is.null(opt[[option_name]])) {
        stop(
            paste("Missing required option:", option_name),
            call. = FALSE
        )
    }
}


# -------------------------------------------------------------------------
# 3. Read the PLINK FAM file
# -------------------------------------------------------------------------

fam <- fread(
    opt$fam,
    header = FALSE,
    data.table = FALSE,
    colClasses = "character"
)

if (ncol(fam) < 6) {
    stop(
        "Malformed FAM file: fewer than six columns were found.",
        call. = FALSE
    )
}

# Retain only the six standard PLINK FAM columns.
fam <- fam[, 1:6, drop = FALSE]

names(fam) <- c(
    "FID",
    "IID",
    "PAT",
    "MAT",
    "SEX",
    "PHENO"
)

if (nrow(fam) == 0) {
    stop(
        "The input FAM file contains no samples.",
        call. = FALSE
    )
}


# -------------------------------------------------------------------------
# 4. Parse semicolon-separated file lists
# -------------------------------------------------------------------------

parse_files <- function(value) {
    if (
        is.null(value) ||
        length(value) == 0 ||
        is.na(value) ||
        trimws(value) == "" ||
        tolower(trimws(value)) %in% c("na", "none")
    ) {
        return(character(0))
    }

    paths <- strsplit(
        value,
        split = ";",
        fixed = TRUE
    )[[1]]

    paths <- trimws(paths)
    paths[nzchar(paths)]
}

keep_files <- parse_files(opt$`keep-files`)
remove_files <- parse_files(opt$`remove-files`)


# -------------------------------------------------------------------------
# 5. Read sample IDs from an external keep/remove file
# -------------------------------------------------------------------------
#
# Supported formats:
#
#   IID
#
# or:
#
#   FID IID
#
# Empty lines, comments beginning with "#", and common headers are ignored.
# Only sample IDs present in the input FAM file are returned.
# -------------------------------------------------------------------------

read_sample_ids <- function(path, fam_data) {
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

    # Remove empty lines and comment lines.
    lines <- lines[
        nzchar(lines) &
        !grepl("^#", lines)
    ]

    if (length(lines) == 0) {
        return(character(0))
    }

    fields <- strsplit(
        lines,
        split = "[[:space:]]+"
    )

    sample_ids <- character(0)

    for (values in fields) {
        if (length(values) == 0) {
            next
        }

        first_value <- toupper(
            sub("^#", "", values[1])
        )

        # Skip common header rows.
        if (first_value %in% c("FID", "IID")) {
            next
        }

        # For two-column files, preferentially interpret the second
        # column as IID. Otherwise, interpret the first column as IID.
        if (
            length(values) >= 2 &&
            values[2] %in% fam_data$IID
        ) {
            sample_ids <- c(
                sample_ids,
                values[2]
            )
        } else if (values[1] %in% fam_data$IID) {
            sample_ids <- c(
                sample_ids,
                values[1]
            )
        }
    }

    unique(sample_ids)
}


# -------------------------------------------------------------------------
# 6. Apply external keep files
# -------------------------------------------------------------------------
#
# If multiple keep files are supplied, their intersection is used:
# a sample must occur in every keep file to remain eligible.
#
# If no keep file is supplied, all samples initially remain eligible.
# -------------------------------------------------------------------------

eligible_ids <- unique(fam$IID)

report_rows <- list()
report_index <- 1L

for (keep_file in keep_files) {
    file_ids <- read_sample_ids(
        keep_file,
        fam
    )

    eligible_ids <- intersect(
        eligible_ids,
        file_ids
    )

    report_rows[[report_index]] <- data.frame(
        file = keep_file,
        operation = "KEEP_INTERSECTION",
        matched_ids = length(file_ids),
        stringsAsFactors = FALSE
    )

    report_index <- report_index + 1L
}


# -------------------------------------------------------------------------
# 7. Apply external removal files
# -------------------------------------------------------------------------
#
# If multiple removal files are supplied, their union is used:
# a sample appearing in any removal file is excluded.
# -------------------------------------------------------------------------

remove_ids <- character(0)

for (remove_file in remove_files) {
    file_ids <- read_sample_ids(
        remove_file,
        fam
    )

    remove_ids <- union(
        remove_ids,
        file_ids
    )

    report_rows[[report_index]] <- data.frame(
        file = remove_file,
        operation = "REMOVE_UNION",
        matched_ids = length(file_ids),
        stringsAsFactors = FALSE
    )

    report_index <- report_index + 1L
}


# -------------------------------------------------------------------------
# 8. Determine the final eligible and excluded samples
# -------------------------------------------------------------------------

final_ids <- setdiff(
    eligible_ids,
    remove_ids
)

kept_samples <- fam[
    fam$IID %in% final_ids,
    c("FID", "IID"),
    drop = FALSE
]

excluded_samples <- fam[
    !fam$IID %in% final_ids,
    c("FID", "IID"),
    drop = FALSE
]

# rep(..., nrow(...)) works correctly even when there are no exclusions.
excluded_samples$REASON <- rep(
    "external_or_eligibility_exclusion",
    nrow(excluded_samples)
)


# -------------------------------------------------------------------------
# 9. Create output directories
# -------------------------------------------------------------------------

output_paths <- c(
    opt$`keep-output`,
    opt$`exclusions-output`,
    opt$report
)

for (output_path in output_paths) {
    output_directory <- dirname(output_path)

    if (
        !dir.exists(output_directory) &&
        !dir.create(
            output_directory,
            recursive = TRUE,
            showWarnings = FALSE
        )
    ) {
        stop(
            paste(
                "Could not create output directory:",
                output_directory
            ),
            call. = FALSE
        )
    }
}


# -------------------------------------------------------------------------
# 10. Write the PLINK keep file
# -------------------------------------------------------------------------
#
# PLINK expects two columns without a header:
#
#   FID IID
# -------------------------------------------------------------------------

write.table(
    kept_samples,
    file = opt$`keep-output`,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
)


# -------------------------------------------------------------------------
# 11. Write the exclusions table
# -------------------------------------------------------------------------

write.table(
    excluded_samples,
    file = opt$`exclusions-output`,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
)


# -------------------------------------------------------------------------
# 12. Build and write the summary report
# -------------------------------------------------------------------------

if (length(report_rows) > 0) {
    report <- do.call(
        rbind,
        report_rows
    )
} else {
    report <- data.frame(
        file = character(0),
        operation = character(0),
        matched_ids = integer(0),
        stringsAsFactors = FALSE
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
    ),
    stringsAsFactors = FALSE
)

report <- rbind(
    report,
    total_report
)

write.table(
    report,
    file = opt$report,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
)