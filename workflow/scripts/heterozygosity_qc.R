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
        "--het",
        type = "character",
        help = "PLINK .het input file"
    ),
    make_option(
        "--n-sd",
        type = "double",
        default = 3,
        help = paste(
            "Number of standard deviations used to define",
            "heterozygosity outliers [default: %default]"
        )
    ),
    make_option(
        "--plot",
        type = "character",
        help = "Output PNG histogram"
    ),
    make_option(
        "--suggested-exclusions",
        type = "character",
        help = "Output PLINK-compatible sample exclusion file"
    ),
    make_option(
        "--table",
        type = "character",
        help = "Output table containing all samples"
    ),
    make_option(
        "--summary",
        type = "character",
        help = "Output QC summary table"
    )
)

opt <- parse_args(
    OptionParser(option_list = option_list)
)


# -------------------------------------------------------------------------
# 2. Validate arguments
# -------------------------------------------------------------------------

required_options <- c(
    "het",
    "plot",
    "suggested-exclusions",
    "table",
    "summary"
)

for (option_name in required_options) {
    if (
        is.null(opt[[option_name]]) ||
        !nzchar(opt[[option_name]])
    ) {
        stop(
            paste("Missing required option:", option_name),
            call. = FALSE
        )
    }
}

if (!file.exists(opt$het)) {
    stop(
        paste("Input .het file not found:", opt$het),
        call. = FALSE
    )
}

if (
    !is.finite(opt$`n-sd`) ||
    opt$`n-sd` <= 0
) {
    stop(
        "--n-sd must be a positive number.",
        call. = FALSE
    )
}


# -------------------------------------------------------------------------
# 3. Create output directories
# -------------------------------------------------------------------------

output_paths <- c(
    opt$plot,
    opt$`suggested-exclusions`,
    opt$table,
    opt$summary
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
# 4. Read the PLINK heterozygosity report
# -------------------------------------------------------------------------

het_data <- fread(
    opt$het,
    data.table = FALSE
)

if (nrow(het_data) == 0) {
    stop(
        "The input .het file contains no samples.",
        call. = FALSE
    )
}

# PLINK may name the family-ID column either FID or #FID.
if ("#FID" %in% names(het_data)) {
    fid_column <- "#FID"
} else if ("FID" %in% names(het_data)) {
    fid_column <- "FID"
} else {
    # Some PLINK datasets may not contain family IDs.
    fid_column <- NA_character_
}

required_columns <- c("IID", "F")

missing_columns <- setdiff(
    required_columns,
    names(het_data)
)

if (length(missing_columns) > 0) {
    stop(
        paste(
            "Unexpected .het format. Missing columns:",
            paste(missing_columns, collapse = ", ")
        ),
        call. = FALSE
    )
}


# -------------------------------------------------------------------------
# 5. Calculate heterozygosity thresholds
# -------------------------------------------------------------------------

het_data$F <- suppressWarnings(
    as.numeric(het_data$F)
)

valid_f <- is.finite(het_data$F)

if (!any(valid_f)) {
    stop(
        "The .het file contains no valid finite F values.",
        call. = FALSE
    )
}

mean_f <- mean(
    het_data$F[valid_f]
)

sd_f <- sd(
    het_data$F[valid_f]
)

if (!is.finite(sd_f)) {
    stop(
        "The standard deviation of F could not be calculated.",
        call. = FALSE
    )
}

lower_threshold <- mean_f - opt$`n-sd` * sd_f
upper_threshold <- mean_f + opt$`n-sd` * sd_f


# -------------------------------------------------------------------------
# 6. Identify heterozygosity outliers
# -------------------------------------------------------------------------
#
# A sample is classified as an outlier when its F value is below or above:
#
#   mean(F) ± n_sd × sd(F)
#
# Samples with missing or non-finite F values are not classified as
# heterozygosity outliers.
# -------------------------------------------------------------------------

het_data$HET_OUTLIER <- (
    valid_f &
    (
        het_data$F < lower_threshold |
        het_data$F > upper_threshold
    )
)

outliers <- het_data[
    het_data$HET_OUTLIER,
    ,
    drop = FALSE
]


# -------------------------------------------------------------------------
# 7. Write the complete per-sample table
# -------------------------------------------------------------------------

write.table(
    het_data,
    file = opt$table,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
)


# -------------------------------------------------------------------------
# 8. Write the PLINK-compatible exclusion file
# -------------------------------------------------------------------------
#
# PLINK --remove accepts a two-column sample file containing FID and IID.
# If FID is unavailable, IID is used for both identifiers.
#
# rep(..., length = 0) ensures that an empty but valid exclusion file is
# generated when no outliers are found.
# -------------------------------------------------------------------------

if (is.na(fid_column)) {
    exclusion_fid <- outliers$IID
} else {
    exclusion_fid <- outliers[[fid_column]]
}

exclusion_table <- data.frame(
    FID = exclusion_fid,
    IID = outliers$IID,
    stringsAsFactors = FALSE
)

write.table(
    exclusion_table,
    file = opt$`suggested-exclusions`,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
)


# -------------------------------------------------------------------------
# 9. Create the heterozygosity histogram
# -------------------------------------------------------------------------

png(
    filename = opt$plot,
    width = 1200,
    height = 750,
    res = 150
)

hist(
    het_data$F[valid_f],
    breaks = 80,
    main = "Autosomal heterozygosity F",
    xlab = "Inbreeding coefficient (F)",
    col = "steelblue",
    border = "white"
)

abline(
    v = c(
        lower_threshold,
        upper_threshold
    ),
    col = "red",
    lty = 2,
    lwd = 2
)

legend(
    "topright",
    legend = paste0(
        "Mean \u00b1 ",
        opt$`n-sd`,
        " SD"
    ),
    col = "red",
    lty = 2,
    lwd = 2,
    bty = "n"
)

dev.off()


# -------------------------------------------------------------------------
# 10. Write the QC summary
# -------------------------------------------------------------------------

summary_table <- data.frame(
    metric = c(
        "n_samples",
        "n_valid_F",
        "n_missing_F",
        "mean_F",
        "sd_F",
        "n_sd",
        "lower_threshold",
        "upper_threshold",
        "n_outliers"
    ),
    value = c(
        nrow(het_data),
        sum(valid_f),
        sum(!valid_f),
        mean_f,
        sd_f,
        opt$`n-sd`,
        lower_threshold,
        upper_threshold,
        nrow(outliers)
    ),
    stringsAsFactors = FALSE
)

write.table(
    summary_table,
    file = opt$summary,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
)