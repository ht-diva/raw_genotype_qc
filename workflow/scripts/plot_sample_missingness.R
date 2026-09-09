#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})

options <- list(
    make_option("--smiss", type = "character"),
    make_option("--threshold", type = "double"),
    make_option("--plot", type = "character"),
    make_option("--summary", type = "character")
)

args <- parse_args(
    OptionParser(option_list = options)
)

data <- fread(
    args$smiss,
    data.table = FALSE
)

missingness_columns <- names(data)[
    grepl("^F_MISS", toupper(names(data)))
]

if (!length(missingness_columns)) {
    stop(
        paste(
            "No F_MISS column in",
            args$smiss
        )
    )
}

missingness <- as.numeric(
    data[[missingness_columns[1]]]
)

dir.create(
    dirname(args$plot),
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    dirname(args$summary),
    recursive = TRUE,
    showWarnings = FALSE
)

png(
    args$plot,
    width = 1200,
    height = 750,
    res = 150
)

hist(
    missingness,
    breaks = 80,
    main = "Sample missingness",
    xlab = "Sample missing genotype rate"
)

abline(
    v = args$threshold,
    lty = 2
)

dev.off()

summary <- data.frame(
    metric = c(
        "n_samples",
        "threshold",
        "n_above_threshold",
        "mean_missing_rate",
        "max_missing_rate"
    ),
    value = c(
        nrow(data),
        args$threshold,
        sum(
            missingness > args$threshold,
            na.rm = TRUE
        ),
        mean(
            missingness,
            na.rm = TRUE
        ),
        max(
            missingness,
            na.rm = TRUE
        )
    )
)

write.table(
    summary,
    args$summary,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)