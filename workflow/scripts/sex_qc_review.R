#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

# Command-line arguments --------------------------------------------------------

option_list <- list(
  make_option(
    "--sexcheck",
    type = "character",
    help = "PLINK sex-check file"
  ),
  make_option(
    "--fam",
    type = "character",
    help = "PLINK FAM file"
  ),
  make_option(
    "--threshold",
    type = "double",
    help = paste(
      "Fixed chrX F threshold:",
      "F <= threshold is female;",
      "F > threshold is male"
    )
  ),
  make_option(
    "--plot",
    type = "character",
    help = "Output diagnostic PDF"
  ),
  make_option(
    "--table",
    type = "character",
    help = "Output sample-level classification TSV"
  ),
  make_option(
    "--candidate-threshold",
    type = "character",
    help = "Output candidate-threshold YAML"
  ),
  make_option(
    "--summary",
    type = "character",
    help = "Output summary TSV"
  )
)

args <- parse_args(
  OptionParser(option_list = option_list)
)

# Validate arguments ------------------------------------------------------------

required_args <- c(
  "sexcheck",
  "fam",
  "threshold",
  "plot",
  "table",
  "candidate-threshold",
  "summary"
)

missing_args <- required_args[
  vapply(
    required_args,
    function(argument) {
      is.null(args[[argument]])
    },
    logical(1)
  )
]

if (length(missing_args)) {
  stop(
    "Missing required option(s): ",
    paste0("--", missing_args, collapse = ", "),
    call. = FALSE
  )
}

threshold <- args$threshold

if (!is.finite(threshold)) {
  stop(
    "--threshold must be a finite number.",
    call. = FALSE
  )
}

# Create output directories -----------------------------------------------------

output_paths <- c(
  args$plot,
  args$table,
  args[["candidate-threshold"]],
  args$summary
)

invisible(
  lapply(
    unique(dirname(output_paths)),
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  )
)

# Read PLINK sex-check results --------------------------------------------------

sexcheck <- fread(
  args$sexcheck,
  data.table = FALSE
)

if (!"IID" %in% names(sexcheck)) {
  stop(
    "Sex-check file is missing the IID column.",
    call. = FALSE
  )
}

f_column <- intersect(
  c("F", "XF", "X_F"),
  names(sexcheck)
)

if (!length(f_column)) {
  stop(
    "Cannot find a chrX F column. Available columns: ",
    paste(names(sexcheck), collapse = ", "),
    call. = FALSE
  )
}

sexcheck_id_columns <- if ("FID" %in% names(sexcheck)) {
  c("FID", "IID")
} else {
  "IID"
}

sexcheck <- sexcheck[
  ,
  c(sexcheck_id_columns, f_column[1]),
  drop = FALSE
]

names(sexcheck)[ncol(sexcheck)] <- "F"

sexcheck$F <- suppressWarnings(
  as.numeric(sexcheck$F)
)

if (anyDuplicated(sexcheck[sexcheck_id_columns])) {
  stop(
    "Sex-check file contains duplicate sample identifiers.",
    call. = FALSE
  )
}

# Read FAM file -----------------------------------------------------------------

fam <- fread(
  args$fam,
  header = FALSE,
  data.table = FALSE,
  colClasses = "character"
)

if (ncol(fam) < 6L) {
  stop(
    "Malformed FAM file: expected at least six columns.",
    call. = FALSE
  )
}

fam <- fam[, seq_len(6), drop = FALSE]

names(fam) <- c(
  "FID",
  "IID",
  "PAT",
  "MAT",
  "SEX",
  "PHENO"
)

if (anyDuplicated(fam[c("FID", "IID")])) {
  stop(
    "FAM file contains duplicate FID/IID pairs.",
    call. = FALSE
  )
}

# Join FAM and sex-check records ------------------------------------------------

fam$.FAM_ORDER <- seq_len(nrow(fam))

join_columns <- if ("FID" %in% names(sexcheck)) {
  c("FID", "IID")
} else {
  "IID"
}

samples <- merge(
  fam[, c("FID", "IID", "SEX", ".FAM_ORDER")],
  sexcheck,
  by = join_columns,
  all.x = TRUE,
  sort = FALSE
)

samples <- samples[
  order(samples$.FAM_ORDER),
  ,
  drop = FALSE
]

samples$.FAM_ORDER <- NULL

if (nrow(samples) != nrow(fam)) {
  stop(
    paste(
      "Internal join error:",
      "the classification table and FAM file",
      "contain different numbers of samples."
    ),
    call. = FALSE
  )
}

# Standardize reported sex ------------------------------------------------------

# PLINK sex codes:
#   1 = male
#   2 = female
#   0 = unknown

samples$REPORTED_SEX <- suppressWarnings(
  as.integer(samples$SEX)
)

samples$REPORTED_SEX[
  is.na(samples$REPORTED_SEX) |
    !samples$REPORTED_SEX %in% c(1L, 2L)
] <- 0L

# Classify genetic sex using the fixed threshold --------------------------------

# Classification:
#   F <= threshold = female
#   F > threshold  = male
#   missing F      = unknown

samples$GENETIC_SEX <- 0L

samples$GENETIC_SEX[
  is.finite(samples$F) &
    samples$F <= threshold
] <- 2L

samples$GENETIC_SEX[
  is.finite(samples$F) &
    samples$F > threshold
] <- 1L

# Add readable sex labels -------------------------------------------------------

sex_labels <- c(
  `0` = "unknown",
  `1` = "male",
  `2` = "female"
)

samples$REPORTED_SEX_LABEL <- unname(
  sex_labels[
    as.character(samples$REPORTED_SEX)
  ]
)

samples$GENETIC_SEX_LABEL <- unname(
  sex_labels[
    as.character(samples$GENETIC_SEX)
  ]
)

# Compare reported and genetic sex ---------------------------------------------

has_reported_sex <-
  samples$REPORTED_SEX %in% c(1L, 2L)

has_genetic_sex <-
  samples$GENETIC_SEX %in% c(1L, 2L)

has_f_value <- is.finite(samples$F)

samples$REPORTED_GENETIC_SEX_MATCH <-
  has_reported_sex &
  has_genetic_sex &
  samples$REPORTED_SEX == samples$GENETIC_SEX

# Assign QC status --------------------------------------------------------------

samples$STATUS <- "OK"

samples$STATUS[
  !has_f_value
] <- "MISSING_GENETIC_SEX"

samples$STATUS[
  has_reported_sex &
    has_genetic_sex &
    samples$REPORTED_SEX != samples$GENETIC_SEX
] <- "DISCORDANT"

# Unknown reported sex takes precedence.
samples$STATUS[
  samples$REPORTED_SEX == 0L
] <- "REPORTED_SEX_UNKNOWN"

# With one threshold, every finite F value is classified.
# Therefore, there is no ambiguous interval.
samples$EXCLUDE <- samples$STATUS == "DISCORDANT"

# Write the sample-level comparison table --------------------------------------

output_columns <- c(
  "FID",
  "IID",
  "F",
  "REPORTED_SEX",
  "REPORTED_SEX_LABEL",
  "GENETIC_SEX",
  "GENETIC_SEX_LABEL",
  "REPORTED_GENETIC_SEX_MATCH",
  "STATUS",
  "EXCLUDE"
)

write.table(
  samples[, output_columns, drop = FALSE],
  file = args$table,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  na = "NA"
)

# Write fixed candidate threshold ----------------------------------------------

writeLines(
  c(
    "# Fixed candidate threshold; manual approval is required.",
    sprintf(
      "sex_f_threshold: %.8g",
      threshold
    )
  ),
  args[["candidate-threshold"]]
)

# Create diagnostic plots -------------------------------------------------------

reported_female_f <- samples$F[
  samples$REPORTED_SEX == 2L &
    is.finite(samples$F)
]

reported_male_f <- samples$F[
  samples$REPORTED_SEX == 1L &
    is.finite(samples$F)
]

finite_f <- samples$F[
  is.finite(samples$F)
]

pdf(
  args$plot,
  width = 9,
  height = 7
)

if (length(finite_f)) {
  # First plot: distribution of all finite F values.
  hist(
    finite_f,
    breaks = 80,
    col = "grey80",
    border = "white",
    main = "chrX F distribution",
    xlab = "chrX inbreeding coefficient (F)"
  )

  abline(
    v = threshold,
    col = "black",
    lty = 2,
    lwd = 2
  )

  legend(
    "topright",
    legend = sprintf(
      "Fixed threshold = %.4g",
      threshold
    ),
    col = "black",
    lty = 2,
    lwd = 2,
    bty = "n"
  )

  # Second plot: F values grouped by reported sex.
  female_plot_values <- if (length(reported_female_f)) {
    reported_female_f
  } else {
    NA_real_
  }

  male_plot_values <- if (length(reported_male_f)) {
    reported_male_f
  } else {
    NA_real_
  }

  boxplot(
    list(
      "Reported female" = female_plot_values,
      "Reported male" = male_plot_values
    ),
    col = c(
      rgb(1, 0, 0, 0.35),
      rgb(0, 0, 1, 0.35)
    ),
    ylab = "chrX inbreeding coefficient (F)",
    main = "chrX F by reported sex"
  )

  abline(
    h = threshold,
    col = "black",
    lty = 2,
    lwd = 2
  )
} else {
  plot.new()

  title(
    "chrX F distribution"
  )

  text(
    0.5,
    0.5,
    "No finite chrX F values available"
  )
}

dev.off()

# Write summary -----------------------------------------------------------------

summary_table <- data.frame(
  metric = c(
    "sex_f_threshold",
    "n_samples",
    "n_samples_with_f",
    "n_genetic_female",
    "n_genetic_male",
    "n_reported_genetic_match",
    "n_discordant",
    "n_missing_genetic_sex",
    "n_reported_sex_unknown",
    "n_excluded"
  ),
  value = c(
    threshold,
    nrow(samples),
    sum(has_f_value),
    sum(samples$GENETIC_SEX == 2L),
    sum(samples$GENETIC_SEX == 1L),
    sum(samples$REPORTED_GENETIC_SEX_MATCH),
    sum(samples$STATUS == "DISCORDANT"),
    sum(samples$STATUS == "MISSING_GENETIC_SEX"),
    sum(samples$STATUS == "REPORTED_SEX_UNKNOWN"),
    sum(samples$EXCLUDE)
  ),
  stringsAsFactors = FALSE
)

write.table(
  summary_table,
  file = args$summary,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)