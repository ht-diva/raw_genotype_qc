#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

# Command-line options ----------------------------------------------------------

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
    "--female-max-f",
    type = "double",
    help = "Maximum chrX F for genetic females"
  ),
  make_option(
    "--male-min-f",
    type = "double",
    help = "Minimum chrX F for genetic males"
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
    "--candidate-thresholds",
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
  "female-max-f",
  "male-min-f",
  "plot",
  "table",
  "candidate-thresholds",
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

female_max_f <- args[["female-max-f"]]
male_min_f <- args[["male-min-f"]]

if (
  !is.finite(female_max_f) ||
  !is.finite(male_min_f)
) {
  stop(
    "Sex-QC thresholds must be finite numbers.",
    call. = FALSE
  )
}

if (female_max_f >= male_min_f) {
  stop(
    "female-max-f must be lower than male-min-f.",
    call. = FALSE
  )
}

# Create output directories -----------------------------------------------------

output_paths <- c(
  args$plot,
  args$table,
  args[["candidate-thresholds"]],
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

# Read FAM ----------------------------------------------------------------------

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

# Join FAM and sex-check results ------------------------------------------------

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

# PLINK codes:
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

# Classify genetic sex ----------------------------------------------------------

# F <= female_max_f = female
# F >= male_min_f   = male
# Values between the thresholds remain ambiguous.

samples$GENETIC_SEX <- 0L

samples$GENETIC_SEX[
  is.finite(samples$F) &
    samples$F <= female_max_f
] <- 2L

samples$GENETIC_SEX[
  is.finite(samples$F) &
    samples$F >= male_min_f
] <- 1L

# Add readable labels -----------------------------------------------------------

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

has_f_value <- is.finite(samples$F)

has_reported_sex <-
  samples$REPORTED_SEX %in% c(1L, 2L)

has_genetic_sex <-
  samples$GENETIC_SEX %in% c(1L, 2L)

samples$REPORTED_GENETIC_SEX_MATCH <-
  has_reported_sex &
  has_genetic_sex &
  samples$REPORTED_SEX == samples$GENETIC_SEX

# Assign status -----------------------------------------------------------------

samples$STATUS <- "OK"

samples$STATUS[
  !has_f_value
] <- "MISSING_GENETIC_SEX"

samples$STATUS[
  has_f_value &
    !has_genetic_sex
] <- "AMBIGUOUS_GENETIC_SEX"

samples$STATUS[
  has_reported_sex &
    has_genetic_sex &
    samples$REPORTED_SEX != samples$GENETIC_SEX
] <- "DISCORDANT"

# Unknown reported sex takes precedence.
samples$STATUS[
  samples$REPORTED_SEX == 0L
] <- "REPORTED_SEX_UNKNOWN"

# Candidate exclusions:
#   - reported/genetic sex discordance
#   - F inside the ambiguous interval
samples$EXCLUDE <-
  samples$STATUS %in% c(
    "DISCORDANT",
    "AMBIGUOUS_GENETIC_SEX"
  )

# Write sample-level classification --------------------------------------------

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

# Write candidate thresholds ---------------------------------------------------

# These values come directly from config.yaml.
# No empirical threshold is calculated.

writeLines(
  c(
    "# Fixed candidate thresholds; manual approval required.",
    sprintf(
      "female_max_f: %.8g",
      female_max_f
    ),
    sprintf(
      "male_min_f: %.8g",
      male_min_f
    )
  ),
  args[["candidate-thresholds"]]
)

# Create diagnostic plots -------------------------------------------------------

finite_f <- samples$F[
  is.finite(samples$F)
]

reported_female_f <- samples$F[
  samples$REPORTED_SEX == 2L &
    is.finite(samples$F)
]

reported_male_f <- samples$F[
  samples$REPORTED_SEX == 1L &
    is.finite(samples$F)
]

pdf(
  args$plot,
  width = 9,
  height = 7
)

if (length(finite_f)) {
  # Page 1: overall F distribution.
  hist(
    finite_f,
    breaks = 80,
    col = "grey80",
    border = "white",
    main = "chrX F distribution",
    xlab = "chrX inbreeding coefficient (F)"
  )

  abline(
    v = c(female_max_f, male_min_f),
    col = c("red", "blue"),
    lty = 2,
    lwd = 2
  )

  legend(
    "topright",
    legend = c(
      sprintf(
        "Female maximum F = %.4g",
        female_max_f
      ),
      sprintf(
        "Male minimum F = %.4g",
        male_min_f
      )
    ),
    col = c("red", "blue"),
    lty = 2,
    lwd = 2,
    bty = "n"
  )

  # Page 2: F distribution by reported sex.
  plot_groups <- list()

  if (length(reported_female_f)) {
    plot_groups[["Reported female"]] <-
      reported_female_f
  }

  if (length(reported_male_f)) {
    plot_groups[["Reported male"]] <-
      reported_male_f
  }

  if (length(plot_groups)) {
    boxplot(
      plot_groups,
      col = c(
        rgb(1, 0, 0, 0.35),
        rgb(0, 0, 1, 0.35)
      )[seq_along(plot_groups)],
      main = "chrX F by reported sex",
      ylab = "chrX inbreeding coefficient (F)"
    )

    abline(
      h = c(female_max_f, male_min_f),
      col = c("red", "blue"),
      lty = 2,
      lwd = 2
    )
  } else {
    plot.new()
    title("chrX F by reported sex")
    text(
      0.5,
      0.5,
      "No samples with known reported sex"
    )
  }
} else {
  plot.new()
  title("chrX F distribution")
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
    "female_max_f",
    "male_min_f",
    "n_samples",
    "n_samples_with_f",
    "n_genetic_female",
    "n_genetic_male",
    "n_reported_genetic_match",
    "n_discordant",
    "n_ambiguous",
    "n_missing_genetic_sex",
    "n_reported_sex_unknown",
    "n_excluded"
  ),
  value = c(
    female_max_f,
    male_min_f,
    nrow(samples),
    sum(has_f_value),
    sum(samples$GENETIC_SEX == 2L),
    sum(samples$GENETIC_SEX == 1L),
    sum(samples$REPORTED_GENETIC_SEX_MATCH),
    sum(samples$STATUS == "DISCORDANT"),
    sum(
      samples$STATUS ==
        "AMBIGUOUS_GENETIC_SEX"
    ),
    sum(
      samples$STATUS ==
        "MISSING_GENETIC_SEX"
    ),
    sum(
      samples$STATUS ==
        "REPORTED_SEX_UNKNOWN"
    ),
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