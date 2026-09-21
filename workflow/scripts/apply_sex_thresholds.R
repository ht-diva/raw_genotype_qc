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
    help = "PLINK sex-check results"
  ),
  make_option(
    "--fam",
    type = "character",
    help = "PLINK FAM file"
  ),
  make_option(
    "--threshold-file",
    type = "character",
    help = "YAML file containing the accepted sex_f_threshold"
  ),
  make_option(
    "--exclusions",
    type = "character",
    help = "Output PLINK-compatible FID/IID exclusion list"
  ),
  make_option(
    "--classification",
    type = "character",
    help = "Output sample-level classification table"
  ),
  make_option(
    "--summary",
    type = "character",
    help = "Output summary table"
  )
)

args <- parse_args(
  OptionParser(option_list = option_list)
)

# Validate command-line arguments ----------------------------------------------

required_args <- c(
  "sexcheck",
  "fam",
  "threshold-file",
  "exclusions",
  "classification",
  "summary"
)

missing_args <- required_args[
  vapply(
    required_args,
    function(argument) {
      is.null(args[[argument]]) ||
        !nzchar(args[[argument]])
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

# Read the manually accepted threshold -----------------------------------------

parse_threshold <- function(path) {
  if (!file.exists(path)) {
    stop(
      "Accepted threshold file does not exist: ",
      path,
      call. = FALSE
    )
  }

  lines <- readLines(
    path,
    warn = FALSE
  )

  # Remove comments and empty lines.
  lines <- sub("#.*$", "", lines)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]

  threshold_line <- lines[
    grepl(
      "^sex_f_threshold[[:space:]]*:",
      lines
    )
  ]

  if (length(threshold_line) != 1L) {
    stop(
      paste(
        "Accepted threshold file must contain exactly one",
        "'sex_f_threshold:' entry."
      ),
      call. = FALSE
    )
  }

  threshold <- suppressWarnings(
    as.numeric(
      trimws(
        sub(
          "^[^:]+:",
          "",
          threshold_line
        )
      )
    )
  )

  if (!is.finite(threshold)) {
    stop(
      "sex_f_threshold must be a finite number.",
      call. = FALSE
    )
  }

  threshold
}

threshold <- parse_threshold(
  args[["threshold-file"]]
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

classification <- merge(
  fam[, c("FID", "IID", "SEX", ".FAM_ORDER")],
  sexcheck,
  by = join_columns,
  all.x = TRUE,
  sort = FALSE
)

classification <- classification[
  order(classification$.FAM_ORDER),
  ,
  drop = FALSE
]

classification$.FAM_ORDER <- NULL

if (nrow(classification) != nrow(fam)) {
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

classification$REPORTED_SEX <- suppressWarnings(
  as.integer(classification$SEX)
)

classification$REPORTED_SEX[
  is.na(classification$REPORTED_SEX) |
    !classification$REPORTED_SEX %in% c(1L, 2L)
] <- 0L

# Classify genetic sex ----------------------------------------------------------

# Classification rule:
#   F <= threshold = female
#   F > threshold  = male
#   missing F      = unknown

classification$GENETIC_SEX <- 0L

classification$GENETIC_SEX[
  is.finite(classification$F) &
    classification$F <= threshold
] <- 2L

classification$GENETIC_SEX[
  is.finite(classification$F) &
    classification$F > threshold
] <- 1L

# Add readable labels -----------------------------------------------------------

sex_labels <- c(
  `0` = "unknown",
  `1` = "male",
  `2` = "female"
)

classification$REPORTED_SEX_LABEL <- unname(
  sex_labels[
    as.character(classification$REPORTED_SEX)
  ]
)

classification$GENETIC_SEX_LABEL <- unname(
  sex_labels[
    as.character(classification$GENETIC_SEX)
  ]
)

# Compare reported and genetic sex ---------------------------------------------

has_reported_sex <-
  classification$REPORTED_SEX %in% c(1L, 2L)

has_genetic_sex <-
  classification$GENETIC_SEX %in% c(1L, 2L)

has_f_value <- is.finite(classification$F)

classification$REPORTED_GENETIC_SEX_MATCH <-
  has_reported_sex &
  has_genetic_sex &
  classification$REPORTED_SEX ==
    classification$GENETIC_SEX

# Assign sample status ----------------------------------------------------------

classification$STATUS <- "OK"

classification$STATUS[
  !has_f_value
] <- "MISSING_GENETIC_SEX"

classification$STATUS[
  has_reported_sex &
    has_genetic_sex &
    classification$REPORTED_SEX !=
      classification$GENETIC_SEX
] <- "DISCORDANT"

# Unknown reported sex takes precedence.
classification$STATUS[
  classification$REPORTED_SEX == 0L
] <- "REPORTED_SEX_UNKNOWN"

# With one threshold there is no ambiguous interval.
# Only discordant samples are excluded.
classification$EXCLUDE <-
  classification$STATUS == "DISCORDANT"

# Create output directories -----------------------------------------------------

output_paths <- c(
  args$exclusions,
  args$classification,
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
  classification[, output_columns, drop = FALSE],
  file = args$classification,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  na = "NA"
)

# Write PLINK-compatible exclusion list ----------------------------------------

exclusions <- classification[
  classification$EXCLUDE,
  c("FID", "IID"),
  drop = FALSE
]

names(exclusions)[1] <- "#FID"

write.table(
  exclusions,
  file = args$exclusions,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

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
    nrow(classification),
    sum(has_f_value),
    sum(classification$GENETIC_SEX == 2L),
    sum(classification$GENETIC_SEX == 1L),
    sum(classification$REPORTED_GENETIC_SEX_MATCH),
    sum(classification$STATUS == "DISCORDANT"),
    sum(
      classification$STATUS ==
        "MISSING_GENETIC_SEX"
    ),
    sum(
      classification$STATUS ==
        "REPORTED_SEX_UNKNOWN"
    ),
    sum(classification$EXCLUDE)
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