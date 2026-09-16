#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--sexcheck", type = "character",
              help = "PLINK sex-check results"),
  make_option("--fam", type = "character",
              help = "PLINK FAM file"),
  make_option("--thresholds", type = "character",
              help = "Accepted threshold file"),
  make_option("--exclude-ambiguous", type = "character", default = "true",
              help = "Exclude ambiguous genetic-sex calls [default: %default]"),
  make_option("--exclusions", type = "character",
              help = "Output FID/IID exclusion list"),
  make_option("--classification", type = "character",
              help = "Output sample-level classification table"),
  make_option("--summary", type = "character",
              help = "Output summary table")
)

args <- parse_args(OptionParser(option_list = option_list))

# Validate command-line arguments ----------------------------------------------

required_args <- c(
  "sexcheck",
  "fam",
  "thresholds",
  "exclusions",
  "classification",
  "summary"
)

missing_args <- required_args[
  vapply(required_args, function(x) {
    is.null(args[[x]]) || !nzchar(args[[x]])
  }, logical(1))
]

if (length(missing_args)) {
  stop(
    "Missing required option(s): ",
    paste0("--", missing_args, collapse = ", "),
    call. = FALSE
  )
}

parse_boolean <- function(value, option_name) {
  normalized <- tolower(trimws(value))

  if (normalized %in% c("true", "t", "1", "yes", "y")) {
    return(TRUE)
  }

  if (normalized %in% c("false", "f", "0", "no", "n")) {
    return(FALSE)
  }

  stop(
    option_name,
    " must be true/false, yes/no, or 1/0; received: ",
    value,
    call. = FALSE
  )
}

exclude_ambiguous <- parse_boolean(
  args[["exclude-ambiguous"]],
  "--exclude-ambiguous"
)

# Read and validate accepted thresholds ----------------------------------------

parse_thresholds <- function(path) {
  if (!file.exists(path)) {
    stop("Threshold file does not exist: ", path, call. = FALSE)
  }

  lines <- readLines(path, warn = FALSE)
  lines <- sub("#.*$", "", lines)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines) & grepl(":", lines, fixed = TRUE)]

  keys <- trimws(sub(":.*$", "", lines))
  values <- trimws(sub("^[^:]*:", "", lines))

  if (anyDuplicated(keys)) {
    duplicated_keys <- unique(keys[duplicated(keys)])
    stop(
      "Threshold file contains duplicate key(s): ",
      paste(duplicated_keys, collapse = ", "),
      call. = FALSE
    )
  }

  thresholds <- setNames(
    suppressWarnings(as.numeric(values)),
    keys
  )

  required_keys <- c("female_max_f", "male_min_f")
  missing_keys <- setdiff(required_keys, names(thresholds))

  if (length(missing_keys)) {
    stop(
      "Threshold file must define: ",
      paste(missing_keys, collapse = ", "),
      call. = FALSE
    )
  }

  female_max_f <- unname(thresholds["female_max_f"])
  male_min_f <- unname(thresholds["male_min_f"])

  if (!is.finite(female_max_f) || !is.finite(male_min_f)) {
    stop("Sex-classification thresholds must be finite numbers.", call. = FALSE)
  }

  if (female_max_f >= male_min_f) {
    stop(
      "female_max_f must be lower than male_min_f.",
      call. = FALSE
    )
  }

  c(
    female_max_f = female_max_f,
    male_min_f = male_min_f
  )
}

thresholds <- parse_thresholds(args$thresholds)
female_max_f <- unname(thresholds["female_max_f"])
male_min_f <- unname(thresholds["male_min_f"])

# Read sex-check results --------------------------------------------------------

sexcheck <- fread(args$sexcheck, data.table = FALSE)

if (!"IID" %in% names(sexcheck)) {
  stop("Sex-check file is missing the IID column.", call. = FALSE)
}

f_column <- intersect(c("F", "XF", "X_F"), names(sexcheck))

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
sexcheck$F <- suppressWarnings(as.numeric(sexcheck$F))

if (anyDuplicated(sexcheck[sexcheck_id_columns])) {
  stop(
    "Sex-check file contains duplicate sample identifiers.",
    call. = FALSE
  )
}

# Read the FAM file -------------------------------------------------------------

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
names(fam) <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENO")

if (anyDuplicated(fam[c("FID", "IID")])) {
  stop("FAM file contains duplicate FID/IID pairs.", call. = FALSE)
}

# Join records while preserving FAM order --------------------------------------

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
    "Internal join error: classification row count differs from FAM row count.",
    call. = FALSE
  )
}

# Classify reported and genetic sex --------------------------------------------

classification$REPORTED_SEX <- suppressWarnings(
  as.integer(classification$SEX)
)

classification$REPORTED_SEX[
  is.na(classification$REPORTED_SEX) |
    !classification$REPORTED_SEX %in% c(1L, 2L)
] <- 0L

classification$GENETIC_SEX <- 0L

classification$GENETIC_SEX[
  is.finite(classification$F) &
    classification$F <= female_max_f
] <- 2L

classification$GENETIC_SEX[
  is.finite(classification$F) &
    classification$F >= male_min_f
] <- 1L

has_reported_sex <- classification$REPORTED_SEX %in% c(1L, 2L)
has_genetic_sex <- classification$GENETIC_SEX %in% c(1L, 2L)
has_f_value <- is.finite(classification$F)

classification$STATUS <- "OK"

classification$STATUS[!has_f_value] <- "MISSING_GENETIC_SEX"

classification$STATUS[
  has_reported_sex &
    has_f_value &
    !has_genetic_sex
] <- "AMBIGUOUS_GENETIC_SEX"

classification$STATUS[
  has_reported_sex &
    has_genetic_sex &
    classification$REPORTED_SEX != classification$GENETIC_SEX
] <- "DISCORDANT"

# Give unknown reported sex precedence over genetic-sex classifications.
classification$STATUS[
  classification$REPORTED_SEX == 0L
] <- "REPORTED_SEX_UNKNOWN"

classification$EXCLUDE <-
  classification$STATUS == "DISCORDANT" |
  (
    exclude_ambiguous &
      classification$STATUS == "AMBIGUOUS_GENETIC_SEX"
  )

# Write sample-level classification --------------------------------------------

write.table(
  classification,
  file = args$classification,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
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

# Write summary ----------------------------------------------------------------

summary_table <- data.frame(
  metric = c(
    "female_max_f",
    "male_min_f",
    "exclude_ambiguous",
    "n_samples",
    "n_samples_with_f",
    "n_missing_genetic_sex",
    "n_discordant",
    "n_ambiguous",
    "n_reported_unknown",
    "n_excluded"
  ),
  value = c(
    female_max_f,
    male_min_f,
    exclude_ambiguous,
    nrow(classification),
    sum(is.finite(classification$F)),
    sum(classification$STATUS == "MISSING_GENETIC_SEX"),
    sum(classification$STATUS == "DISCORDANT"),
    sum(classification$STATUS == "AMBIGUOUS_GENETIC_SEX"),
    sum(classification$STATUS == "REPORTED_SEX_UNKNOWN"),
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