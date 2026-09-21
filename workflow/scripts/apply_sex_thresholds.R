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
    "--unmatched-genotypes",
    type = "character",
    help = "Genotype samples not matched to phenotype metadata"
  ),
  make_option(
    "--thresholds",
    type = "character",
    help = "Accepted threshold YAML file"
  ),
  make_option(
    "--exclude-ambiguous",
    type = "character",
    default = "true",
    help = paste(
      "Exclude ambiguous genetic-sex calls",
      "[default: %default]"
    )
  ),
  make_option(
    "--exclusions",
    type = "character",
    help = "Output PLINK-compatible FID/IID exclusion list"
  ),
  make_option(
    "--classification",
    type = "character",
    help = "Output sample-level classification TSV"
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
  "unmatched-genotypes",
  "thresholds",
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

# Parse Boolean option ----------------------------------------------------------

parse_boolean <- function(value, option_name) {
  normalized <- tolower(
    trimws(value)
  )

  if (normalized %in% c(
    "true",
    "t",
    "1",
    "yes",
    "y"
  )) {
    return(TRUE)
  }

  if (normalized %in% c(
    "false",
    "f",
    "0",
    "no",
    "n"
  )) {
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

# Read manually accepted thresholds --------------------------------------------

parse_thresholds <- function(path) {
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

  lines <- sub("#.*$", "", lines)
  lines <- trimws(lines)

  lines <- lines[
    nzchar(lines) &
      grepl(":", lines, fixed = TRUE)
  ]

  keys <- trimws(
    sub(":.*$", "", lines)
  )

  values <- trimws(
    sub("^[^:]*:", "", lines)
  )

  if (anyDuplicated(keys)) {
    duplicated_keys <- unique(
      keys[duplicated(keys)]
    )

    stop(
      "Threshold file contains duplicate key(s): ",
      paste(duplicated_keys, collapse = ", "),
      call. = FALSE
    )
  }

  thresholds <- setNames(
    suppressWarnings(
      as.numeric(values)
    ),
    keys
  )

  required_keys <- c(
    "female_max_f",
    "male_min_f"
  )

  missing_keys <- setdiff(
    required_keys,
    names(thresholds)
  )

  if (length(missing_keys)) {
    stop(
      "Threshold file must define: ",
      paste(missing_keys, collapse = ", "),
      call. = FALSE
    )
  }

  female_max_f <- unname(
    thresholds["female_max_f"]
  )

  male_min_f <- unname(
    thresholds["male_min_f"]
  )

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
      paste(
        "female_max_f must be lower",
        "than male_min_f."
      ),
      call. = FALSE
    )
  }

  c(
    female_max_f = female_max_f,
    male_min_f = male_min_f
  )
}

thresholds <- parse_thresholds(
  args$thresholds
)

female_max_f <- unname(
  thresholds["female_max_f"]
)

male_min_f <- unname(
  thresholds["male_min_f"]
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

# Read genotype samples not matched to phenotype metadata ----------------------

unmatched_genotypes <- fread(
  args[["unmatched-genotypes"]],
  data.table = FALSE,
  colClasses = "character"
)

if (
  !all(c("FID", "IID") %in% names(unmatched_genotypes))
) {
  stop(
    paste(
      "The unmatched-genotype file must contain",
      "FID and IID columns."
    ),
    call. = FALSE
  )
}

if (
  anyDuplicated(
    unmatched_genotypes[c("FID", "IID")]
  )
) {
  stop(
    paste(
      "The unmatched-genotype file contains",
      "duplicate FID/IID pairs."
    ),
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

classification$REPORTED_SEX <- suppressWarnings(
  as.integer(classification$SEX)
)

classification$REPORTED_SEX[
  is.na(classification$REPORTED_SEX) |
    !classification$REPORTED_SEX %in% c(1L, 2L)
] <- 0L

# Determine phenotype-metadata matching status ---------------------------------

classification_key <- paste(
  classification$FID,
  classification$IID,
  sep = "\r"
)

unmatched_key <- paste(
  unmatched_genotypes$FID,
  unmatched_genotypes$IID,
  sep = "\r"
)

classification$METADATA_MATCHED <-
  !classification_key %in% unmatched_key

classification$REPORTED_SEX_STATUS <- "AVAILABLE"

classification$REPORTED_SEX_STATUS[
  classification$REPORTED_SEX == 0L &
    !classification$METADATA_MATCHED
] <- "GENOTYPE_NOT_MATCHED_TO_PHENOTYPE"

classification$REPORTED_SEX_STATUS[
  classification$REPORTED_SEX == 0L &
    classification$METADATA_MATCHED
] <- "MISSING_OR_UNRECOGNISED_PHENOTYPE_SEX"

# Classify genetic sex ----------------------------------------------------------

classification$GENETIC_SEX <- 0L

classification$GENETIC_SEX[
  is.finite(classification$F) &
    classification$F <= female_max_f
] <- 2L

classification$GENETIC_SEX[
  is.finite(classification$F) &
    classification$F >= male_min_f
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

has_f_value <- is.finite(classification$F)

has_reported_sex <-
  classification$REPORTED_SEX %in% c(1L, 2L)

has_genetic_sex <-
  classification$GENETIC_SEX %in% c(1L, 2L)

classification$REPORTED_GENETIC_SEX_MATCH <-
  has_reported_sex &
    has_genetic_sex &
    classification$REPORTED_SEX ==
      classification$GENETIC_SEX

# Assign QC status --------------------------------------------------------------

classification$STATUS <- "OK"

classification$STATUS[
  !has_f_value
] <- "MISSING_GENETIC_SEX"

classification$STATUS[
  has_f_value &
    !has_genetic_sex
] <- "AMBIGUOUS_GENETIC_SEX"

classification$STATUS[
  has_reported_sex &
    has_genetic_sex &
    classification$REPORTED_SEX !=
      classification$GENETIC_SEX
] <- "DISCORDANT"

classification$STATUS[
  classification$REPORTED_SEX == 0L
] <- "REPORTED_SEX_UNKNOWN"

# Determine exclusions ----------------------------------------------------------

classification$EXCLUDE <-
  classification$STATUS == "DISCORDANT" |
    (
      exclude_ambiguous &
        classification$STATUS ==
          "AMBIGUOUS_GENETIC_SEX"
    )

# Validate reported-sex counts --------------------------------------------------

n_genotype_not_matched <- sum(
  classification$REPORTED_SEX_STATUS ==
    "GENOTYPE_NOT_MATCHED_TO_PHENOTYPE"
)

n_missing_or_unrecognised_sex <- sum(
  classification$REPORTED_SEX_STATUS ==
    "MISSING_OR_UNRECOGNISED_PHENOTYPE_SEX"
)

n_reported_sex_unknown <- sum(
  classification$REPORTED_SEX == 0L
)

if (
  n_genotype_not_matched +
    n_missing_or_unrecognised_sex !=
    n_reported_sex_unknown
) {
  stop(
    paste(
      "Internal consistency error:",
      "reported-sex unknown categories do not sum",
      "to the total unknown count."
    ),
    call. = FALSE
  )
}

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
  "METADATA_MATCHED",
  "REPORTED_SEX",
  "REPORTED_SEX_LABEL",
  "REPORTED_SEX_STATUS",
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
    "female_max_f",
    "male_min_f",
    "exclude_ambiguous",
    "n_samples",
    "n_samples_with_f",
    "n_genetic_female",
    "n_genetic_male",
    "n_reported_genetic_match",
    "n_discordant",
    "n_ambiguous",
    "n_missing_genetic_sex",
    "n_genotype_not_matched_to_phenotype",
    "n_missing_or_unrecognised_phenotype_sex",
    "n_reported_sex_unknown",
    "n_excluded"
  ),
  value = c(
    female_max_f,
    male_min_f,
    exclude_ambiguous,
    nrow(classification),
    sum(has_f_value),
    sum(classification$GENETIC_SEX == 2L),
    sum(classification$GENETIC_SEX == 1L),
    sum(classification$REPORTED_GENETIC_SEX_MATCH),
    sum(classification$STATUS == "DISCORDANT"),
    sum(
      classification$STATUS ==
        "AMBIGUOUS_GENETIC_SEX"
    ),
    sum(!has_f_value),
    n_genotype_not_matched,
    n_missing_or_unrecognised_sex,
    n_reported_sex_unknown,
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