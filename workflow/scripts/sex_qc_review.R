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
    "--unmatched-genotypes",
    type = "character",
    help = "Genotype samples not matched to phenotype metadata"
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
    "--exclude-ambiguous",
    type = "character",
    help = "Exclude ambiguous genetic-sex calls"
  ),
  make_option(
    "--exclude-discordant",
    type = "character",
    help = "Exclude discordant reported/genetic-sex calls"
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
  "unmatched-genotypes",
  "female-max-f",
  "male-min-f",
  "exclude-ambiguous",
  "exclude-discordant",
  "plot",
  "table",
  "candidate-thresholds",
  "summary"
)

missing_args <- required_args[
  vapply(
    required_args,
    function(argument) {
      is.null(args[[argument]]) ||
        length(args[[argument]]) == 0L
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

parse_boolean <- function(value, option_name) {
  normalized <- tolower(trimws(value))
  if (normalized %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (normalized %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop(option_name, " must be true/false, yes/no, or 1/0; received: ", value,
       call. = FALSE)
}

exclude_ambiguous <- parse_boolean(
  args[["exclude-ambiguous"]], "--exclude-ambiguous"
)
exclude_discordant <- parse_boolean(
  args[["exclude-discordant"]], "--exclude-discordant"
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

# Determine phenotype-metadata matching status ---------------------------------

sample_key <- paste(
  samples$FID,
  samples$IID,
  sep = "\r"
)

unmatched_key <- paste(
  unmatched_genotypes$FID,
  unmatched_genotypes$IID,
  sep = "\r"
)

samples$METADATA_MATCHED <-
  !sample_key %in% unmatched_key

samples$REPORTED_SEX_STATUS <- "AVAILABLE"

samples$REPORTED_SEX_STATUS[
  samples$REPORTED_SEX == 0L &
    !samples$METADATA_MATCHED
] <- "GENOTYPE_NOT_MATCHED_TO_PHENOTYPE"

samples$REPORTED_SEX_STATUS[
  samples$REPORTED_SEX == 0L &
    samples$METADATA_MATCHED
] <- "MISSING_OR_UNRECOGNISED_PHENOTYPE_SEX"

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

# Assign QC status --------------------------------------------------------------

samples$STATUS <- "OK"
samples$REVIEW_FLAG <- ""

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

samples$STATUS[
  has_genetic_sex &
    samples$REPORTED_SEX_STATUS ==
      "GENOTYPE_NOT_MATCHED_TO_PHENOTYPE"
] <- "GENOTYPE_NOT_MATCHED_TO_PHENOTYPE"

samples$REVIEW_FLAG[
  samples$REPORTED_SEX_STATUS ==
    "GENOTYPE_NOT_MATCHED_TO_PHENOTYPE"
] <- "WARNING / investigate"

samples$STATUS[
  has_genetic_sex &
    samples$REPORTED_SEX_STATUS ==
      "MISSING_OR_UNRECOGNISED_PHENOTYPE_SEX"
] <- "MISSING_OR_UNRECOGNISED_PHENOTYPE_SEX"

# Candidate exclusions:
#   - discordant reported and genetic sex
#   - ambiguous genetic sex
samples$EXCLUDE <-
  (exclude_discordant & samples$STATUS == "DISCORDANT") |
  (exclude_ambiguous & samples$STATUS == "AMBIGUOUS_GENETIC_SEX")

# Validate reported-sex counts --------------------------------------------------

n_genotype_not_matched <- sum(
  samples$REPORTED_SEX_STATUS ==
    "GENOTYPE_NOT_MATCHED_TO_PHENOTYPE"
)

n_missing_or_unrecognised_sex <- sum(
  samples$REPORTED_SEX_STATUS ==
    "MISSING_OR_UNRECOGNISED_PHENOTYPE_SEX"
)

n_reported_sex_unknown <- sum(
  samples$REPORTED_SEX == 0L
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
  "REVIEW_FLAG",
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

  # Page 2: F values grouped by reported sex.
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
    "exclude_ambiguous",
    "exclude_discordant",
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
    exclude_discordant,
    nrow(samples),
    sum(has_f_value),
    sum(samples$GENETIC_SEX == 2L),
    sum(samples$GENETIC_SEX == 1L),
    sum(samples$REPORTED_GENETIC_SEX_MATCH),
    sum(samples$STATUS == "DISCORDANT"),
    sum(samples$STATUS == "AMBIGUOUS_GENETIC_SEX"),
    sum(!has_f_value),
    n_genotype_not_matched,
    n_missing_or_unrecognised_sex,
    n_reported_sex_unknown,
    sum(samples$EXCLUDE)
  ),
  stringsAsFactors = FALSE
)

summary_table$warning <- ""
summary_table$warning[
  summary_table$metric ==
    "n_genotype_not_matched_to_phenotype"
] <- "WARNING / investigate; not automatically excluded"

write.table(
  summary_table,
  file = args$summary,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
