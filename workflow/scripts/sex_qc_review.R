#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

options <- list(
  make_option("--sexcheck", type = "character",
              help = "PLINK sex-check file"),
  make_option("--fam", type = "character",
              help = "PLINK FAM file"),
  make_option("--candidate-female-max-f", type = "double",
              help = "Maximum chrX F for candidate females"),
  make_option("--candidate-male-min-f", type = "double",
              help = "Minimum chrX F for candidate males"),
  make_option("--tail-quantile", type = "double", default = 0.005,
              help = "Tail quantile used for empirical thresholds [default: %default]"),
  make_option("--plot", type = "character",
              help = "Output diagnostic PDF"),
  make_option("--table", type = "character",
              help = "Output sample-level TSV"),
  make_option("--suggested-thresholds", type = "character",
              help = "Output suggested-threshold YAML"),
  make_option("--summary", type = "character",
              help = "Output summary TSV")
)

args <- parse_args(OptionParser(option_list = options))

required_args <- c(
  "sexcheck",
  "fam",
  "candidate-female-max-f",
  "candidate-male-min-f",
  "plot",
  "table",
  "suggested-thresholds",
  "summary"
)

missing_args <- required_args[
  vapply(required_args, function(x) is.null(args[[x]]), logical(1))
]

if (length(missing_args)) {
  stop(
    "Missing required option(s): ",
    paste0("--", missing_args, collapse = ", "),
    call. = FALSE
  )
}

female_max_f <- args[["candidate-female-max-f"]]
male_min_f   <- args[["candidate-male-min-f"]]

if (!is.finite(female_max_f) || !is.finite(male_min_f)) {
  stop("Candidate F thresholds must be finite numbers.", call. = FALSE)
}

if (female_max_f >= male_min_f) {
  stop(
    "--candidate-female-max-f must be less than --candidate-male-min-f.",
    call. = FALSE
  )
}

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

id_columns <- if ("FID" %in% names(sexcheck)) c("FID", "IID") else "IID"
sexcheck <- sexcheck[, c(id_columns, f_column[1]), drop = FALSE]
names(sexcheck)[ncol(sexcheck)] <- "F"
sexcheck$F <- suppressWarnings(as.numeric(sexcheck$F))

if (anyDuplicated(sexcheck[id_columns])) {
  stop(
    "Sex-check file contains duplicate sample identifiers.",
    call. = FALSE
  )
}

# Read FAM ---------------------------------------------------------------------

fam <- fread(
  args$fam,
  header = FALSE,
  data.table = FALSE,
  colClasses = "character"
)

if (ncol(fam) < 6L) {
  stop("Malformed FAM file: expected at least six columns.", call. = FALSE)
}

fam <- fam[, seq_len(6), drop = FALSE]
names(fam) <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENO")

if (anyDuplicated(fam[c("FID", "IID")])) {
  stop("FAM file contains duplicate FID/IID pairs.", call. = FALSE)
}

# Join while preserving FAM order ----------------------------------------------

fam$.FAM_ORDER <- seq_len(nrow(fam))
join_columns <- if ("FID" %in% names(sexcheck)) c("FID", "IID") else "IID"

samples <- merge(
  fam[, c("FID", "IID", "SEX", ".FAM_ORDER")],
  sexcheck,
  by = join_columns,
  all.x = TRUE,
  sort = FALSE
)

samples <- samples[order(samples$.FAM_ORDER), ]
samples$.FAM_ORDER <- NULL

# Classify candidate genetic sex -----------------------------------------------

samples$REPORTED_SEX <- suppressWarnings(as.integer(samples$SEX))
samples$REPORTED_SEX[
  is.na(samples$REPORTED_SEX) |
    !samples$REPORTED_SEX %in% c(1L, 2L)
] <- 0L

samples$CANDIDATE_GENETIC_SEX <- 0L
samples$CANDIDATE_GENETIC_SEX[
  is.finite(samples$F) & samples$F <= female_max_f
] <- 2L
samples$CANDIDATE_GENETIC_SEX[
  is.finite(samples$F) & samples$F >= male_min_f
] <- 1L

samples$CANDIDATE_STATUS <- "OK"

samples$CANDIDATE_STATUS[!is.finite(samples$F)] <-
  "MISSING_GENETIC_SEX"

samples$CANDIDATE_STATUS[samples$REPORTED_SEX == 0L] <-
  "REPORTED_SEX_UNKNOWN"

known_reported_sex <- samples$REPORTED_SEX %in% c(1L, 2L)
known_genetic_sex  <- samples$CANDIDATE_GENETIC_SEX %in% c(1L, 2L)

samples$CANDIDATE_STATUS[
  known_reported_sex &
    is.finite(samples$F) &
    !known_genetic_sex
] <- "AMBIGUOUS_GENETIC_SEX"

samples$CANDIDATE_STATUS[
  known_reported_sex &
    known_genetic_sex &
    samples$REPORTED_SEX != samples$CANDIDATE_GENETIC_SEX
] <- "DISCORDANT"

write.table(
  samples,
  file = args$table,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Calculate empirical threshold suggestions -----------------------------------

tail_quantile <- max(0.0001, min(0.1, args[["tail-quantile"]]))

reported_female_f <- samples$F[
  samples$REPORTED_SEX == 2L & is.finite(samples$F)
]

reported_male_f <- samples$F[
  samples$REPORTED_SEX == 1L & is.finite(samples$F)
]

empirical_female_max <- if (length(reported_female_f) >= 20L) {
  as.numeric(quantile(
    reported_female_f,
    probs = 1 - tail_quantile,
    names = FALSE,
    na.rm = TRUE
  ))
} else {
  NA_real_
}

empirical_male_min <- if (length(reported_male_f) >= 20L) {
  as.numeric(quantile(
    reported_male_f,
    probs = tail_quantile,
    names = FALSE,
    na.rm = TRUE
  ))
} else {
  NA_real_
}

clusters_separated <-
  is.finite(empirical_female_max) &&
  is.finite(empirical_male_min) &&
  empirical_female_max < empirical_male_min

suggested_female_max <- if (clusters_separated) {
  empirical_female_max
} else {
  female_max_f
}

suggested_male_min <- if (clusters_separated) {
  empirical_male_min
} else {
  male_min_f
}

threshold_lines <- c(
  "# Automatic diagnostic suggestion only; user approval is required.",
  sprintf("female_max_f: %.8g", suggested_female_max),
  sprintf("male_min_f: %.8g", suggested_male_min),
  sprintf("# empirical_tail_quantile: %.8g", tail_quantile),
  sprintf(
    "# empirical_clusters_separated: %s",
    tolower(as.character(clusters_separated))
  )
)

writeLines(threshold_lines, args[["suggested-thresholds"]])

# Create diagnostic plots ------------------------------------------------------

finite_f <- samples$F[is.finite(samples$F)]

pdf(args$plot, width = 8, height = 6)
on.exit(dev.off(), add = TRUE)

if (length(finite_f)) {
  plot_range <- range(finite_f)

  hist(
    reported_female_f,
    breaks = 80,
    col = rgb(1, 0, 0, 0.35),
    main = "chrX F by reported sex",
    xlab = "chrX inbreeding coefficient (F)",
    xlim = plot_range
  )

  hist(
    reported_male_f,
    breaks = 80,
    col = rgb(0, 0, 1, 0.35),
    add = TRUE
  )

  abline(v = c(female_max_f, male_min_f), lty = 2)

  legend(
    "topright",
    legend = c(
      "Reported female",
      "Reported male",
      "Candidate cutoffs"
    ),
    fill = c(
      rgb(1, 0, 0, 0.35),
      rgb(0, 0, 1, 0.35),
      NA
    ),
    lty = c(NA, NA, 2),
    bty = "n"
  )

  known_sex_f <- samples$F[
    samples$REPORTED_SEX %in% c(1L, 2L) &
      is.finite(samples$F)
  ]

  hist(
    known_sex_f,
    breaks = 80,
    main = "Empirical diagnostic suggestion",
    xlab = "chrX inbreeding coefficient (F)"
  )

  abline(
    v = c(suggested_female_max, suggested_male_min),
    lty = 2
  )
} else {
  plot.new()
  title("chrX F by reported sex")
  text(0.5, 0.5, "No finite chrX F values available")
}

dev.off()
on.exit(NULL, add = FALSE)

# Write summary ---------------------------------------------------------------

summary_table <- data.frame(
  metric = c(
    "n_samples",
    "n_samples_with_f",
    "candidate_female_max_f",
    "candidate_male_min_f",
    "empirical_suggested_female_max_f",
    "empirical_suggested_male_min_f",
    "empirical_clusters_separated",
    "n_candidate_discordant",
    "n_candidate_ambiguous",
    "n_missing_genetic_sex",
    "n_reported_sex_unknown"
  ),
  value = c(
    nrow(samples),
    sum(is.finite(samples$F)),
    female_max_f,
    male_min_f,
    suggested_female_max,
    suggested_male_min,
    clusters_separated,
    sum(samples$CANDIDATE_STATUS == "DISCORDANT"),
    sum(samples$CANDIDATE_STATUS == "AMBIGUOUS_GENETIC_SEX"),
    sum(samples$CANDIDATE_STATUS == "MISSING_GENETIC_SEX"),
    sum(samples$REPORTED_SEX == 0L)
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