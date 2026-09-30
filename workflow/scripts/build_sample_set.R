#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

# -------------------------------------------------------------------------
# Command-line arguments
# -------------------------------------------------------------------------

option_list <- list(
  make_option(
    "--fam",
    type = "character",
    help = "Input PLINK FAM file"
  ),
  make_option(
    "--exclude-files",
    type = "character",
    help = "Semicolon-separated list of reviewed sample-exclusion files"
  ),
  make_option(
    "--exclusions-output",
    type = "character",
    help = "Output table containing excluded samples"
  ),
  make_option(
    "--keep-output",
    type = "character",
    help = "Output PLINK keep file containing retained FID/IID pairs"
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

# -------------------------------------------------------------------------
# Read FAM file
# -------------------------------------------------------------------------

fam <- fread(
  args$fam,
  header = FALSE,
  data.table = FALSE,
  colClasses = "character"
)

if (ncol(fam) < 6) {
  stop("FAM file must contain at least 6 columns.")
}

fam <- fam[, 1:6]

colnames(fam) <- c(
  "FID",
  "IID",
  "PAT",
  "MAT",
  "SEX",
  "PHENO"
)

# -------------------------------------------------------------------------
# Parse exclusion-file list
# -------------------------------------------------------------------------

exclude_files <- strsplit(
  args$`exclude-files`,
  ";",
  fixed = TRUE
)[[1]]

exclude_files <- trimws(exclude_files)
exclude_files <- exclude_files[nzchar(exclude_files)]

# -------------------------------------------------------------------------
# Read one exclusion file
# -------------------------------------------------------------------------

read_exclusions <- function(path, valid_iids) {

  if (!file.exists(path)) {
    stop(
      "Required reviewed exclusion file does not exist: ",
      path
    )
  }

  lines <- readLines(
    path,
    warn = FALSE
  )

  lines <- trimws(lines)

  # Remove empty lines and comments.
  lines <- lines[
    nzchar(lines) &
      !grepl("^#", lines)
  ]

  if (length(lines) == 0L) {
    return(character())
  }

  fields <- strsplit(
    lines,
    "[[:space:]]+"
  )

  excluded_iids <- character()

  for (row in fields) {

    # Skip possible header lines.
    first_field <- toupper(
      sub("^#", "", row[1])
    )

    if (first_field %in% c("FID", "IID")) {
      next
    }

    # Prefer the second column when the file contains FID IID.
    if (
      length(row) >= 2L &&
      row[2] %in% valid_iids
    ) {
      excluded_iids <- c(
        excluded_iids,
        row[2]
      )

    # Otherwise allow files containing IID only.
    } else if (row[1] %in% valid_iids) {
      excluded_iids <- c(
        excluded_iids,
        row[1]
      )
    }
  }

  unique(excluded_iids)
}

# -------------------------------------------------------------------------
# Combine reviewed QC exclusions
# -------------------------------------------------------------------------

requested_exclusions <- unique(
  unlist(
    lapply(
      exclude_files,
      read_exclusions,
      valid_iids = fam$IID
    )
  )
)

unknown_exclusions <- setdiff(
  requested_exclusions,
  fam$IID
)

excluded_iids <- intersect(
  requested_exclusions,
  fam$IID
)

# -------------------------------------------------------------------------
# Build keep and exclusion tables
# -------------------------------------------------------------------------

keep_samples <- fam[
  !fam$IID %in% excluded_iids,
  c("FID", "IID")
]

excluded_samples <- fam[
  fam$IID %in% excluded_iids,
  c("FID", "IID")
]

excluded_samples$REASON <- "accepted_qc_exclusion"

# -------------------------------------------------------------------------
# Write outputs
# -------------------------------------------------------------------------

# PLINK --keep compatible file: no header.
write.table(
  keep_samples,
  file = args$`keep-output`,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

write.table(
  excluded_samples,
  file = args$`exclusions-output`,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# -------------------------------------------------------------------------
# Summary
# -------------------------------------------------------------------------

summary_table <- data.frame(
  metric = c(
    "n_input_samples",
    "n_excluded",
    "n_kept",
    "n_requested_exclusions_not_in_fam"
  ),
  value = c(
    nrow(fam),
    length(excluded_iids),
    nrow(keep_samples),
    length(unknown_exclusions)
  )
)

write.table(
  summary_table,
  file = args$summary,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)