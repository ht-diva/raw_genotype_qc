#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(bigsnpr)
  library(bigstatsr)
})

option_list <- list(
  make_option(
    "--bed-prefix",
    type = "character",
    help = "Prefix of the PLINK BED/BIM/FAM files"
  ),
  make_option(
    "--n-pcs",
    type = "integer",
    default = 5,
    help = "Number of principal components [default: %default]"
  ),
  make_option(
    "--autosvd-max-iter",
    type = "integer",
    default = 1,
    help = "Maximum number of bed_autoSVD iterations [default: %default]"
  ),
  make_option(
    "--autosvd-min-maf",
    type = "double",
    default = 0.01,
    help = "Minimum MAF used by bed_autoSVD [default: %default]"
  ),
  make_option(
    "--keep-adjusted-p-above",
    type = "double",
    default = 0.05,
    help = paste(
      "Keep variants with adjusted pcadapt p-value above this threshold",
      "[default: %default]"
    )
  ),
  make_option(
    "--keep-variants",
    type = "character",
    help = "Output file containing variants retained for KING"
  ),
  make_option(
    "--excluded-variants",
    type = "character",
    help = "Output table containing variants excluded by pcadapt"
  ),
  make_option(
    "--plot-output",
    type = "character",
    help = "Output PDF containing pcadapt diagnostics"
  ),
  make_option(
    "--summary-output",
    type = "character",
    help = "Output summary TSV"
  ),
  make_option(
    "--threads",
    type = "integer",
    default = 1,
    help = "Number of threads [default: %default]"
  )
)

args <- parse_args(OptionParser(option_list = option_list))

# Load the PLINK BED dataset.
geno <- bed(paste0(args$`bed-prefix`, ".bed"))

n_samples <- nrow(geno)
n_variants <- ncol(geno)

n_pcs <- min(
  args$`n-pcs`,
  n_samples - 1L,
  n_variants - 1L
)

if (n_pcs < 1L) {
  stop("Too few samples or variants to run pcadapt.")
}

# Historical BELIEVE procedure:
# 1. Run a preliminary autoSVD.
# 2. Use those PCs as covariates in bed_pcadapt.
svd <- bed_autoSVD(
  geno,
  k = n_pcs,
  min.maf = args$`autosvd-min-maf`,
  max.iter = args$`autosvd-max-iter`,
  ncores = args$threads
)

pcadapt_fit <- bed_pcadapt(
  geno,
  svd$u[, seq_len(n_pcs), drop = FALSE],
  ncores = args$threads
)

adjusted_p <- as.numeric(
  predict(pcadapt_fit, log10 = FALSE)
)

if (length(adjusted_p) != n_variants) {
  stop(
    "Unexpected number of pcadapt adjusted p-values: ",
    length(adjusted_p),
    " values returned for ",
    n_variants,
    " variants."
  )
}

variant_ids <- as.character(geno$map$marker.ID)

keep_idx <- which(
  is.finite(adjusted_p) &
    adjusted_p > args$`keep-adjusted-p-above`
)

if (length(keep_idx) < 100L) {
  stop(
    "pcadapt filtering retained fewer than 100 variants. ",
    "Review the filtering threshold and input data before running KING."
  )
}

writeLines(
  variant_ids[keep_idx],
  args$`keep-variants`
)

variant_table <- data.frame(
  CHR = geno$map$chromosome,
  POS = geno$map$physical.pos,
  ID = variant_ids,
  PCADAPT_ADJUSTED_P = adjusted_p,
  EXCLUDED = !(seq_len(n_variants) %in% keep_idx)
)

write.table(
  variant_table[variant_table$EXCLUDED, ],
  file = args$`excluded-variants`,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

pdf(
  args$`plot-output`,
  width = 10,
  height = 7
)

plot(svd, type = "scores", scores = seq_len(n_pcs))
plot(svd, type = "loadings", loadings = seq_len(n_pcs))
plot(pcadapt_fit, type = "Manhattan")
plot(pcadapt_fit, type = "Q-Q")

dev.off()

summary_table <- data.frame(
  metric = c(
    "n_samples",
    "n_variants_input",
    "n_pcs_preliminary",
    "autosvd_max_iter",
    "autosvd_min_maf",
    "pcadapt_keep_adjusted_p_above",
    "n_variants_kept_for_KING",
    "n_variants_excluded_for_population_structure"
  ),
  value = c(
    n_samples,
    n_variants,
    n_pcs,
    args$`autosvd-max-iter`,
    args$`autosvd-min-maf`,
    args$`keep-adjusted-p-above`,
    length(keep_idx),
    n_variants - length(keep_idx)
  )
)

write.table(
  summary_table,
  file = args$`summary-output`,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
