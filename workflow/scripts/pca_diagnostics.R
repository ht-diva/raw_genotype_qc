#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file)), "pca_helpers.R"))

parsed <- parse_args(OptionParser(option_list = list(
  make_option("--method", type = "character", default = "autosvd"),
  make_option("--bed-prefix", type = "character"),
  make_option("--unrelated-file", type = "character"),
  make_option("--keep-file", type = "character"),
  make_option("--n-pcs", type = "integer", default = 20),
  make_option("--robust-distance-n-pcs", type = "integer", default = 6),
  make_option("--robust-mad-multiplier", type = "double", default = 6),
  make_option("--autosvd-max-iter", type = "integer", default = 5),
  make_option("--autosvd-min-maf", type = "double", default = 0.01),
  make_option("--seed", type = "integer", default = 1),
  make_option("--pcs-output", type = "character"),
  make_option("--scores-plot", type = "character"),
  make_option("--loadings-plot", type = "character"),
  make_option("--loadings-output", type = "character"),
  make_option("--lrldr-output", type = "character"),
  make_option("--robust-output", type = "character"),
  make_option("--suggested-outliers", type = "character"),
  make_option("--retained-variants", type = "character"),
  make_option("--summary-output", type = "character"),
  make_option("--threads", type = "integer", default = 1)
)), positional_arguments = TRUE)
o <- parsed$options
if (!o$method %in% c("plink2", "autosvd")) stop("Expected plink2 or autosvd for --method")
expected_args <- if (o$method == "plink2") 5L else 0L
if (length(parsed$args) != expected_args) stop("Unexpected PCA input arguments")
outputs <- c("pcs-output", "scores-plot", "loadings-plot", "loadings-output",
             "robust-output", "suggested-outliers", "retained-variants", "summary-output")
require_options(o, c("bed-prefix", "unrelated-file", "keep-file", outputs))
check_number(o$`n-pcs`, "n-pcs", 2, integer = TRUE)
check_number(o$`robust-distance-n-pcs`, "robust-distance-n-pcs", 2, integer = TRUE)
check_number(o$`robust-mad-multiplier`, "robust-mad-multiplier", 0)
check_number(o$threads, "threads", 1, integer = TRUE)
check_number(o$seed, "seed", 0, .Machine$integer.max, integer = TRUE)
set.seed(o$seed)
check_paths(paste0(o$`bed-prefix`, c(".bed", ".bim", ".fam")))
for (path in unlist(o[outputs])) dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

# Read IDs as character strings, preserving leading zeroes and FID/IID identity.
fam <- read_fam_ids(o$`bed-prefix`)
ind.keep <- read_id_indices(o$`keep-file`, fam)
ind.unrel <- read_id_indices(o$`unrelated-file`, fam)
if (any(!ind.unrel %in% ind.keep)) stop("Unrelated samples must be a subset of the retained samples")
if (length(ind.keep) < 10L || length(ind.unrel) < 10L) stop("Too few retained/unrelated samples for PCA (minimum 10)")
proj.rows <- setdiff(ind.keep, ind.unrel)
if (o$method == "plink2") {
  # Use the same score calculation for every retained sample. Scale each axis
  # to unit length in the training set and align its sign to the fitted PC.
  check_paths(parsed$args)
  read_plink_table <- function(path, id_columns) {
    header <- strsplit(readLines(path, n = 1L), "\t", fixed = TRUE)[[1]]
    header <- sub("^#", "", header)
    classes <- ifelse(header %in% id_columns, "character", "numeric")
    read.delim(path, col.names = header, colClasses = classes, check.names = FALSE,
               quote = "", na.strings = character())
  }
  fit <- read_plink_table(parsed$args[1], c("FID", "IID"))
  eigenval <- scan(parsed$args[2], quiet = TRUE)
  weights <- read_plink_table(parsed$args[3], c("CHROM", "ID", "REF", "ALT", "A1"))
  scored <- read_plink_table(parsed$args[4], c("FID", "IID"))
  variants <- readLines(parsed$args[5])
  if (!length(variants) || anyDuplicated(variants) || anyDuplicated(weights$ID) ||
      !setequal(variants, weights$ID)) {
    stop("PLINK projection did not score the full fitted marker set exactly once")
  }
  fit_rows <- match_id_rows(fam[ind.unrel, , drop = FALSE], fit, "PLINK eigenvectors")
  score_rows <- match_id_rows(fam[ind.keep, , drop = FALSE], scored, "PLINK scores")
  if (nrow(fit) != length(ind.unrel) || nrow(scored) != length(ind.keep)) {
    stop("PLINK fit or score sample counts differ from the requested sample sets")
  }
  if (!length(eigenval) || any(!is.finite(eigenval)) || max(eigenval) <= 0) {
    stop("PLINK PCA returned invalid eigenvalues")
  }
  axes <- which(eigenval > max(eigenval) * sqrt(.Machine$double.eps))
  k <- length(axes)
  if (k < 2L) stop("PLINK PCA returned fewer than two nondegenerate components")
  pc_columns <- paste0("PC", axes)
  fitted <- check_matrix(fit[fit_rows, pc_columns], length(ind.unrel), k, "PLINK eigenvectors")
  sums <- check_matrix(scored[score_rows, paste0(pc_columns, "_SUM")], length(ind.keep), k,
                       "PLINK projected scores")
  training <- sums[match(ind.unrel, ind.keep), , drop = FALSE]
  scale <- sqrt(colSums(training^2))
  agreement <- vapply(seq_len(k), function(j) cor(fitted[, j], training[, j]), numeric(1))
  if (any(!is.finite(scale) | scale == 0) || any(!is.finite(agreement) | agreement == 0)) {
    stop("PLINK projection returned a degenerate training axis")
  }
  if (any(abs(agreement) < 0.95)) {
    warning("Some PLINK projections agree poorly with the fitted eigenvectors. ",
            "Review the per-PC correlations in the PCA summary.")
  }
  PC <- sweep(sums, 2, scale * sign(agreement), "/")

  bim <- read.table(paste0(o$`bed-prefix`, ".bim"), header = FALSE, quote = "",
                    comment.char = "", colClasses = c("character", "character", "NULL",
                                                      "integer", "NULL", "NULL"))
  map_rows <- match(weights$ID, bim[[2]])
  if (anyNA(map_rows) || anyDuplicated(bim[[2]])) stop("Invalid marker IDs in the PCA BIM")
  load <- data.frame(CHR = bim[[1]][map_rows], POS = bim[[3]][map_rows], ID = weights$ID,
                     sweep(as.matrix(weights[pc_columns]), 2, sign(agreement), "*"),
                     check.names = FALSE)
  loadings <- check_matrix(load[pc_columns], nrow(load), k, "PLINK ALT weights")
  load$ALT_ALLELE <- weights$A1
  names(load)[match(pc_columns, names(load))] <- paste0("LOADING_PC", seq_len(k))
  retained <- weights$ID
  approximate <- length(ind.unrel) > 5000L &&
    length(retained) >= 2 * length(eigenval) * (length(eigenval) + 1)
  projection_label <- "PLINK 2 projection using unrelated training frequencies"
  variant_axis <- "LD-pruned variant index (input order)"
  method_summary <- c(list(
    n_variants_input = nrow(bim), n_variants_plink_retained = length(retained),
    method = "PLINK 2 PCA + variance-standardized scoring of all retained samples",
    plink_pca_algorithm = if (approximate) "approx" else "exact_mean_imputed",
    score_scale = "unit L2 norm in unrelated samples, signs aligned to fitted PCs"
  ), as.list(setNames(abs(agreement), paste0("fit_projection_correlation_PC", seq_len(k)))),
     as.list(setNames(eigenval[axes], paste0("eigenvalue_PC", seq_len(k)))))
} else {
  require_options(o, "lrldr-output")
  suppressPackageStartupMessages({
    library(bigsnpr)
    library(bigstatsr)
  })
  if (packageVersion("bigsnpr") < "1.12.14") stop("bigsnpr >= 1.12.14 is required for min.maf")
  obj <- bed(paste0(o$`bed-prefix`, ".bed"))
  if (nrow(fam) != nrow(obj)) stop("FAM/bed sample count mismatch")
  check_number(o$`autosvd-max-iter`, "autosvd-max-iter", 1, integer = TRUE)
  check_number(o$`autosvd-min-maf`, "autosvd-min-maf", 0, 0.5)
  k <- min(o$`n-pcs`, length(ind.unrel) - 1L, ncol(obj) - 1L)
  if (k < 2L) stop("Too few dimensions for PCA")

  svd <- bed_autoSVD(obj, ind.row = ind.unrel, k = k,
                     min.maf = o$`autosvd-min-maf`, max.iter = o$`autosvd-max-iter`, ncores = o$threads)
  retained <- attr(svd, "subset")
  if (is.null(retained) || !length(retained) || anyNA(retained) ||
      any(retained < 1 | retained > ncol(obj)) || anyDuplicated(retained)) {
    stop("autoSVD returned invalid retained-variant indices")
  }
  loadings <- check_matrix(svd$v, length(retained), k, "PCA loadings")
  lrldr <- attr(svd, "lrldr")
  if (is.null(lrldr)) lrldr <- data.frame(Chr = integer(), Start = integer(), Stop = integer(), Iter = integer())
  lrldr <- as.data.frame(lrldr)

  PC <- matrix(NA_real_, nrow = length(ind.keep), ncol = k)
  PC[match(ind.unrel, ind.keep), ] <- check_matrix(predict(svd), length(ind.unrel), k, "reference scores")
  if (length(proj.rows)) {
    pr <- bed_projectSelfPCA(svd, obj, ind.row = proj.rows, ind.col = retained, ncores = o$threads)
    PC[match(proj.rows, ind.keep), ] <- check_matrix(pr$OADP_proj, length(proj.rows), k, "OADP scores")
  }

  map <- obj$map
  load <- data.frame(CHR = map$chromosome[retained], POS = map$physical.pos[retained],
                     ID = as.character(map$marker.ID[retained]), loadings, check.names = FALSE)
  names(load)[seq.int(4L, ncol(load))] <- paste0("LOADING_PC", seq_len(k))
  write_tsv(lrldr, o$`lrldr-output`)
  projection_label <- "OADP projection"
  variant_axis <- "autoSVD-retained variant index (input order)"
  method_summary <- list(
    n_variants_input = ncol(obj), n_variants_autoSVD_retained = length(retained),
    n_lrld_regions = nrow(lrldr), autosvd_min_maf = o$`autosvd-min-maf`,
    autosvd_max_iter = o$`autosvd-max-iter`,
    bigsnpr_version = as.character(packageVersion("bigsnpr")),
    bigstatsr_version = as.character(packageVersion("bigstatsr")),
    method = "bed_autoSVD on unrelated + bed_projectSelfPCA OADP projection; review only"
  )
}
PC <- check_matrix(PC, length(ind.keep), k, "combined scores")
colnames(PC) <- paste0("PC", seq_len(k))
out <- data.frame(fam[ind.keep, , drop = FALSE],
                  PCA_ROLE = ifelse(ind.keep %in% ind.unrel, "UNRELATED_REFERENCE", "PROJECTED_RELATED"),
                  PC, check.names = FALSE)

# Preserve the supplied diagnostic: OGK is fitted to all retained samples.
# dist_ogk returns squared Mahalanobis distances, not ordinary distances.
if (!requireNamespace("bigutilsr", quietly = TRUE)) stop("Install bigutilsr in the R environment")
rk <- min(o$`robust-distance-n-pcs`, k)
distance_squared <- tryCatch(
  bigutilsr::dist_ogk(PC[, seq_len(rk), drop = FALSE]),
  error = function(e) stop("OGK diagnostic failed; review PC rank/sample composition: ", conditionMessage(e))
)
if (length(distance_squared) != nrow(PC)) stop("Unexpected OGK distance length")
r <- robust_threshold(distance_squared, o$`robust-mad-multiplier`)
rob <- data.frame(FID = out$FID, IID = out$IID,
                  ROBUST_LOG_DISTANCE = r$log_distance, ROBUST_DISTANCE_SQUARED = distance_squared,
                  THRESHOLD = r$threshold, SUGGESTED_OUTLIER = r$outlier)

write_tsv(out, o$`pcs-output`)
write_tsv(rob, o$`robust-output`)
write_tsv(rob[rob$SUGGESTED_OUTLIER, , drop = FALSE], o$`suggested-outliers`)
write_tsv(load, o$`loadings-output`)
writeLines(load$ID, o$`retained-variants`)

with_pdf(o$`scores-plot`, 8, 7, function() {
  for (i in pc_plot_pairs(colnames(PC), min(k, 10L))) {
    plot(PC[, i], PC[, i + 1L], pch = 19, cex = 0.45,
         col = ifelse(out$PCA_ROLE == "UNRELATED_REFERENCE", 1, 2),
         xlab = paste0("PC", i), ylab = paste0("PC", i + 1L),
         main = paste("Study PCA:", projection_label))
    legend("topright", legend = c("Unrelated reference", "Projected related"), col = c(1, 2), pch = 19, cex = 0.7)
  }
  hist(r$log_distance, breaks = "FD", main = "Robust PCA log-distance (diagnostic)",
       xlab = "Log squared robust Mahalanobis distance")
  abline(v = r$threshold, lty = 2)
})
with_pdf(o$`loadings-plot`, 11, 7, function() {
  for (i in seq_len(min(k, 10L))) {
    plot(seq_along(retained), loadings[, i], pch = 20, cex = 0.25,
         xlab = variant_axis, ylab = "Loading",
         main = paste0("PC", i, " loadings: inspect residual peaks"))
  }
})
write_summary(c(list(
  n_samples = length(ind.keep), n_unrelated_reference = length(ind.unrel),
  n_projected_related = length(proj.rows),
  n_pcs_requested = o$`n-pcs`, n_pcs = k, robust_distance_n_pcs = rk,
  robust_distance_threshold = r$threshold, robust_threshold_method = r$method,
  robust_log_distance_median = r$median, robust_log_distance_raw_MAD = r$mad,
  robust_mad_multiplier = o$`robust-mad-multiplier`, robust_fit_samples = "all_retained",
  n_zero_robust_distances = sum(distance_squared == 0), n_suggested_outliers = sum(r$outlier),
  seed = o$seed, bigutilsr_version = as.character(packageVersion("bigutilsr"))
), method_summary), o$`summary-output`)
