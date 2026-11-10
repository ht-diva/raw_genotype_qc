#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(bigsnpr)
  library(bigstatsr)
})
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file)), "pca_helpers.R"))

o <- parse_args(OptionParser(option_list = list(
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
)))
outputs <- c("pcs-output", "scores-plot", "loadings-plot", "loadings-output", "lrldr-output",
             "robust-output", "suggested-outliers", "retained-variants", "summary-output")
require_options(o, c("bed-prefix", "unrelated-file", "keep-file", outputs))
check_number(o$`n-pcs`, "n-pcs", 2, integer = TRUE)
check_number(o$`robust-distance-n-pcs`, "robust-distance-n-pcs", 2, integer = TRUE)
check_number(o$`robust-mad-multiplier`, "robust-mad-multiplier", 0)
check_number(o$`autosvd-max-iter`, "autosvd-max-iter", 1, integer = TRUE)
check_number(o$`autosvd-min-maf`, "autosvd-min-maf", 0, 0.5)
check_number(o$threads, "threads", 1, integer = TRUE)
check_number(o$seed, "seed", 0, .Machine$integer.max, integer = TRUE)
if (packageVersion("bigsnpr") < "1.12.14") stop("bigsnpr >= 1.12.14 is required for min.maf")
if (!requireNamespace("bigutilsr", quietly = TRUE)) stop("Install bigutilsr in the R environment")
set.seed(o$seed)
check_paths(paste0(o$`bed-prefix`, c(".bed", ".bim", ".fam")))
for (path in unlist(o[outputs])) dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

# Read IDs as character strings, preserving leading zeroes and FID/IID identity.
obj <- bed(paste0(o$`bed-prefix`, ".bed"))
fam <- read_fam_ids(o$`bed-prefix`)
if (nrow(fam) != nrow(obj)) stop("FAM/bed sample count mismatch")
ind.keep <- read_id_indices(o$`keep-file`, fam)
ind.unrel <- read_id_indices(o$`unrelated-file`, fam)
if (any(!ind.unrel %in% ind.keep)) stop("Unrelated samples must be a subset of the retained samples")
if (length(ind.keep) < 10L || length(ind.unrel) < 10L) stop("Too few retained/unrelated samples for PCA (minimum 10)")
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
proj.rows <- setdiff(ind.keep, ind.unrel)
if (length(proj.rows)) {
  pr <- bed_projectSelfPCA(svd, obj, ind.row = proj.rows, ind.col = retained, ncores = o$threads)
  PC[match(proj.rows, ind.keep), ] <- check_matrix(pr$OADP_proj, length(proj.rows), k, "OADP scores")
}
PC <- check_matrix(PC, length(ind.keep), k, "combined scores")
colnames(PC) <- paste0("PC", seq_len(k))
out <- data.frame(fam[ind.keep, , drop = FALSE],
                  PCA_ROLE = ifelse(ind.keep %in% ind.unrel, "UNRELATED_REFERENCE", "PROJECTED_RELATED"),
                  PC, check.names = FALSE)

# Preserve the supplied diagnostic: OGK is fitted to all retained samples.
# dist_ogk returns squared Mahalanobis distances, not ordinary distances.
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

map <- obj$map
load <- data.frame(CHR = map$chromosome[retained], POS = map$physical.pos[retained],
                   ID = as.character(map$marker.ID[retained]), loadings, check.names = FALSE)
names(load)[seq.int(4L, ncol(load))] <- paste0("LOADING_PC", seq_len(k))
write_tsv(out, o$`pcs-output`)
write_tsv(rob, o$`robust-output`)
write_tsv(rob[rob$SUGGESTED_OUTLIER, , drop = FALSE], o$`suggested-outliers`)
write_tsv(load, o$`loadings-output`)
writeLines(load$ID, o$`retained-variants`)
write_tsv(lrldr, o$`lrldr-output`)

with_pdf(o$`scores-plot`, 8, 7, function() {
  for (i in pc_plot_pairs(colnames(PC), min(k, 10L))) {
    plot(PC[, i], PC[, i + 1L], pch = 19, cex = 0.45,
         col = ifelse(out$PCA_ROLE == "UNRELATED_REFERENCE", 1, 2),
         xlab = paste0("PC", i), ylab = paste0("PC", i + 1L),
         main = "Study PCA: unrelated reference + OADP projection")
    legend("topright", legend = c("Unrelated reference", "Projected related"), col = c(1, 2), pch = 19, cex = 0.7)
  }
  hist(r$log_distance, breaks = "FD", main = "Robust PCA log-distance (diagnostic)",
       xlab = "Log squared robust Mahalanobis distance")
  abline(v = r$threshold, lty = 2)
})
with_pdf(o$`loadings-plot`, 11, 7, function() {
  for (i in seq_len(min(k, 10L))) {
    plot(seq_along(retained), loadings[, i], pch = 20, cex = 0.25,
         xlab = "autoSVD-retained variant index (input order)", ylab = "Loading",
         main = paste0("PC", i, " loadings: inspect residual peaks"))
  }
})
write_summary(list(
  n_samples = length(ind.keep), n_unrelated_reference = length(ind.unrel),
  n_projected_related = length(proj.rows), n_variants_input = ncol(obj),
  n_variants_autoSVD_retained = length(retained), n_lrld_regions = nrow(lrldr),
  n_pcs_requested = o$`n-pcs`, n_pcs = k, autosvd_min_maf = o$`autosvd-min-maf`,
  autosvd_max_iter = o$`autosvd-max-iter`, robust_distance_n_pcs = rk,
  robust_distance_threshold = r$threshold, robust_threshold_method = r$method,
  robust_log_distance_median = r$median, robust_log_distance_raw_MAD = r$mad,
  robust_mad_multiplier = o$`robust-mad-multiplier`, robust_fit_samples = "all_retained",
  n_zero_robust_distances = sum(distance_squared == 0), n_suggested_outliers = sum(r$outlier),
  seed = o$seed, bigsnpr_version = as.character(packageVersion("bigsnpr")),
  bigstatsr_version = as.character(packageVersion("bigstatsr")),
  bigutilsr_version = as.character(packageVersion("bigutilsr")),
  method = "bed_autoSVD on unrelated + bed_projectSelfPCA OADP projection; review only"
), o$`summary-output`)
