# Shared validation and ID handling. Uses base R so helpers can be tested alone.
check_paths <- function(paths) {
  if (anyNA(paths) || any(!nzchar(paths)) || any(!file.exists(paths))) {
    stop("Missing input file(s): ", paste(paths[is.na(paths) | !nzchar(paths) | !file.exists(paths)], collapse = ", "))
  }
}

require_options <- function(o, keys) {
  for (key in keys) {
    value <- o[[key]]
    if (is.null(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
      stop("Missing required option --", key)
    }
  }
}

check_number <- function(x, label, lower, upper = Inf, integer = FALSE) {
  if (length(x) != 1L || !is.finite(x) || x < lower || x > upper ||
      (integer && x != floor(x))) stop("Invalid ", label)
}

id_key <- function(x) {
  if (!all(c("FID", "IID") %in% names(x)) || anyNA(x[, c("FID", "IID")]) ||
      any(!nzchar(x$FID)) || any(!nzchar(x$IID))) stop("Missing FID/IID values")
  paste(x$FID, x$IID, sep = "\t")
}

unique_ids <- function(x, label) {
  key <- id_key(x)
  if (anyDuplicated(key)) stop("Duplicate FID/IID pairs in ", label)
  key
}

read_fam_ids <- function(prefix) {
  path <- paste0(prefix, ".fam")
  check_paths(path)
  x <- read.table(path, header = FALSE, colClasses = "character", quote = "",
                  comment.char = "", na.strings = character(), stringsAsFactors = FALSE)
  if (ncol(x) != 6L || !nrow(x)) stop("Expected a nonempty six-column FAM: ", path)
  ids <- data.frame(FID = x[[1]], IID = x[[2]], stringsAsFactors = FALSE)
  unique_ids(ids, path)
  ids
}

read_id_indices <- function(path, fam) {
  check_paths(path)
  x <- read.table(path, header = FALSE, colClasses = "character", quote = "",
                  comment.char = "", na.strings = character(), stringsAsFactors = FALSE)
  if (!nrow(x)) stop("Empty sample-ID file: ", path)
  first <- sub("^#", "", as.character(x[1, ]))
  has_header <- first[1] %in% c("FID", "IID") && "IID" %in% first
  if (has_header) {
    if (anyDuplicated(first)) stop("Duplicate header columns: ", path)
    names(x) <- first
    x <- x[-1, , drop = FALSE]
    iid <- x$IID
    fid <- if ("FID" %in% names(x)) x$FID else NULL
  } else {
    if (!ncol(x) %in% 1:2) stop("Expected IID or FID/IID columns: ", path)
    iid <- x[[ncol(x)]]
    fid <- if (ncol(x) == 2L) x[[1]] else NULL
  }
  if (!length(iid)) stop("No samples in ID file: ", path)
  if (is.null(fid)) {
    if (anyDuplicated(fam$IID)) stop("IID-only list is ambiguous; provide FID and IID: ", path)
    if (anyDuplicated(iid)) stop("Duplicate IID in ", path)
    idx <- match(iid, fam$IID)
  } else {
    ids <- data.frame(FID = fid, IID = iid, stringsAsFactors = FALSE)
    idx <- match(unique_ids(ids, path), unique_ids(fam, "FAM"))
  }
  if (anyNA(idx)) stop("Sample IDs not found in FAM: ", path)
  sort(idx) # Output and PCA fit use FAM row order, regardless of list order.
}

match_id_rows <- function(query, target, label) {
  idx <- match(unique_ids(query, paste(label, "query")), unique_ids(target, label))
  if (anyNA(idx)) stop("Unmatched FID/IID pairs in ", label)
  idx
}

reference_label_rows <- function(fam, labels, fid_col, iid_col, group_col) {
  if (!all(c(iid_col, group_col) %in% names(labels))) stop("Reference labels lack IID/group columns")
  if (nzchar(fid_col) && fid_col %in% names(labels)) {
    ids <- data.frame(FID = labels[[fid_col]], IID = labels[[iid_col]], stringsAsFactors = FALSE)
    match(unique_ids(fam, "reference FAM"), unique_ids(ids, "reference labels"))
  } else {
    if (anyDuplicated(fam$IID) || anyDuplicated(labels[[iid_col]])) {
      stop("Reference IID matching is ambiguous; supply FID in the labels and configure fid_col")
    }
    match(fam$IID, labels[[iid_col]])
  }
}

normalize_build <- function(x) {
  alias <- c(b37 = "hg19", hg19 = "hg19", grch37 = "hg19", "37" = "hg19",
             b38 = "hg38", hg38 = "hg38", grch38 = "hg38", "38" = "hg38")
  value <- unname(alias[tolower(x)])
  if (length(value) != 1L || is.na(value)) stop("Unsupported genome build: ", x)
  value
}

check_matrix <- function(x, nr, nc, label) {
  x <- as.matrix(x)
  if (!identical(dim(x), c(as.integer(nr), as.integer(nc))) || any(!is.finite(x))) {
    stop("Unexpected dimensions or nonfinite values in ", label)
  }
  x
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
}

write_summary <- function(values, path) {
  write_tsv(data.frame(metric = names(values), value = as.character(unlist(values, use.names = FALSE))), path)
}

with_pdf <- function(path, width, height, draw) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(path, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  draw()
}

robust_threshold <- function(distance_squared, multiplier) {
  if (!length(distance_squared) || any(!is.finite(distance_squared)) || any(distance_squared < 0)) {
    stop("Invalid OGK squared distances; review PC rank and sample composition")
  }
  # A zero distance means a point at the estimated center. Clamp only for log/plotting.
  ldist <- log(pmax(distance_squared, .Machine$double.xmin))
  med <- median(ldist)
  mad_raw <- mad(ldist, center = med, constant = 1)
  fallback <- mad_raw == 0
  threshold <- if (fallback) unname(quantile(ldist, 0.999)) else med + multiplier * mad_raw
  list(log_distance = ldist, threshold = threshold, median = med, mad = mad_raw,
       method = if (fallback) "empirical_0.999_quantile_zero_MAD" else "median_plus_multiplier_raw_MAD",
       outlier = ldist > threshold)
}

pc_plot_pairs <- function(column_names, max_pcs = 6L) {
  available <- as.integer(sub("^PC", "", grep("^PC[0-9]+$", column_names, value = TRUE)))
  starts <- seq.int(1L, max_pcs - 1L, by = 2L)
  starts[starts %in% available & (starts + 1L) %in% available]
}

centroid_assignments <- function(ref_pc, new_pc, groups, probability, multiplier) {
  if (length(groups) != nrow(ref_pc) || anyNA(groups) || any(!nzchar(groups)) ||
      any(!is.finite(ref_pc)) || any(!is.finite(new_pc)) || ncol(ref_pc) != ncol(new_pc)) {
    stop("Invalid reference scores, target scores or group labels")
  }
  if ("UNASSIGNED" %in% groups) stop("UNASSIGNED is reserved; rename this reference group")
  counts <- table(groups)
  if (any(counts < 2L)) stop("Each reference group needs at least two samples to estimate a radius")
  levels <- sort(unique(groups))
  centers <- lapply(levels, function(g) colMeans(ref_pc[groups == g, , drop = FALSE]))
  radii <- vapply(seq_along(levels), function(j) {
    x <- ref_pc[groups == levels[j], , drop = FALSE]
    d <- sqrt(rowSums(sweep(x, 2, centers[[j]], "-")^2))
    unname(quantile(d, probability)) * multiplier
  }, numeric(1))
  if (any(!is.finite(radii)) || any(radii <= 0)) stop("Reference group radius is zero or nonfinite; review group sample sizes and PCA")
  # Keep memory proportional to sample count; no sample-by-group distance matrix.
  best <- rep(Inf, nrow(new_pc))
  index <- integer(nrow(new_pc))
  for (j in seq_along(levels)) {
    d <- sqrt(rowSums(sweep(new_pc, 2, centers[[j]], "-")^2))
    replace <- d < best
    index[replace] <- j
    best[replace] <- d[replace]
  }
  if (any(index == 0L)) stop("Nonfinite distance to all reference centroids")
  assigned <- best <= radii[index]
  list(rows = data.frame(REFERENCE_GROUP = ifelse(assigned, levels[index], "UNASSIGNED"),
                         NEAREST_REFERENCE_GROUP = levels[index], REFERENCE_DISTANCE = best,
                         REFERENCE_RADIUS = radii[index], ASSIGNED = assigned),
       groups = data.frame(GROUP = levels, N = as.integer(counts[levels]), RADIUS = radii))
}
