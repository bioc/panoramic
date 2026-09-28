#' Internal: add directional colocalization labels
#'
#' @return A data frame with `coloc_source`, `coloc_target`, and
#'   `coloc_direction` columns added when the required cell-type columns exist.
#' @keywords internal
#' @noRd
.panoramic_add_coloc_direction <- function(df) {
  if (!is.data.frame(df)) return(df)
  if (!all(c("ct1", "ct2") %in% colnames(df))) return(df)
  if (!"coloc_source" %in% colnames(df)) {
    df$coloc_source <- as.character(df$ct2)
  }
  if (!"coloc_target" %in% colnames(df)) {
    df$coloc_target <- as.character(df$ct1)
  }
  if (!"coloc_direction" %in% colnames(df)) {
    df$coloc_direction <- paste0(df$coloc_source, " -> ", df$coloc_target)
  }
  df
}

#' Internal: flatten spatialstats assays to a long table
#'
#' @return A data frame containing feature identifiers and the `yi` and `vi`
#'   spatial-statistic assays.
#' @keywords internal
#' @noRd
.panoramic_extract_spatialstats_table <- function(se, drop_na = FALSE) {
  if (!inherits(se, "SummarizedExperiment")) {
    stop("`se` must be a SummarizedExperiment.", call. = FALSE)
  }
  if (!all(c("yi", "vi") %in% names(SummarizedExperiment::assays(se)))) {
    stop("`se` must contain assays `yi` and `vi`.", call. = FALSE)
  }

  yi <- SummarizedExperiment::assay(se, "yi")
  vi <- SummarizedExperiment::assay(se, "vi")
  rd <- S4Vectors::as.data.frame(SummarizedExperiment::rowData(se))
  cd <- S4Vectors::as.data.frame(SummarizedExperiment::colData(se))

  n_feat <- nrow(yi)
  n_samp <- ncol(yi)
  if (n_feat == 0L || n_samp == 0L) {
    out <- cbind(rd[0, , drop = FALSE], cd[0, , drop = FALSE])
    out$yi <- numeric(0)
    out$vi <- numeric(0)
    return(out)
  }

  rid <- rep(seq_len(n_feat), times = n_samp)
  cid <- rep(seq_len(n_samp), each = n_feat)

  feat_df <- rd[rid, , drop = FALSE]
  samp_df <- cd[cid, , drop = FALSE]

  out <- cbind(
    feat_df,
    samp_df,
    yi = as.numeric(yi),
    vi = as.numeric(vi)
  )

  if (isTRUE(drop_na)) {
    out <- out[is.finite(out$yi) & is.finite(out$vi), , drop = FALSE]
  }

  out <- .panoramic_add_coloc_direction(out)
  rownames(out) <- NULL
  out
}

#' Internal: extract feature-level meta-analysis table from rowData
#'
#' @return A data frame containing feature identifiers and meta-analysis
#'   statistics stored in `rowData(se)`.
#' @keywords internal
#' @noRd
.panoramic_extract_meta_table <- function(
    se,
    feature_cols = c("ct1", "ct2", "radius_um", "stat"),
    drop_na = FALSE
) {
  if (!inherits(se, "SummarizedExperiment")) {
    stop("`se` must be a SummarizedExperiment.", call. = FALSE)
  }

  rd <- S4Vectors::as.data.frame(SummarizedExperiment::rowData(se))
  if (nrow(rd) == 0L) return(rd)

  feature_cols <- intersect(feature_cols, colnames(rd))
  if (length(feature_cols) == 0L) {
    feature_cols <- character(0)
  }

  meta_cols <- setdiff(colnames(rd), feature_cols)
  out <- rd[, c(feature_cols, meta_cols), drop = FALSE]

  if (isTRUE(drop_na) && length(meta_cols) > 0L) {
    keep <- rowSums(!is.na(out[, meta_cols, drop = FALSE])) > 0L
    out <- out[keep, , drop = FALSE]
  }

  out <- .panoramic_add_coloc_direction(out)
  rownames(out) <- NULL
  out
}
