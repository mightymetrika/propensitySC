#' Compute match-frequency matrix and record PS values over k random splits
#' @param data data.frame of covariates (no response column)
#' @param rhs_formula one-sided formula for propensity score model, e.g. ~ X1 + X2
#' @param k number of iterations (default = min(100, sqrt(n)))
#' @param match_method matching method for MatchIt (e.g., "nearest")
#' @param distance_metric distance metric for MatchIt (e.g., "logit")
#' @param ... pass additional arguments to MatchIt::matchit()
#'
#' @keywords internal
#'
#' @return list(freq = n x n match-frequency matrix, ps_list = list of PS numeric vectors)
compute_match_frequency <- function(data,
                                    rhs_formula,
                                    k = min(100, floor(sqrt(nrow(data)))),
                                    match_method = "nearest",
                                    distance_metric = "logit",
                                    ...) {
  n <- nrow(data)
  freq_mat <- matrix(0, n, n)
  ps_list  <- vector("list", k)

  for (i in seq_len(k)) {
    # 1. random binary assignment as pseudo-treatment
    data$part_k <- sample(c(rep(0, floor(n/2)), rep(1, ceiling(n/2))), n, replace = FALSE)
    # 2. build propensity-score formula by updating one-sided rhs_formula
    ps_formula <- stats::update(rhs_formula, part_k ~ .)
    # 3. run matching
    m <- MatchIt::matchit(ps_formula,
                          data = data,
                          method = match_method,
                          distance = distance_metric,
                          ...)
    # extract propensity scores and subclasses
    ps  <- m$distance            # ps score per unit
    sub <- m$subclass            # subclass label for each unit
    ps_list[[i]] <- ps
    # accumulate match frequencies: all units within same subclass get a +1
    labs <- unique(sub[!is.na(sub)])
    for (lab in labs) {
      idx <- which(sub == lab)
      freq_mat[idx, idx] <- freq_mat[idx, idx] + 1
    }
  }
  diag(freq_mat) <- 0
  list(freq = freq_mat, ps_list = ps_list)
}

#' Compute mean propensity score per observation
#'
#' @param ps_list A list of numeric vectors, each containing the propensity scores
#'   for all observations from one matching iteration.
#' @return A numeric vector of length \(n\) giving the mean propensity score for
#'   each observation (row-wise mean across all iterations).
#' @keywords internal
compute_mean_ps <- function(ps_list) {
  ps_mat <- do.call(cbind, ps_list)
  rowMeans(ps_mat, na.rm = TRUE)
}

#' Compute pairwise absolute-distance matrix from a numeric vector
#'
#' @param x A numeric vector of length \eqn{n}, e.g. mean propensity scores.
#' @return A numeric \eqn{n \times n} matrix where entry \eqn{[i,j]} is the Manhattan
#'   distance \eqn{|x[i] - x[j]|}.
#' @keywords internal
compute_dist_matrix <- function(x) {
  as.matrix(stats::dist(x, method = "manhattan"))
}

#' Select initial anchors via farthest-first traversal on mean PS
#' @param meanPS numeric vector of length n
#' @param g_max number of anchors to select
#' @return integer vector of anchor indices
#' @keywords internal
select_initial_anchors <- function(meanPS, g_max) {
  n <- length(meanPS)
  dist_mat <- compute_dist_matrix(meanPS)
  anchors <- integer(g_max)
  # start with the point farthest from the overall mean
  anchors[1] <- which.max(abs(meanPS - mean(meanPS)))
  for (i in 2:g_max) {
    remaining <- setdiff(seq_len(n), anchors[1:(i-1)])
    d_min <- sapply(remaining, function(j) min(dist_mat[j, anchors[1:(i-1)]], na.rm = TRUE))
    anchors[i] <- remaining[which.max(d_min)]
  }
  anchors
}

#' Assign clusters based on match-frequency and tie-break by distance
#' @param freq_mat n x n match-frequency matrix
#' @param dist_mat n x n distance matrix (e.g., from compute_dist_matrix)
#' @param anchors integer vector of anchor indices
#' @return integer vector length n: cluster label (1..g) indicating which anchor
#' @keywords internal
assign_clusters <- function(freq_mat, dist_mat, anchors) {
  n <- nrow(freq_mat)
  g <- length(anchors)
  assignments <- integer(n)
  for (i in seq_len(n)) {
    freqs <- freq_mat[i, anchors]
    best <- which(freqs == max(freqs))
    if (length(best) > 1) {
      dists <- dist_mat[i, anchors[best]]
      best <- best[which.min(dists)]
      if (length(best) > 1) best <- sample(best, 1)
    }
    assignments[i] <- best
  }
  assignments
}

#' Compute clustering quality: within-group and between-group metrics
#' @param meanPS numeric vector of length n
#' @param anchors integer vector of chosen anchor indices
#' @param assignments integer vector of cluster labels (1..g)
#' @return list(within, between, overall = between/within)
#' @keywords internal
compute_quality <- function(meanPS, anchors, assignments) {
  # within: avg |PS_i - PS_anchor_i|
  within_vals <- abs(meanPS - meanPS[anchors[assignments]])
  within <- mean(within_vals, na.rm = TRUE)
  # between: min distance among anchor PS values
  a_ps <- meanPS[anchors]
  pairwise <- abs(outer(a_ps, a_ps, "-"))
  between <- min(pairwise[lower.tri(pairwise)])
  list(within = within, between = between, overall = between / within)
}


#' Calculate clustering accuracy using best label matching
#' @param true_labels vector of true cluster labels
#' @param pred_labels vector of predicted cluster labels
#' @return numeric accuracy score
#' @keywords internal
calculate_clustering_accuracy <- function(true_labels, pred_labels) {
  confusion <- table(true_labels, pred_labels)

  # Hungarian algorithm would be ideal, but for simplicity use greedy matching
  accuracy <- 0
  used_pred <- c()

  for (true_cluster in unique(true_labels)) {
    available_pred <- setdiff(unique(pred_labels), used_pred)
    if (length(available_pred) > 0) {
      matches <- confusion[as.character(true_cluster), as.character(available_pred)]
      best_pred <- available_pred[which.max(matches)]
      accuracy <- accuracy + confusion[as.character(true_cluster), as.character(best_pred)]
      used_pred <- c(used_pred, best_pred)
    }
  }

  accuracy / length(true_labels)
}

#' Create plot data for visualization of clustering results
#'
#' @param data A data.frame of the original covariates used for clustering.
#' @param rhs_formula A one-sided formula object specifying which covariates
#'   were used in the propensity‐score model (e.g. `~ x1 + x2 + factor1`).
#' @param cluster_labels An integer vector or factor of cluster assignments
#'   for each observation (as returned by `propensity_score_clustering()`).
#' @param meanPS Numeric vector of mean propensity scores for each observation.
#' @return A list with two elements:
#'   \describe{
#'     \item{data}{A data.frame with columns `PC1`, `PC2`, `Cluster` (factor),
#'       and `MeanPS`, ready for plotting.}
#'     \item{pca_result}{The full `prcomp` object, in case you want custom charts.}
#'   }
#' @keywords internal
create_plot_data <- function(data, rhs_formula, cluster_labels, meanPS) {
  # Get model matrix for PCA
  model_matrix <- stats::model.matrix(rhs_formula, data)

  # Remove intercept if present
  if (colnames(model_matrix)[1] == "(Intercept)") {
    model_matrix <- model_matrix[, -1, drop = FALSE]
  }

  # Perform PCA
  pca_result <- stats::prcomp(model_matrix, scale. = TRUE, center = TRUE)

  # Create plot data
  list(
    data = data.frame(
      PC1 = pca_result$x[,1],
      PC2 = pca_result$x[,2],
      Cluster = factor(cluster_labels),
      MeanPS = meanPS
    ),
    pca_result = pca_result
  )
}
