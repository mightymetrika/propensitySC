#' Propensity Score Clustering
#'
#' Performs clustering by iteratively creating propensity scores from random balanced
#' binary partitions and building a match frequency matrix. The algorithm randomly
#' splits observations into two equal-sized groups, fits propensity score models,
#' and records matching patterns across iterations. Cluster anchors are identified
#' using farthest-first traversal on mean propensity scores, then observations are
#' assigned to clusters based on matching frequency and propensity score distance.
#' The number of clusters is determined by optimizing a quality metric that balances
#' within-cluster compactness and between-cluster separation.
#'
#' @param data A data frame containing covariates for clustering.
#' @param rhs_formula A one-sided formula specifying the right-hand side of the
#'   propensity score model (e.g., \code{~ X1 + X2 + X3}).
#' @param k Number of random matching iterations (default: \code{min(100, floor(sqrt(nrow(data))))}).
#' @param max_groups Maximum number of initial clusters (default: \code{floor(sqrt(nrow(data)))}).
#' @param min_improve Early stopping threshold for quality improvement (default: 0.01).
#' @param verbose Logical indicating whether to print progress messages (default: FALSE).
#' @param ... Additional arguments passed to MatchIt::matchit()
#'
#' @return An object of class \code{ps_clustering} containing:
#'   \item{cluster_assignments}{Integer vector of cluster assignments for each observation.}
#'   \item{n_clusters}{Number of clusters found.}
#'   \item{anchors}{Indices of anchor observations for each cluster.}
#'   \item{quality_scores}{Matrix of quality metrics for each number of clusters tested.}
#'   \item{best_score}{The highest quality score achieved.}
#'   \item{mean_propensity_scores}{Average propensity score for each observation across iterations.}
#'   \item{match_frequency_matrix}{Matrix of match frequencies between observations.}
#'   \item{plot_data}{Data frame containing PCA coordinates and cluster assignments for visualization.}
#'   \item{pca_result}{PCA results used for visualization.}
#'   \item{data}{Original input data.}
#'   \item{formula}{Formula used for propensity score modeling.}
#'
#' @examples
#' # Cluster the iris dataset using all measurements
#' result <- propensity_score_clustering(
#'   iris[, -5],
#'   ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width,
#'   verbose = FALSE
#' )
#' print(result)
#'
#' # Compare with known species
#' summary(result, true_labels = iris$Species)
#'
#' # Visualize clusters
#' plot(result, which = 1)  # PCA plot
#'
#' @export
propensity_score_clustering <- function(data,
                                        rhs_formula,
                                        k = min(100, floor(sqrt(nrow(data)))),
                                        max_groups = floor(sqrt(nrow(data))),
                                        min_improve = 0.01,
                                        verbose = FALSE, ...) {
  n <- nrow(data)

  if (verbose) cat("Starting propensity score clustering with", k, "iterations\n")

  # 1) Build match-frequency and store PS draws
  tmp <- compute_match_frequency(data, rhs_formula, k,
                                 match_method = "nearest",
                                 distance_metric = "logit",
                                 ...)
  freq_mat <- tmp$freq
  meanPS <- compute_mean_ps(tmp$ps_list)
  dist_mat <- compute_dist_matrix(meanPS)

  if (verbose) {
    cat("\nMatch matrix properties:\n")
    cat("Proportion of zero entries:", mean(freq_mat == 0), "\n")
    cat("Max matches between any pair:", max(freq_mat), "\n")
    cat("Average matches per observation:", mean(rowSums(freq_mat > 0)), "\n")
  }

  # 2) Select initial anchors
  anchors <- select_initial_anchors(meanPS, max_groups)

  if (verbose) cat("\nSelected", length(anchors), "initial anchors\n")

  # 3) Iterative anchor removal
  quality_results <- list()
  best_overall <- -Inf
  best_assign <- NULL
  best_g <- length(anchors)
  best_anchors <- anchors
  current_anchors <- anchors

  for (g in seq(length(current_anchors), 2)) {
    if (verbose) cat("Clustering with", g, "groups\n")

    # assign clusters
    assign <- assign_clusters(freq_mat, dist_mat, current_anchors)

    # compute quality
    q <- compute_quality(meanPS, current_anchors, assign)
    quality_results[[as.character(g)]] <- c(g = g,
                                            within = q$within,
                                            between = q$between,
                                            overall = q$overall)

    if (q$overall > best_overall) {
      best_overall <- q$overall
      best_assign <- assign
      best_g <- g
      best_anchors <- current_anchors
    }

    # try removing each anchor and see quality
    if (g > 2) {
      cand_scores <- sapply(seq_along(current_anchors), function(idx) {
        trial_anchors <- current_anchors[-idx]
        trial_assign <- assign_clusters(freq_mat, dist_mat, trial_anchors)
        q2 <- compute_quality(meanPS, trial_anchors, trial_assign)
        q2$overall
      })
      best_idx <- which.max(cand_scores)

      if (verbose) {
        improvement <- cand_scores[best_idx] - q$overall
        cat("  Best removal improves quality by:", round(improvement, 4), "\n")
      }

      if ((cand_scores[best_idx] - q$overall) < -min_improve) {
        if (verbose) cat("Early stopping: no significant improvement\n")
        break
      }
      current_anchors <- current_anchors[-best_idx]
    }
  }

  # Convert cluster assignments to consecutive integers
  cluster_labels <- as.integer(factor(best_assign))

  # Create plot data for later visualization
  plot_info <- create_plot_data(data, rhs_formula, cluster_labels, meanPS)

  # Create result object with S3 class
  result <- structure(
    list(
      cluster_assignments = cluster_labels,
      n_clusters = best_g,
      anchors = best_anchors,
      quality_scores = do.call(rbind, quality_results),
      best_score = best_overall,
      mean_propensity_scores = meanPS,
      match_frequency_matrix = freq_mat,
      plot_data = plot_info$data,
      pca_result = plot_info$pca_result,
      data = data,
      formula = rhs_formula
    ),
    class = "ps_clustering"
  )

  return(result)
}
