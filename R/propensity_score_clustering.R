#' Enhanced propensity score clustering
#' @param data data.frame of covariates
#' @param rhs_formula one-sided formula for PS (e.g. ~ X1 + X2 + X3)
#' @param k number of random-matchit iterations
#' @param max_groups maximum initial clusters (default = sqrt(n))
#' @param min_improve early stopping threshold on quality drop
#' @param verbose logical, whether to print progress
#' @return ps_clustering object with clustering results
#' @export
propensity_score_clustering <- function(data,
                                        rhs_formula,
                                        k = min(100, floor(sqrt(nrow(data)))),
                                        max_groups = floor(sqrt(nrow(data))),
                                        min_improve = 0.01,
                                        verbose = TRUE) {
  n <- nrow(data)

  if (verbose) cat("Starting propensity score clustering with", k, "iterations\n")

  # 1) Build match-frequency and store PS draws
  tmp <- compute_match_frequency(data, rhs_formula, k,
                                 match_method = "nearest",
                                 distance_metric = "logit")
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

# File: R/methods.R (or could be in the same file)

#' Print method for ps_clustering objects
#' @param x ps_clustering object
#' @param ... additional arguments (unused)
#' @export
#' @method print ps_clustering
print.ps_clustering <- function(x, ...) {
  cat("Propensity Score Clustering Results\n")
  cat("==================================\n")
  cat("Number of clusters found:", x$n_clusters, "\n")
  cat("Best quality score:", round(x$best_score, 4), "\n")
  cat("\nCluster sizes:\n")
  print(table(x$cluster_assignments))
  cat("\nAnchors (indices):", x$anchors, "\n")

  invisible(x)
}

#' Summary method for ps_clustering objects
#' @param object ps_clustering object
#' @param true_labels optional vector of true cluster labels for comparison
#' @param ... additional arguments (unused)
#' @export
#' @method summary ps_clustering
summary.ps_clustering <- function(object, true_labels = NULL, ...) {
  cat("Propensity Score Clustering Summary\n")
  cat("==================================\n")
  cat("Number of observations:", length(object$cluster_assignments), "\n")
  cat("Number of clusters:", object$n_clusters, "\n")
  cat("Best quality score:", round(object$best_score, 4), "\n")

  cat("\nCluster sizes:\n")
  print(table(object$cluster_assignments))

  cat("\nQuality metrics by number of clusters:\n")
  print(round(object$quality_scores, 4))

  if (!is.null(true_labels)) {
    cat("\nComparison with true labels:\n")
    confusion <- table(True = true_labels, Predicted = object$cluster_assignments)
    print(confusion)

    # Calculate adjusted Rand index if mclust is available
    if (requireNamespace("mclust", quietly = TRUE)) {
      ari <- mclust::adjustedRandIndex(true_labels, object$cluster_assignments)
      cat("\nAdjusted Rand Index:", round(ari, 4), "\n")
    }

    # Add clustering accuracy (best matching)
    accuracy <- calculate_clustering_accuracy(true_labels, object$cluster_assignments)
    cat("Clustering accuracy (best matching):", round(accuracy, 4), "\n")
  }

  invisible(object)
}

#' Plot method for ps_clustering objects
#' @param x ps_clustering object
#' @param which integer vector specifying which plots to create (1-4)
#' @param ... additional arguments passed to plotting functions
#' @export
#' @method plot ps_clustering
plot.ps_clustering <- function(x, which = 1:4, ...) {
  # Save current par settings and restore on exit
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  # Determine layout based on number of plots
  n_plots <- length(which)
  if (n_plots == 1) {
    par(mfrow = c(1, 1))
  } else if (n_plots == 2) {
    par(mfrow = c(1, 2))
  } else if (n_plots <= 4) {
    par(mfrow = c(2, 2))
  } else {
    par(mfrow = c(ceiling(n_plots/2), 2))
  }

  # Create requested plots
  if (1 %in% which) {
    # 1. PCA plot colored by clusters
    plot(x$plot_data$PC1, x$plot_data$PC2,
         col = as.numeric(x$plot_data$Cluster),
         pch = 19,
         xlab = paste0("PC1 (", round(summary(x$pca_result)$importance[2,1]*100, 1), "%)"),
         ylab = paste0("PC2 (", round(summary(x$pca_result)$importance[2,2]*100, 1), "%)"),
         main = "PCA: Clusters",
         ...)

    # Add cluster centers (anchors)
    anchor_pcs <- x$plot_data[x$anchors, c("PC1", "PC2")]
    points(anchor_pcs$PC1, anchor_pcs$PC2, pch = 8, cex = 2, lwd = 2)

    # Add legend
    legend("topright",
           legend = paste("Cluster", levels(x$plot_data$Cluster)),
           col = 1:length(levels(x$plot_data$Cluster)),
           pch = 19,
           cex = 0.8)
  }

  if (2 %in% which) {
    # 2. Mean propensity score by cluster
    boxplot(MeanPS ~ Cluster, data = x$plot_data,
            main = "Mean Propensity Score by Cluster",
            xlab = "Cluster", ylab = "Mean PS",
            col = 1:x$n_clusters,
            ...)
  }

  if (3 %in% which) {
    # 3. Quality scores over iterations
    plot(x$quality_scores[, "g"], x$quality_scores[, "overall"],
         type = "b", pch = 19,
         xlab = "Number of Clusters",
         ylab = "Quality Score (Between/Within)",
         main = "Quality Score by Number of Clusters",
         ...)
    abline(v = x$n_clusters, col = "red", lty = 2)

    # Add text annotation for best
    best_idx <- which(x$quality_scores[, "g"] == x$n_clusters)
    text(x$n_clusters, x$quality_scores[best_idx, "overall"],
         labels = paste("Best:", x$n_clusters),
         pos = 3, col = "red")
  }

  if (4 %in% which) {
    # 4. Cluster sizes
    cluster_sizes <- table(x$cluster_assignments)
    barplot(cluster_sizes,
            main = "Cluster Sizes",
            xlab = "Cluster", ylab = "Number of Observations",
            col = 1:x$n_clusters,
            names.arg = paste("Cluster", names(cluster_sizes)),
            ...)
  }

  invisible(x)
}


# propensity_score_clustering <- function(data,
#                                         rhs_formula,
#                                         k = min(100, floor(sqrt(nrow(data)))),
#                                         max_groups = floor(sqrt(nrow(data))),
#                                         min_improve = 0.01) {
#   n <- nrow(data)
#   # 1) Build match-frequency and store PS draws
#   tmp <- compute_match_frequency(data, rhs_formula, k)
#   freq_mat <- tmp$freq
#   meanPS  <- compute_mean_ps(tmp$ps_list)
#   dist_mat <- compute_dist_matrix(meanPS)
#
#   # 2) Select initial anchors
#   anchors <- select_initial_anchors(meanPS, max_groups)
#
#   # 3) Iterative anchor removal
#   quality_results <- list()
#   best_overall <- -Inf
#   best_assign  <- NULL
#   best_g       <- length(anchors)
#   current_anchors <- anchors
#
#   for (g in seq(length(current_anchors), 2)) {
#     # assign clusters
#     assign <- assign_clusters(freq_mat, dist_mat, current_anchors)
#     # compute quality
#     q <- compute_quality(meanPS, current_anchors, assign)
#     quality_results[[as.character(g)]] <- c(g = g,
#                                             within = q$within,
#                                             between = q$between,
#                                             overall = q$overall)
#     if (q$overall > best_overall) {
#       best_overall <- q$overall
#       best_assign  <- assign
#       best_g       <- g
#     }
#     # try removing each anchor and see quality
#     if (g > 2) {
#       cand_scores <- sapply(seq_along(current_anchors), function(idx) {
#         trial_anchors <- current_anchors[-idx]
#         trial_assign  <- assign_clusters(freq_mat, dist_mat, trial_anchors)
#         q2 <- compute_quality(meanPS, trial_anchors, trial_assign)
#         q2$overall
#       })
#       best_idx <- which.max(cand_scores)
#       if ((cand_scores[best_idx] - q$overall) < -min_improve) break
#       current_anchors <- current_anchors[-best_idx]
#     }
#   }
#
#   list(
#     assignments   = best_assign,
#     quality_table = do.call(rbind, quality_results),
#     best_g        = best_g
#   )
# }


# propensity_score_clustering <- function(data,
#                                         formula,
#                                         k = NULL,
#                                         max_groups = NULL,
#                                         min_improvement = 0.01,
#                                         seed = NULL,
#                                         verbose = TRUE) {
#
#   # Set seed if provided
#   if (!is.null(seed)) set.seed(seed)
#
#   # Set defaults
#   n <- nrow(data)
#   if (is.null(k)) k <- min(100, ceiling(sqrt(n)))
#   if (is.null(max_groups)) max_groups <- ceiling(sqrt(n))
#
#   if (verbose) cat("Starting propensity score clustering with", k, "iterations\n")
#
#   # Step 1-2: Iteration Phase
#   matching_results <- perform_iterations(data, formula, k, verbose)
#
#   # Step 3: Anchor Selection
#   anchors <- select_anchors(matching_results$match_matrix,
#                             matching_results$distance_matrix,
#                             max_groups)
#
#   if (verbose) cat("Selected", length(anchors), "anchors\n")
#
#   # Step 4: Iterative Clustering
#   clustering_results <- iterative_clustering(matching_results,
#                                              anchors,
#                                              min_improvement,
#                                              verbose)
#
#   # Step 5: Final Selection and Visualization
#   final_results <- finalize_results(data,
#                                     formula,
#                                     clustering_results,
#                                     matching_results)
#
#   return(final_results)
# }
