#' Print method for ps_clustering objects
#'
#' Prints a concise summary of propensity score clustering results, including
#' the number of clusters found, quality score, cluster sizes, and anchor indices.
#'
#' @param x ps_clustering object
#' @param ... additional arguments (unused)
#'
#' @return Invisibly returns the input object.
#'
#' @examples
#' # Run clustering on iris data
#' result <- propensity_score_clustering(
#'   iris[, -5],
#'   ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width,
#'   verbose = FALSE
#' )
#'
#' # Print method is called automatically when object name is typed
#' result
#'
#' # Or explicitly call print
#' print(result)
#'
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
#'
#' Provides a detailed summary of propensity score clustering results, including
#' cluster sizes, quality metrics across different numbers of clusters, and
#' optional comparison with true labels if provided.
#'
#' @param object ps_clustering object
#' @param true_labels optional vector of true cluster labels for comparison
#' @param ... additional arguments (unused)
#'
#' @return Invisibly returns the input object.
#'
#' @details
#' When true labels are provided, the summary includes a confusion matrix,
#' Adjusted Rand Index (if the mclust package is available), and clustering
#' accuracy based on best label matching.
#'
#' @examples
#' # Run clustering on iris data
#' result <- propensity_score_clustering(
#'   iris[, -5],
#'   ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width,
#'   verbose = FALSE
#' )
#'
#' # Basic summary without true labels
#' summary(result)
#'
#' # Summary with comparison to true species labels
#' summary(result, true_labels = iris$Species)
#'
#' # Summary with custom true labels
#' custom_labels <- rep(c("A", "B"), each = 75)
#' summary(result, true_labels = custom_labels)
#'
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
#'
#' Creates visualizations for propensity score clustering results, including
#' PCA plots with cluster assignments, propensity score distributions,
#' quality score trajectories, and cluster size bar charts.
#'
#' @param x ps_clustering object
#' @param which integer vector specifying which plots to create (1-4):
#'   \itemize{
#'     \item 1: PCA plot with clusters and anchors
#'     \item 2: Boxplot of mean propensity scores by cluster
#'     \item 3: Quality score by number of clusters
#'     \item 4: Bar chart of cluster sizes
#'   }
#' @param ... additional arguments passed to plotting functions
#'
#' @return Invisibly returns the input object.
#'
#' @examples
#' # Run clustering on iris data
#' result <- propensity_score_clustering(
#'   iris[, -5],
#'   ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width,
#'   verbose = FALSE
#' )
#'
#' # Create only the PCA plot
#' plot(result, which = 1)
#'
#' @export
#' @method plot ps_clustering
plot.ps_clustering <- function(x, which = 1:4, ...) {
  # Save current par settings and restore on exit
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))

  # Determine layout based on number of plots
  n_plots <- length(which)
  if (n_plots == 1) {
    graphics::par(mfrow = c(1, 1))
  } else if (n_plots == 2) {
    graphics::par(mfrow = c(1, 2))
  } else if (n_plots <= 4) {
    graphics::par(mfrow = c(2, 2))
  } else {
    graphics::par(mfrow = c(ceiling(n_plots/2), 2))
  }

  # Create requested plots
  if (1 %in% which) {
    # 1. PCA plot colored by clusters
    # plot(x$plot_data$PC1, x$plot_data$PC2,
    #      col = as.numeric(x$plot_data$Cluster),
    #      pch = 19,
    #      xlab = paste0("PC1 (", round(summary(x$pca_result)$importance[2,1]*100, 1), "%)"),
    #      ylab = paste0("PC2 (", round(summary(x$pca_result)$importance[2,2]*100, 1), "%)"),
    #      main = "PCA: Clusters",
    #      ...)

    # Add cluster centers (anchors)
    # anchor_pcs <- x$plot_data[x$anchors, c("PC1", "PC2")]
    # graphics::points(anchor_pcs$PC1, anchor_pcs$PC2, pch = 8, cex = 2, lwd = 2)

    # anchor_pcs <- x$plot_data[x$anchors, c("PC1", "PC2")]
    # anchor_clusters <- x$plot_data$Cluster[x$anchors]
    # graphics::points(anchor_pcs$PC1, anchor_pcs$PC2,
    #                  pch = 8,
    #                  cex = 2,
    #                  lwd = 2,
    #                  col = as.numeric(anchor_clusters))  # Add color by cluster

    # anchor_pcs <- x$plot_data[x$anchors, c("PC1", "PC2")]
    # anchor_colors <- as.numeric(x$plot_data$Cluster)[x$anchors]
    # graphics::points(anchor_pcs$PC1, anchor_pcs$PC2,
    #                  pch = 8,
    #                  cex = 2,
    #                  lwd = 2,
    #                  col = anchor_colors)

    # anchor_pcs <- x$plot_data[x$anchors, c("PC1", "PC2")]
    # graphics::points(anchor_pcs$PC1, anchor_pcs$PC2,
    #                  pch = 8,
    #                  cex = 2,
    #                  lwd = 2,
    #                  col = as.numeric(x$plot_data$Cluster[x$anchors]))

    # anchor_pcs <- x$plot_data[x$anchors, c("PC1", "PC2")]
    # anchor_clusters <- x$plot_data$Cluster[x$anchors]
    #
    # # Let's see what values we're getting
    # print(anchor_clusters)  # This will show the factor levels
    # print(as.numeric(anchor_clusters))  # This will show the numeric codes
    #
    # # The fix: ensure we use the same color mapping as the main plot
    # graphics::points(anchor_pcs$PC1, anchor_pcs$PC2,
    #                  pch = 8,
    #                  cex = 2,
    #                  lwd = 2,
    #                  col = as.integer(factor(anchor_clusters)))

    # 1. PCA plot colored by clusters
    plot(x$plot_data$PC1, x$plot_data$PC2,
         col = as.numeric(x$plot_data$Cluster),
         pch = 19,
         xlab = paste0("PC1 (", round(summary(x$pca_result)$importance[2,1]*100, 1), "%)"),
         ylab = paste0("PC2 (", round(summary(x$pca_result)$importance[2,2]*100, 1), "%)"),
         main = "PCA: Clusters",
         ...)

    # Add cluster centers (anchors) with matching colors
    anchor_pcs <- x$plot_data[x$anchors, c("PC1", "PC2")]
    anchor_clusters <- x$plot_data$Cluster[x$anchors]
    graphics::points(anchor_pcs$PC1, anchor_pcs$PC2,
                     pch = 8,
                     cex = 2,
                     lwd = 2,
                     col = as.numeric(anchor_clusters))

    # Add legend
    graphics::legend("topright",
                     legend = paste("Cluster", levels(x$plot_data$Cluster)),
                     col = 1:length(levels(x$plot_data$Cluster)),
                     pch = 19,
                     cex = 0.8)
  }

  if (2 %in% which) {
    # 2. Mean propensity score by cluster
    graphics::boxplot(MeanPS ~ Cluster, data = x$plot_data,
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
    graphics::abline(v = x$n_clusters, col = "red", lty = 2)

    # Add text annotation for best
    best_idx <- which(x$quality_scores[, "g"] == x$n_clusters)
    graphics::text(x$n_clusters, x$quality_scores[best_idx, "overall"],
                   labels = paste("Best:", x$n_clusters),
                   pos = 3, col = "red")
  }

  if (4 %in% which) {
    # 4. Cluster sizes
    cluster_sizes <- table(x$cluster_assignments)
    graphics::barplot(cluster_sizes,
                      main = "Cluster Sizes",
                      xlab = "Cluster", ylab = "Number of Observations",
                      col = 1:x$n_clusters,
                      names.arg = paste("Cluster", names(cluster_sizes)),
                      ...)
  }

  invisible(x)
}
