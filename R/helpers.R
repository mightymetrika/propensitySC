#' Compute match-frequency matrix and record PS values over k random splits
#' @param data data.frame of covariates (no response column)
#' @param rhs_formula one-sided formula for propensity score model, e.g. ~ X1 + X2
#' @param k number of iterations (default = min(100, sqrt(n)))
#' @param match_method matching method for MatchIt (e.g., "nearest")
#' @param distance_metric distance metric for MatchIt (e.g., "logit")
#' @return list(freq = n x n match-frequency matrix, ps_list = list of PS numeric vectors)
compute_match_frequency <- function(data,
                                    rhs_formula,
                                    k = min(100, floor(sqrt(nrow(data)))),
                                    match_method = "nearest",
                                    distance_metric = "logit") {
  n <- nrow(data)
  freq_mat <- matrix(0, n, n)
  ps_list  <- vector("list", k)

  for (i in seq_len(k)) {
    # 1. random binary assignment as pseudo-treatment
    data$part_k <- sample(c(rep(0, floor(n/2)), rep(1, ceiling(n/2))), n, replace = FALSE)
    # 2. build propensity-score formula by updating one-sided rhs_formula
    ps_formula <- update(rhs_formula, part_k ~ .)
    # 3. run matching
    m <- MatchIt::matchit(ps_formula,
                          data = data,
                          method = match_method,
                          distance = distance_metric)
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
compute_mean_ps <- function(ps_list) {
  ps_mat <- do.call(cbind, ps_list)
  rowMeans(ps_mat, na.rm = TRUE)
}

#' Compute pairwise absolute-distance matrix from a numeric vector
compute_dist_matrix <- function(x) {
  as.matrix(dist(x, method = "manhattan"))
}

#' Select initial anchors via farthest-first traversal on mean PS
#' @param meanPS numeric vector of length n
#' @param g_max number of anchors to select
#' @return integer vector of anchor indices
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

#' Create plot data for visualization (internal function)
#' @keywords internal
create_plot_data <- function(data, rhs_formula, cluster_labels, meanPS) {
  # Get model matrix for PCA
  model_matrix <- model.matrix(rhs_formula, data)

  # Remove intercept if present
  if (colnames(model_matrix)[1] == "(Intercept)") {
    model_matrix <- model_matrix[, -1, drop = FALSE]
  }

  # Perform PCA
  pca_result <- prcomp(model_matrix, scale. = TRUE, center = TRUE)

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

# # Perform iteration phase of propensity score matching
# perform_iterations <- function(data, formula, k, verbose) {
#   n <- nrow(data)
#
#   # Initialize storage for results
#   match_matrix <- matrix(0, n, n)
#   distance_matrix <- matrix(0, n, n)
#   match_counts <- matrix(0, n, n)
#
#   for (i in 1:k) {
#     if (verbose && i %% 10 == 0) cat("Iteration", i, "of", k, "\n")
#
#     # Random binary assignment
#     data$part_k <- sample(0:1, n, replace = TRUE)
#
#     # Run propensity score matching
#     ps_match <- MatchIt::matchit(as.formula(paste("part_k ~",
#                                                   as.character(formula)[2])),
#                                  data = data,
#                                  method = "nearest",
#                                  distance = "logit")
#
#     # Extract match data
#     match_data <- MatchIt::get_matches(ps_match)
#
#     # Update matrices based on matches
#     for (j in 1:nrow(match_data)) {
#       if (!is.na(match_data$subclass[j])) {
#         # Find matched pairs within same subclass
#         same_subclass <- which(match_data$subclass == match_data$subclass[j])
#         if (length(same_subclass) == 2) {
#           idx1 <- as.numeric(rownames(match_data)[same_subclass[1]])
#           idx2 <- as.numeric(rownames(match_data)[same_subclass[2]])
#
#           match_matrix[idx1, idx2] <- match_matrix[idx1, idx2] + 1
#           match_matrix[idx2, idx1] <- match_matrix[idx2, idx1] + 1
#
#           # Store distance (propensity score difference)
#           dist_val <- abs(match_data$distance[same_subclass[1]] -
#                             match_data$distance[same_subclass[2]])
#           distance_matrix[idx1, idx2] <- distance_matrix[idx1, idx2] + dist_val
#           distance_matrix[idx2, idx1] <- distance_matrix[idx2, idx1] + dist_val
#
#           match_counts[idx1, idx2] <- match_counts[idx1, idx2] + 1
#           match_counts[idx2, idx1] <- match_counts[idx2, idx1] + 1
#         }
#       }
#     }
#   }
#
#   # Average distances
#   distance_matrix <- distance_matrix / (match_counts + 1e-10)
#   distance_matrix[match_counts == 0] <- Inf
#
#   return(list(match_matrix = match_matrix,
#               distance_matrix = distance_matrix,
#               match_counts = match_counts))
# }
#
# # Select anchors that maximize minimum pairwise distance
# select_anchors <- function(match_matrix, distance_matrix, max_groups) {
#   n <- nrow(match_matrix)
#
#   # Find observations that were never matched to each other
#   never_matched <- which(rowSums(match_matrix == 0) > 0)
#
#   if (length(never_matched) < 2) {
#     stop("Not enough unmatched observations to form anchors")
#   }
#
#   # If we have more candidates than needed, select those with maximum min distance
#   if (length(never_matched) > max_groups) {
#     # Compute pairwise distances among candidates
#     candidate_dist <- distance_matrix[never_matched, never_matched]
#     diag(candidate_dist) <- Inf
#
#     # Greedy selection: start with pair with maximum distance
#     selected <- numeric(0)
#     remaining <- never_matched
#
#     # Find the pair with maximum distance
#     max_idx <- which(candidate_dist == max(candidate_dist[is.finite(candidate_dist)]),
#                      arr.ind = TRUE)[1,]
#     selected <- c(never_matched[max_idx[1]], never_matched[max_idx[2]])
#     remaining <- setdiff(remaining, selected)
#
#     # Iteratively add points that maximize minimum distance to selected set
#     while (length(selected) < max_groups && length(remaining) > 0) {
#       min_dists <- apply(distance_matrix[remaining, selected, drop = FALSE], 1, min)
#       next_idx <- remaining[which.max(min_dists)]
#       selected <- c(selected, next_idx)
#       remaining <- setdiff(remaining, selected)
#     }
#
#     return(selected)
#   } else {
#     return(never_matched)
#   }
# }
#
# # Perform iterative clustering
# iterative_clustering <- function(matching_results, anchors, min_improvement, verbose) {
#   best_score <- -Inf
#   best_clustering <- NULL
#   scores <- numeric(0)
#
#   current_anchors <- anchors
#   g <- length(anchors)
#
#   while (g >= 2) {
#     if (verbose) cat("Clustering with", g, "groups\n")
#
#     # Assign observations to groups
#     assignments <- assign_to_groups(matching_results, current_anchors)
#
#     # Compute quality score
#     score <- compute_quality_score(assignments,
#                                    matching_results$distance_matrix,
#                                    matching_results$match_matrix)
#     scores <- c(scores, score)
#
#     if (score > best_score) {
#       best_score <- score
#       best_clustering <- list(assignments = assignments,
#                               anchors = current_anchors,
#                               n_groups = g)
#     }
#
#     # Try removing each anchor and compute resulting quality
#     if (g > 2) {
#       removal_scores <- numeric(length(current_anchors))
#
#       for (i in 1:length(current_anchors)) {
#         temp_anchors <- current_anchors[-i]
#         temp_assignments <- assign_to_groups(matching_results, temp_anchors)
#         removal_scores[i] <- compute_quality_score(temp_assignments,
#                                                    matching_results$distance_matrix,
#                                                    matching_results$match_matrix)
#       }
#
#       # Find best removal
#       best_removal <- which.max(removal_scores)
#       improvement <- removal_scores[best_removal] - score
#
#       if (improvement < -min_improvement) {
#         if (verbose) cat("Early stopping: no significant improvement\n")
#         break
#       }
#
#       current_anchors <- current_anchors[-best_removal]
#       g <- g - 1
#     } else {
#       break
#     }
#   }
#
#   return(list(best_clustering = best_clustering,
#               scores = scores))
# }
#
# # Assign observations to groups based on matching frequency
# assign_to_groups <- function(matching_results, anchors) {
#   n <- nrow(matching_results$match_matrix)
#   assignments <- rep(NA, n)
#
#   # Anchors are assigned to themselves
#   for (i in 1:length(anchors)) {
#     assignments[anchors[i]] <- i
#   }
#
#   # Assign non-anchor observations
#   non_anchors <- setdiff(1:n, anchors)
#
#   for (obs in non_anchors) {
#     # Count matches with each anchor
#     match_counts <- matching_results$match_matrix[obs, anchors]
#
#     if (max(match_counts) == 0) {
#       # If no matches, assign to closest anchor by distance
#       distances <- matching_results$distance_matrix[obs, anchors]
#       assignments[obs] <- which.min(distances)
#     } else {
#       # Find anchors with maximum match count
#       max_matches <- which(match_counts == max(match_counts))
#
#       if (length(max_matches) == 1) {
#         assignments[obs] <- max_matches
#       } else {
#         # Break ties using distance
#         distances <- matching_results$distance_matrix[obs, anchors[max_matches]]
#         tie_breaker <- max_matches[which.min(distances)]
#         assignments[obs] <- tie_breaker
#       }
#     }
#   }
#
#   return(assignments)
# }
#
# # Compute clustering quality score
# compute_quality_score <- function(assignments, distance_matrix, match_matrix) {
#   groups <- unique(assignments[!is.na(assignments)])
#   n_groups <- length(groups)
#
#   if (n_groups < 2) return(0)
#
#   # Within-group quality (average distance within groups)
#   within_group_dist <- 0
#   n_within_pairs <- 0
#
#   for (g in groups) {
#     group_members <- which(assignments == g)
#     if (length(group_members) > 1) {
#       for (i in 1:(length(group_members)-1)) {
#         for (j in (i+1):length(group_members)) {
#           if (match_matrix[group_members[i], group_members[j]] > 0) {
#             within_group_dist <- within_group_dist +
#               distance_matrix[group_members[i], group_members[j]]
#             n_within_pairs <- n_within_pairs + 1
#           }
#         }
#       }
#     }
#   }
#
#   avg_within <- if (n_within_pairs > 0) within_group_dist / n_within_pairs else 1e-10
#
#   # Between-group quality (minimum distance between groups)
#   min_between <- Inf
#
#   for (g1 in 1:(n_groups-1)) {
#     for (g2 in (g1+1):n_groups) {
#       members1 <- which(assignments == groups[g1])
#       members2 <- which(assignments == groups[g2])
#
#       for (i in members1) {
#         for (j in members2) {
#           if (match_matrix[i, j] > 0) {
#             min_between <- min(min_between, distance_matrix[i, j])
#           }
#         }
#       }
#     }
#   }
#
#   if (!is.finite(min_between)) min_between <- 1
#
#   # Quality score: between-group separation / within-group compactness
#   score <- min_between / avg_within
#
#   return(score)
# }
#
# # Finalize results and create visualization
# finalize_results <- function(data, formula, clustering_results, matching_results) {
#   # Extract model matrix for PCA
#   model_matrix <- model.matrix(formula, data)
#
#   # Remove intercept if present
#   if (colnames(model_matrix)[1] == "(Intercept)") {
#     model_matrix <- model_matrix[, -1, drop = FALSE]
#   }
#
#   # Perform PCA
#   pca_result <- prcomp(model_matrix, scale. = TRUE, center = TRUE)
#
#   # Create plot
#   plot_data <- data.frame(
#     PC1 = pca_result$x[,1],
#     PC2 = pca_result$x[,2],
#     Cluster = factor(clustering_results$best_clustering$assignments)
#   )
#
#   # Create base R plot
#   plot(plot_data$PC1, plot_data$PC2,
#        col = as.numeric(plot_data$Cluster),
#        pch = 19,
#        xlab = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
#        ylab = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)"),
#        main = "Propensity Score Clustering Results")
#
#   # Add cluster centers
#   centers <- aggregate(cbind(PC1, PC2) ~ Cluster, data = plot_data, mean)
#   points(centers$PC1, centers$PC2, pch = 8, cex = 2, col = 1:nrow(centers))
#
#   # Add legend
#   legend("topright",
#          legend = paste("Cluster", levels(plot_data$Cluster)),
#          col = 1:length(levels(plot_data$Cluster)),
#          pch = 19)
#
#   # Return results
#   return(list(
#     cluster_assignments = clustering_results$best_clustering$assignments,
#     n_clusters = clustering_results$best_clustering$n_groups,
#     anchors = clustering_results$best_clustering$anchors,
#     quality_scores = clustering_results$scores,
#     final_score = max(clustering_results$scores),
#     pca_result = pca_result,
#     plot_data = plot_data
#   ))
# }
