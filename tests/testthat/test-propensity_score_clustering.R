test_that("propensity_score_clustering works", {

  # Generate data with true underlying clusters
  set.seed(123)
  # Run clustering
  result <- propensity_score_clustering(
    iris,
    ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width + Species,
    verbose = FALSE,
  ) |> suppressWarnings()

  # print(result)
  #
  # # Compare with known species
  # summary(result, true_labels = iris$Species)
  #
  # # Visualize clusters
  # plot(result, which = 1)  # PCA plot
  # plot(result, which = 2)  # PCA plot
  # plot(result, which = 3)  # PCA plot
  # plot(result, which = 4)  # PCA plot

  # Test 1: Check object structure and class
  expect_s3_class(result, "ps_clustering")

})

test_that("propensity_score_clustering works", {

  # Generate data with true underlying clusters
  set.seed(123)
  n <- 200
  n_groups <- 4
  group_size <- n/n_groups

  # True group assignments
  true_groups <- rep(1:n_groups, each = group_size)

  # Generate data with different characteristics for each group
  data <- data.frame(
    x1 = c(rnorm(group_size, mean = 0, sd = 1),
           rnorm(group_size, mean = 3, sd = 1),
           rnorm(group_size, mean = 0, sd = 3),
           rnorm(group_size, mean = 3, sd = 3)),
    x2 = c(rnorm(group_size, mean = -2, sd = 0.5),
           rnorm(group_size, mean = 2, sd = 0.5),
           rnorm(group_size, mean = 0, sd = 2),
           rnorm(group_size, mean = 1, sd = 2)),
    x3 = factor(c(sample(1:3, group_size, replace = TRUE, prob = c(0.7, 0.2, 0.1)),
                  sample(1:3, group_size, replace = TRUE, prob = c(0.1, 0.7, 0.2)),
                  sample(1:3, group_size, replace = TRUE, prob = c(0.2, 0.1, 0.7)),
                  sample(1:3, group_size, replace = TRUE, prob = c(0.33, 0.33, 0.34)))),
    x4 = c(rbinom(group_size, 1, 0.2),
           rbinom(group_size, 1, 0.8),
           rbinom(group_size, 1, 0.5),
           rbinom(group_size, 1, 0.5)),
    true_group = true_groups
  )

  # Run clustering
  suppressWarnings(
    result <- propensity_score_clustering(data[, -5], ~ x1 + x2 + x3 + x4, verbose = FALSE)
  )

  # Test 1: Check object structure and class
  expect_s3_class(result, "ps_clustering")

  # Test 2: Check that clustering produces valid results
  expect_equal(length(result$cluster_assignments), 200)
  expect_true(result$n_clusters >= 2 && result$n_clusters <= 14)  # Should be between 2 and max_groups

  # Test 3: Check that all observations are assigned to a cluster
  expect_true(all(result$cluster_assignments %in% 1:result$n_clusters))
})



# test_that("propensity_score_clustering works", {
#
#   # Generate data with true underlying clusters
#   set.seed(123)
#   n <- 200
#   n_groups <- 4
#   group_size <- n/n_groups
#
#   # True group assignments
#   true_groups <- rep(1:n_groups, each = group_size)
#
#   # Generate data with different characteristics for each group
#   data <- data.frame(
#     # Continuous variables with different means/variances by group
#     x1 = c(rnorm(group_size, mean = 0, sd = 1),    # Group 1
#            rnorm(group_size, mean = 3, sd = 1),    # Group 2
#            rnorm(group_size, mean = 0, sd = 3),    # Group 3
#            rnorm(group_size, mean = 3, sd = 3)),   # Group 4
#
#     x2 = c(rnorm(group_size, mean = -2, sd = 0.5), # Group 1
#            rnorm(group_size, mean = 2, sd = 0.5),  # Group 2
#            rnorm(group_size, mean = 0, sd = 2),    # Group 3
#            rnorm(group_size, mean = 1, sd = 2)),   # Group 4
#
#     # Factor variable with different probabilities by group
#     x3 = factor(c(sample(1:3, group_size, replace = TRUE, prob = c(0.7, 0.2, 0.1)), # Group 1
#                   sample(1:3, group_size, replace = TRUE, prob = c(0.1, 0.7, 0.2)), # Group 2
#                   sample(1:3, group_size, replace = TRUE, prob = c(0.2, 0.1, 0.7)), # Group 3
#                   sample(1:3, group_size, replace = TRUE, prob = c(0.33, 0.33, 0.34)))), # Group 4
#
#     # Binary variable with different probabilities by group
#     x4 = c(rbinom(group_size, 1, 0.2),  # Group 1
#            rbinom(group_size, 1, 0.8),  # Group 2
#            rbinom(group_size, 1, 0.5),  # Group 3
#            rbinom(group_size, 1, 0.5)), # Group 4
#
#     true_group = true_groups  # Keep for comparison
#   )
#
#   # Visualize the true structure
#   par(mfrow = c(2, 2))
#   plot(data$x1, data$x2, col = data$true_group, pch = 19, main = "True Groups: x1 vs x2")
#   boxplot(x1 ~ true_group, data = data, main = "x1 by True Group")
#   boxplot(x2 ~ true_group, data = data, main = "x2 by True Group")
#   barplot(table(data$true_group, data$x3), beside = TRUE,
#           main = "x3 by True Group", legend = TRUE)
#   par(mfrow = c(1, 1))
#
#   # Run clustering
#   result <- propensity_score_clustering(data[, -5], ~ x1 + x2 + x3 + x4)
#
#   # Print basic results
#   print(result)
#
#   # Get detailed summary
#   summary(result, true_labels = data$true_group)
#
#   # Create plots (all 4 by default)
#   plot(result)
#
#   # Or specific plots
#   plot(result, which = 1)  # Just PCA plot
#   plot(result, which = c(1, 3))  # PCA and quality score plots
#
#
#   expect_equal(2 * 2, 4)
# })
