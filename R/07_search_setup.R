# Threshold grid: NA disables a filter.
filter_grid <- expand.grid(
  threshold_1 = c(NA_real_, 0.25, 0.5, 0.75),
  threshold_2 = c(NA_real_, 0.25),
  threshold_4 = c(NA_real_, 0.005, 0.01, 0.03),
  threshold_5 = c(NA_real_, 0.01, 0.02, 0.05),
  threshold_6 = c(NA_real_, 0.10, 0.50, 0.75),
  threshold_7 = c(NA_real_, 0.8, 0.85, 0.9),
  threshold_8 = c(NA_real_, 0.8, 0.85, 0.9),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

filter_grid$config_id <- seq_len(nrow(filter_grid))

print(nrow(filter_grid))

# Identify equivalent preprocessing within a fold.
preprocessing_key <- function(prep) {
  paste(
    paste(prep$fit_rows, collapse = ","),
    paste(prep$keep, collapse = ","),
    paste(prep$numeric_features, collapse = ","),
    paste(prep$categorical_features, collapse = ","),
    paste(prep$indicators, collapse = ","),
    sep = "|"
  )
}

# Save the grid to document the search.
write.csv(
  filter_grid,
  "outputs/report/planned_threshold_grid.csv",
  row.names = FALSE,
  na = "OFF"
)

# Define the original baseline thresholds.
baseline_thresholds <- c(
  threshold_1 = 0.5,
  threshold_2 = 0.5,
  threshold_4 = 0.01,
  threshold_5 = 0.02,
  threshold_6 = 0.5,
  threshold_7 = 0.85,
  threshold_8 = 0.85
)

# Original baseline, evaluated separately from the search grid.
baseline_configuration <- as.data.frame(
  as.list(baseline_thresholds)
)

baseline_configuration$config_id <- 0L

print(baseline_configuration)

predict_nb_configuration <- function(fit_data, validation_data, prep) {
  # Use fitting rows retained by filter 2.
  fit <- fit_data[prep$fit_rows, , drop = FALSE]

  x_fit <- apply_preprocessing(fit, prep)
  x_validation <- apply_preprocessing(validation_data, prep)

  category_levels <- fit_categories(x_fit, prep)

  x_fit <- apply_categories(x_fit, prep, category_levels)
  x_validation <- apply_categories(x_validation, prep, category_levels)

  model <- e1071::naiveBayes(
    x = x_fit,
    y = fit$FIBR_PREDS,
    laplace = 1
  )

  prob <- predict(
    model,
    newdata = x_validation,
    type = "raw"
  )[, "Sim"]

  stopifnot(
    length(prob) == nrow(validation_data),
    all(is.finite(prob))
  )

  prob
}

predict_knn_configuration <- function(
  fit_data, validation_data, prep, k_values
) {
  fit <- fit_data[prep$fit_rows, , drop = FALSE]

  x_fit <- apply_preprocessing(fit, prep)
  x_validation <- apply_preprocessing(validation_data, prep)

  # Learn scales using fitting data only.
  for (name in prep$numeric_features) {
    minimum <- min(x_fit[[name]])
    value_range <- max(x_fit[[name]]) - minimum

    stopifnot(is.finite(value_range), value_range > 0)

    x_fit[[name]] <- (x_fit[[name]] - minimum) / value_range
    x_validation[[name]] <- (
      x_validation[[name]] - minimum
    ) / value_range
  }

  distances <- gower_distance(
    x_fit, x_validation, prep$numeric_features
  )

  stopifnot(
    all(k_values >= 1),
    all(k_values <= nrow(fit))
  )

  # Sort neighbors once per patient.
  probabilities <- matrix(
    vapply(seq_len(nrow(validation_data)), function(i) {
      neighbor_order <- order(distances[i, ], seq_len(nrow(fit)))
      positive_labels <- as.integer(fit$FIBR_PREDS[neighbor_order] == "Sim")
      cumsum(positive_labels)[k_values] / k_values
    }, numeric(length(k_values))),
    nrow = nrow(validation_data), ncol = length(k_values), byrow = TRUE
  )

  colnames(probabilities) <- paste0("k", k_values)
  stopifnot(all(is.finite(probabilities)))

  probabilities
}

