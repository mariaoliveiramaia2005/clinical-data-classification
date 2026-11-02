# Train the preprocessing model using only the adjustment data

# transform variabales using log per EDA
log_variables <- c("ALT_BLOOD", "AST_BLOOD", "ROE", "L_BLOOD")
apply_log_transform <- function(data) {
  # Do nothing if the log was already applied to this data
  if (isTRUE(attr(data, "log_applied"))) {
    return(data)
  }

  for (name in intersect(log_variables, names(data))) {
    x <- data[[name]]
    # log is undefined for values <= 0.
    stopifnot(all(x > 0, na.rm = TRUE))
    data[[name]] <- log(x)
  }

  attr(data, "log_applied") <- TRUE # mark the data as transformed
  data
}

# Transform the training/test frames explicitly in setup.R. This file can also
# be sourced by the prediction-only loader, without importing or fitting data.

fit_preprocessing <- function(
  fit_data,
  use_indicators = TRUE,
  threshold_1 = 0.5,
  threshold_2 = 0.5,
  threshold_4 = 0.01,
  threshold_5 = 0.02,
  threshold_6 = 0.5,
  threshold_7 = 0.85,
  threshold_8 = 0.85,
  association_cache = NULL
) {
  candidates <- setdiff(names(fit_data), "FIBR_PREDS")


  # Features to keep initial vector
  keep <- candidates
  # Filter 1: percentage of empty cells by feature (numerical and categorical)
  # Empirical threshold_1 to be optimized (OFF, 0.25, 0.5, 0.75)

  missing_fraction <- sapply(
    fit_data[, candidates, drop = FALSE],
    function(x) mean(is.na(x))
  )

  removed_f1 <- character(0)

  if (!is.na(threshold_1)) {
    removed_f1 <- candidates[missing_fraction > threshold_1]
  }

  keep <- setdiff(candidates, removed_f1)

  numeric_features <- intersect(
    keep,
    continuous_vars
  )

  categorical_features <- setdiff(keep, numeric_features)

  # Filter 2: percentage of empty cells by line (adjustment data ONLY).
  # Computed only on the columns kept by filter 1.
  # They are NOT removed from the validation or test data: apply_preprocessing() never removes rows.
  # Empirical threshold_2 to be optimized (OFF, 0.25, 0.5, 0.75)

  row_missing_fraction <- rowMeans(is.na(fit_data[, keep, drop = FALSE]))

  if (is.na(threshold_2)) {
    fit_rows <- seq_len(nrow(fit_data))
  } else {
    fit_rows <- which(row_missing_fraction <= threshold_2)
  }

  removed_f2 <- nrow(fit_data) - length(fit_rows)

  # From here on, every statistic uses only the patients that were kept.
  fit_data <- fit_data[fit_rows, , drop = FALSE]

  if (length(unique(fit_data$FIBR_PREDS)) < 2) {
    stop("The row filter removed every fitting patient from one class.")
  }


  # Filter 3: constant features (numerical and categorical)
  n_unique <- sapply(
    fit_data[, keep, drop = FALSE],
    function(x) length(unique(x[!is.na(x)]))
  )

  removed_f3 <- keep[n_unique <= 1]

  keep <- keep[n_unique > 1]

  numeric_features <- intersect(
    keep,
    continuous_vars # it's in 01_data.R
  )

  categorical_features <- setdiff(keep, numeric_features)

  # Filter 4: MAD/AMP (numerical):
  # Empirical threshold_4 to be optimized (OFF, 0.005, 0.01, 0.03)

  dispersion_ratio <- sapply(
    fit_data[, numeric_features, drop = FALSE],
    function(x) {
      mad_abs <- mad(x, na.rm = TRUE)

      value_range <- diff(range(x, na.rm = TRUE))

      if (value_range == 0) {
        return(0)
      }

      return(mad_abs / value_range)
    }
  )

  removed_f4 <- character(0)

  if (!is.na(threshold_4)) {
    removed_f4 <- names(dispersion_ratio)[
      is.finite(dispersion_ratio) & dispersion_ratio <= threshold_4
    ]
  }

  numeric_features <- setdiff(numeric_features, removed_f4)

  # Filter 5: Cross Correlation Feature-Target (numerical)
  # Empirical threshold_5 to be optimized (OFF, 0.01, 0.02, 0.05)

  numeric_target <- ifelse(fit_data$FIBR_PREDS == "Sim", 1, 0)

  corr_feature_target <- vapply(numeric_features, function(name) {
    x <- fit_data[[name]]
    complete <- is.finite(x) & is.finite(numeric_target)

    if (
      sum(complete) < 3L ||
        length(unique(x[complete])) < 2L ||
        length(unique(numeric_target[complete])) < 2L
    ) {
      return(NA_real_)
    }

    abs(cor(x[complete], numeric_target[complete]))
  }, numeric(1))

  removed_f5 <- character(0)

  if (!is.na(threshold_5)) {
    removed_f5 <- names(corr_feature_target)[
      is.finite(corr_feature_target) &
        corr_feature_target <= threshold_5
    ]
  }

  numeric_features <- setdiff(numeric_features, removed_f5)

  # Filter 6: feature-target association (categorical)
  # Use the chi-square test, Fisher's exact test, or the chi-square test with simulation,
  # depending on the expected frequencies and the size of the table.
  # Empirical threshold_6 to be optimized (OFF, 0.75, 0.50, 0.10)

  # Preserve the random state of the rest of the project.
  had_seed <- exists(
    ".Random.seed",
    envir = .GlobalEnv,
    inherits = FALSE
  )

  if (had_seed) {
    previous_seed <- get(".Random.seed", envir = .GlobalEnv)
  }

  on.exit(
    {
      if (had_seed) {
        assign(".Random.seed", previous_seed, envir = .GlobalEnv)
      } else if (exists(
        ".Random.seed",
        envir = .GlobalEnv, inherits = FALSE
      )) {
        rm(list = ".Random.seed", envir = .GlobalEnv)
      }
    },
    add = TRUE
  )

  # Make simulation tests reproducible.
  set.seed(2027)

  p_values_cat <- vapply(categorical_features, function(name) {
    table_data <- table(
      fit_data[[name]],
      fit_data$FIBR_PREDS
    )

    if (nrow(table_data) < 2 || ncol(table_data) < 2) {
      return(1)
    }

    expected <- outer(
      rowSums(table_data),
      colSums(table_data)
    ) / sum(table_data)

    adequate_approximation <- all(expected >= 5)

    if (adequate_approximation) {
      return(chisq.test(table_data, correct = FALSE)$p.value)
    }

    if (all(dim(table_data) == c(2, 2))) {
      return(fisher.test(table_data)$p.value)
    }

    chisq.test(
      table_data,
      simulate.p.value = TRUE,
      B = 9999
    )$p.value
  }, numeric(1))

  removed_f6 <- character(0)

  if (!is.na(threshold_6)) {
    removed_f6 <- names(p_values_cat)[
      is.finite(p_values_cat) &
        p_values_cat >= threshold_6
    ]
  }

  categorical_features <- setdiff(categorical_features, removed_f6)

  # Filter 7: Cross correlation feature-feature (numerical)
  # Empirical threshold_7 to be optimized (OFF, 0.8, 0.85, 0.9)

  removed_f7 <- character(0)

  if (!is.na(threshold_7) && length(numeric_features) > 1L) {
    correlation_matrix <- cor(fit_data[, numeric_features, drop = FALSE], use = "pairwise.complete.obs")
    high_correlation_indices <- which(abs(correlation_matrix) > threshold_7 & lower.tri(correlation_matrix), arr.ind = TRUE)

    if (nrow(high_correlation_indices) > 0) {
      for (i in seq_len(nrow(high_correlation_indices))) {
        f1 <- rownames(correlation_matrix)[high_correlation_indices[i, 1]]
        f2 <- colnames(correlation_matrix)[high_correlation_indices[i, 2]]

        if (f1 %in% removed_f7 || f2 %in% removed_f7) {
          next
        }

        if ((f1 == "S_AD_ORIT" && f2 == "D_AD_ORIT") || (f1 == "D_AD_ORIT" && f2 == "S_AD_ORIT")) {
          next
        }

        relevance <- corr_feature_target[c(f1, f2)]

        # Keep both variables when correlations cannot be estimated.
        if (any(!is.finite(relevance))) {
          next
        }

        less_relevant_feature <- if (relevance[[1]] < relevance[[2]]) {
          f1
        } else {
          f2
        }
        removed_f7 <- c(removed_f7, less_relevant_feature)
      }
      removed_f7 <- unique(removed_f7)
    }
  }

  numeric_features <- setdiff(numeric_features, removed_f7)


  # Filter 8: association feature-feature (categorical), bias-corrected Cramer's V.
  # In each highly associated pair, the variable less associated with the
  # target is removed. Pairs are processed from the most to the least
  # associated, skipping pairs where one variable was already removed.
  # Empirical threshold_8 to be optimized (OFF, 0.8, 0.85, 0.9)

  removed_f8 <- character(0)

  if (!is.na(threshold_8) && length(categorical_features) > 1) {
    target_associations <- sapply(categorical_features, function(name) {
      cached_cramers_v(
        fit_data, name, "FIBR_PREDS", association_cache
      )
    })

    pairs <- combn(categorical_features, 2)

    pair_associations <- apply(pairs, 2, function(p) {
      cached_cramers_v(
        fit_data, p[1], p[2], association_cache
      )
    })

    for (j in order(pair_associations, decreasing = TRUE)) {
      if (pair_associations[j] <= threshold_8) {
        break
      }

      f1 <- pairs[1, j]
      f2 <- pairs[2, j]

      if (f1 %in% removed_f8 || f2 %in% removed_f8) {
        next
      }

      less_relevant_feature <- if (target_associations[[f1]] < target_associations[[f2]]) f1 else f2
      removed_f8 <- c(removed_f8, less_relevant_feature)
    }
  }

  categorical_features <- setdiff(categorical_features, removed_f8)


  # Update the features to keep vector
  keep <- c(numeric_features, categorical_features)

  # Stop with a clear message if there are no remaining predictors.
  if (length(keep) == 0L) {
    stop(
      "The filters removed every predictor from the fitting data."
    )
  }

  # imputation
  imputation_values <- setNames(
    lapply(keep, function(name) {
      x <- fit_data[[name]]

      if (name %in% numeric_features) {
        median(x, na.rm = TRUE)
      } else {
        frequencies <- table(x)
        as.numeric(names(frequencies)[which.max(frequencies)])
      }
    }),
    keep
  )

  indicators <- character(0)

  if (use_indicators) {
    fit_indicators <- create_missing_indicators(fit_data)

    indicators <- names(fit_indicators)[
      vapply(
        fit_indicators,
        function(x) length(unique(x)) > 1L,
        logical(1)
      )
    ]

    keep <- c(keep, indicators)
    categorical_features <- c(categorical_features, indicators)

    for (name in indicators) {
      frequencies <- table(fit_indicators[[name]])
      imputation_values[[name]] <- as.integer(
        names(frequencies)[which.max(frequencies)]
      )
    }
  }

  list(
    keep = keep,
    numeric_features = numeric_features,
    categorical_features = categorical_features,
    indicators = indicators,
    imputation_values = imputation_values,
    fit_rows = fit_rows,
    n_removed_rows = removed_f2,
    removed = list(
      filter_1 = removed_f1,
      filter_3 = removed_f3,
      filter_4 = removed_f4,
      filter_5 = removed_f5,
      filter_6 = removed_f6,
      filter_7 = removed_f7,
      filter_8 = removed_f8
    )
  )
}

# Calculate missingness indicators before imputation.
create_missing_indicators <- function(data) {
  ecg <- grep(
    "^(ritm_ecg_p_|n_r_ecg_p_|n_p_ecg_p_)",
    names(data),
    value = TRUE
  )

  stopifnot(
    length(ecg) == 26L,
    all(c("S_AD_ORIT", "D_AD_ORIT") %in% names(data))
  )

  data.frame(
    missing_icu_bp = as.integer(
      is.na(data$S_AD_ORIT) | is.na(data$D_AD_ORIT)
    ),
    missing_ecg = as.integer(
      rowSums(!is.na(data[, ecg, drop = FALSE])) == 0L
    ),
    row.names = row.names(data)
  )
}

# Apply the fill patterns you've learned to the adjustment data.
apply_preprocessing <- function(new_data, prep) {
  inputs <- setdiff(prep$keep, prep$indicators)
  result <- new_data[, inputs, drop = FALSE]

  # Impute using values learned from fitting patients.
  for (name in inputs) {
    x <- result[[name]]
    x[is.na(x)] <- prep$imputation_values[[name]]
    result[[name]] <- x
  }

  # Calculate indicators on original data before imputation.
  if (length(prep$indicators) > 0L) {
    indicators <- create_missing_indicators(new_data)
    result <- cbind(
      result,
      indicators[, prep$indicators, drop = FALSE]
    )
  }

  result[, prep$keep, drop = FALSE]
}

# Prepare categories with levels learned only during the fitting process.
fit_categories <- function(x_fit, prep) {
  lapply(
    x_fit[, prep$categorical_features, drop = FALSE],
    function(x) sort(unique(x))
  )
}

apply_categories <- function(x, prep, category_levels) {
  for (name in prep$categorical_features) {
    unknown <- !(x[[name]] %in% category_levels[[name]])

    x[[name]][unknown] <- prep$imputation_values[[name]]

    x[[name]] <- factor(
      x[[name]],
      levels = category_levels[[name]]
    )
  }

  x
}

calculate_metrics <- function(real, probability, threshold = 0.5) {
  predicted <- factor(
    ifelse(probability >= threshold, "Sim", "Nao"),
    levels = c("Nao", "Sim")
  )

  confusion_matrix <- table(
    Actual = factor(real, levels = c("Nao", "Sim")),
    Predicted = predicted
  )

  tp <- confusion_matrix["Sim", "Sim"]
  fn <- confusion_matrix["Sim", "Nao"]
  tn <- confusion_matrix["Nao", "Nao"]
  fp <- confusion_matrix["Nao", "Sim"]

  sensitivity <- tp / (tp + fn)
  specificity <- tn / (tn + fp)

  c(
    accuracy = (tp + tn) / sum(confusion_matrix),
    error_01 = (fp + fn) / sum(confusion_matrix),
    sensitivity = sensitivity,
    specificity = specificity,
    balanced_accuracy = (sensitivity + specificity) / 2,
    precision = if ((tp + fp) > 0) tp / (tp + fp) else NA_real_,
    F1 = 2 * tp / (2 * tp + fp + fn)
  )
}

# Gower distances with scales learned only during the fitting process.
gower_distance <- function(fit, validation, numeric_features) {
  stopifnot(identical(names(fit), names(validation)))
  stopifnot(!anyNA(fit), !anyNA(validation))

  distances <- matrix(
    0,
    nrow = nrow(validation),
    ncol = nrow(fit)
  )

  for (name in names(fit)) {
    if (name %in% numeric_features) {
      differences <- abs(outer(
        validation[[name]],
        fit[[name]],
        "-"
      ))

      # Limit the contribution to 1 for validation values
      # that exceed the range observed during adjustment.
      differences <- pmin(differences, 1)
    } else {
      differences <- outer(
        validation[[name]],
        fit[[name]],
        "!="
      )
    }

    distances <- distances + differences
  }

  distances / ncol(fit)
}

# k-NN score: proportion of positives among the k nearest neighbors.
predict_knn <- function(distances, fit_response, k) {
  stopifnot(
    k >= 1,
    k <= ncol(distances),
    length(fit_response) == ncol(distances)
  )

  vapply(seq_len(nrow(distances)), function(i) {
    # If there is a tie, break the tie based on the order in the adjustment.
    neighbors <- order(
      distances[i, ],
      seq_len(ncol(distances))
    )[seq_len(k)]

    mean(fit_response[neighbors] == "Sim")
  }, numeric(1))
}

# Convert the exclusion lists into a table.
exclusion_table <- function(prep, stage) {
  reasons <- prep$removed

  data.frame(
    stage = rep(stage, sum(lengths(reasons))),
    reason = rep(names(reasons), lengths(reasons)),
    variable = as.character(unlist(reasons, use.names = FALSE)),
    stringsAsFactors = FALSE
  )
}

# Bias-corrected Cramer's V (Bergsma, 2013) between two categorical variables.
# Missing values are excluded (pairwise complete observations).
corrected_cramers_v <- function(x, y) {
  table_data <- table(x, y)
  table_data <- table_data[rowSums(table_data) > 0, colSums(table_data) > 0, drop = FALSE]

  n <- sum(table_data)
  r <- nrow(table_data)
  k <- ncol(table_data)

  if (n < 2 || r < 2 || k < 2) {
    return(0)
  }

  chi2 <- suppressWarnings(chisq.test(table_data, correct = FALSE)$statistic)

  phi2_corr <- max(0, chi2 / n - (k - 1) * (r - 1) / (n - 1))
  r_corr <- r - (r - 1)^2 / (n - 1)
  k_corr <- k - (k - 1)^2 / (n - 1)
  denominator <- min(k_corr - 1, r_corr - 1)

  if (denominator <= 0) {
    return(0)
  }

  unname(sqrt(phi2_corr / denominator))
}

# Diagnose outliers in continuous features.
# Accept only imputed fitting data.
diagnose_outliers <- function(x_fit, numeric_features) {
  x <- as.matrix(x_fit[, numeric_features, drop = FALSE])

  stopifnot(
    ncol(x) >= 2L,
    nrow(x) > ncol(x),
    all(is.finite(x))
  )

  # Preserve the random state of the rest of the project.
  had_seed <- exists(".Random.seed", envir = .GlobalEnv)
  if (had_seed) {
    previous_seed <- get(".Random.seed", envir = .GlobalEnv)
  }

  on.exit({
    if (had_seed) {
      assign(".Random.seed", previous_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv)) {
      rm(list = ".Random.seed", envir = .GlobalEnv)
    }
  })

  set.seed(2026)
  robust_fit <- robustbase::covMcd(x, alpha = 0.75)

  distance_squared <- mahalanobis(
    x,
    center = robust_fit$center,
    cov = robust_fit$cov
  )

  # Approximate reference: 97.5th percentile of the chi-square distribution.
  cutoff_squared <- qchisq(0.975, df = ncol(x))

  list(
    fit = robust_fit,
    variables = numeric_features,
    table_data = data.frame(
      row_id = rownames(x_fit),
      robust_distance = sqrt(distance_squared),
      cutoff = sqrt(cutoff_squared),
      flagged = distance_squared > cutoff_squared
    )
  )
}

cached_cramers_v <- function(data, name1, name2, cache = NULL) {
  if (is.null(cache)) {
    return(corrected_cramers_v(data[[name1]], data[[name2]]))
  }

  pair <- sort(c(name1, name2))

  key <- paste(
    paste(rownames(data), collapse = ","),
    paste(pair, collapse = ":"),
    sep = "|"
  )

  if (!exists(key, envir = cache, inherits = FALSE)) {
    value <- corrected_cramers_v(
      data[[name1]],
      data[[name2]]
    )

    assign(key, value, envir = cache)
  }

  get(key, envir = cache, inherits = FALSE)
}
