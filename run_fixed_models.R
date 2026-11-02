# Portfolio entry point: refit four existing fixed configurations, without search.
# Original modelling functions are retained in R/ without changes.
required <- c("e1071", "pROC")
missing <- required[!vapply(required, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing)) stop("Install packages first: ", paste(missing, collapse = ", "))
if (!file.exists("data/MI.data")) stop("Place the original course MI.data in data/ first.")
dir.create("outputs/report", recursive = TRUE, showWarnings = FALSE)
source("R/01_data.R", encoding = "UTF-8")
source("R/03_functions.R", encoding = "UTF-8")
source("R/07_search_setup.R", encoding = "UTF-8")
training <- apply_log_transform(training)
test <- apply_log_transform(test)
configurations <- read.csv("results/configurations.csv", na.strings = c("OFF", "NA"))
threshold_names <- names(baseline_thresholds)
rows <- list()
for (i in seq_len(nrow(configurations))) {
  cfg <- configurations[i, ]
  prep <- do.call(fit_preprocessing, c(list(fit_data = training,
    use_indicators = cfg$model == "NB"), as.list(cfg[, threshold_names])))
  prob <- if (cfg$model == "NB") {
    predict_nb_configuration(training, test, prep)
  } else {
    as.numeric(predict_knn_configuration(training, test, prep, cfg$k))
  }
  roc <- pROC::roc(test$FIBR_PREDS, prob, levels = c("Nao", "Sim"), direction = "<", quiet = TRUE)
  for (rule in c("fixed_0.5", "youden")) {
    threshold <- if (rule == "fixed_0.5") 0.5 else cfg$youden_threshold
    rows[[length(rows) + 1L]] <- data.frame(model = cfg$model, scenario = cfg$scenario,
      rule = rule, threshold = threshold, AUC = as.numeric(pROC::auc(roc)),
      as.list(calculate_metrics(test$FIBR_PREDS, prob, threshold)))
  }
}
metrics <- do.call(rbind, rows)
write.csv(metrics, "outputs/report/refitted_test_metrics.csv", row.names = FALSE)
print(metrics)
message("Fixed settings retained. Do not retune on these test results.")
