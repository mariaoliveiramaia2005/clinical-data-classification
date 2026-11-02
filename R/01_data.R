# 01_data.R
# Import data, select predictors and split training/test patients.

data <- read.csv(
  "data/MI.data",
  header = FALSE,
  na.strings = "?"
)

# Check the expected dataset structure.
stopifnot(nrow(data) == 1700, ncol(data) == 124)
stopifnot(
  !anyNA(data[[113]]),
  all(data[[113]] %in% c(0, 1))
)

# Question 1: information available at admission.
excluded_columns <- c(
  93, 94, 95, 100, 101, 102, 103, 104, 105
)

input_columns <- setdiff(2:112, excluded_columns)

var_info  <- read.csv("data/mi_variables.csv")
var_names <- var_info$name

stopifnot(ncol(data) == length(var_names))
names(data) <- var_names

# Blood pressures of zero are treated as missing for every patient.
# Retain the original zero locations for the exploratory audit.
bp_cols <- c("S_AD_KBRIG", "D_AD_KBRIG", "S_AD_ORIT", "D_AD_ORIT")
bp_zero <- !is.na(data[bp_cols]) & data[bp_cols] == 0
data[bp_cols][bp_zero] <- NA
cat("Dataset:", nrow(data), "patients x", ncol(data), "columns\n")
input_names   <- var_names[setdiff(2:112, excluded_columns)]
continuous_vars <- c("AGE", "S_AD_KBRIG", "D_AD_KBRIG", "S_AD_ORIT",
                     "D_AD_ORIT", "K_BLOOD", "NA_BLOOD", "ALT_BLOOD",
                     "AST_BLOOD", "KFK_BLOOD", "L_BLOOD", "ROE")
ordinal_vars    <- c("INF_ANAM", "STENOK_AN", "FK_STENOK", "IBS_POST", "GB",
                     "DLIT_AG", "ZSN_A", "ant_im", "lat_im", "inf_im",
                     "post_im", "TIME_B_S")
binary_vars     <- setdiff(input_names, c(continuous_vars, ordinal_vars))
stopifnot(all(c(continuous_vars, ordinal_vars) %in% input_names))  # sanity check



# Exclude the patient identifier and other outcomes from predictors.
admission_data <- data[, input_columns, drop = FALSE]

admission_data$FIBR_PREDS <- factor(
  data[[113]],
  levels = c(0, 1),
  labels = c("Nao", "Sim")
)

# Stratified split: 80% training, 20% test.
set.seed(42)

indices_by_class <- split(
  seq_len(nrow(admission_data)),
  admission_data$FIBR_PREDS
)

training_indices <- unlist(
  lapply(indices_by_class, function(indices) {
    sample(indices, size = floor(0.8 * length(indices)))
  }),
  use.names = FALSE
)

training <- admission_data[training_indices, , drop = FALSE]
test <- admission_data[-training_indices, , drop = FALSE]

cat("\nTraining dimensions:\n")
print(dim(training))

cat("\nTest dimensions:\n")
print(dim(test))

cat("\nTraining classes:\n")
print(table(training$FIBR_PREDS))

cat("\nTest classes:\n")
print(table(test$FIBR_PREDS))