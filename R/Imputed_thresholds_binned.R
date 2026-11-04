# Read a long-form multiply imputed data file and separate it into multiple datasets,
# and then create a table of binary values for each dataset.
# Each column contains binary values indicating whether or not a pre-selected threshold is met.

# Load data.table package to read large CSV file
library(data.table)

# Turn off scientific notation
options(scipen = 999)

# State path of input file
IMPUTED_LONG_PATH <- "data/Imputed_dataset_long.csv"

# State path of output directory (create if not there already)
BINNED_DIR <- file.path("Python", "Imputed_dataset_binned")
if (!dir.exists(BINNED_DIR)) dir.create(BINNED_DIR, recursive = TRUE)

# ---- Function 1: Read a long-form multiply imputed CSV file ----
load_imputed_long_dataset <- function(long_path, n_imputations = 20L) {
  # Read CSV directly
  mi_dataset <- fread(long_path, na.strings = c("", "NA"))

  # If imputation_id is not already present, infer it from row order
  if (!"imputation_id" %in% names(mi_dataset)) {
    n_rows <- nrow(mi_dataset)
    if (n_rows %% n_imputations != 0) {
      stop(sprintf(
        "Expected the imputed file to be divisible by %d rows, but found %d rows.",
        n_imputations, n_rows
      ))
    }
    n_patients <- as.integer(n_rows / n_imputations)
    mi_dataset[, imputation_id := rep(seq_len(n_imputations), each = n_patients)]
  }

  # If patient_id is not already present, create it
  if (!"patient_id" %in% names(mi_dataset)) {
    n_patients <- nrow(mi_dataset) / n_imputations
    mi_dataset[, patient_id := rep(seq_len(n_patients), times = n_imputations)]
  }

  # Split into separate datasets by imputation
  split(mi_dataset, by = "imputation_id", keep.by = FALSE)
}

# ---- Threshold definitions ----
# Provide threshold names, symbols, and values
reps_of_each <- c(3, 1, 3, 3, 1, 3, 4, 3, 4, 1)
threshold_names <- rep(c("Hb", "Respiratory_distress", "para_ul", "Temperature_C", "Coma",
                         "Lactate", "Glucose_mmol_L", "Bicarbonate_mEq_L", "BUN",
                         "Died_28D"),
                       times = reps_of_each)
threshold_symbols <- rep(c("<", "=", ">", "<", "=", ">", "<", "<", ">", "="),
                         times = reps_of_each)
threshold_values <- c(4, 5, 6, 1, 100000, 250000, 500000, 35, 37.5, 40, 1, 5, 7, 10, 2.2, 3,
                      3.5, 4, 10, 12, 15, 8, 10, 15, 20, 1)
num_thresholds <- length(threshold_names)

# Write full names of threshold columns
threshold_full <- character(num_thresholds)

for (i in seq_len(num_thresholds)) {
  threshold_full[i] <- sprintf("%s:%s%s", threshold_names[i], threshold_symbols[i],
                               threshold_values[i])
}

# ---- Function 4: Create a table of binary values for each dataset ----
bin_one_dataset <- function(dt) {
  dt <- copy(dt)

  # Change respiratory distress from TRUE/FALSE to 1/0
  if ("Respiratory_distress" %in% names(dt)) {
    dt[, `Respiratory_distress` := as.integer(`Respiratory_distress`)]
  }

  # Reduce to threshold columns
  keep_raw <- unique(c("patient_id", "imputation_id", threshold_names))
  keep_raw <- keep_raw[keep_raw %in% names(dt)]
  dt <- dt[, ..keep_raw]

  # Store patient IDs if available, otherwise use row numbers
  patient_ids <- if ("patient_id" %in% names(dt)) {
    dt$patient_id
  } else {
    seq_len(nrow(dt))
  }

  # Store imputation IDs if available, otherwise assign 1
  imputation_ids <- if ("imputation_id" %in% names(dt)) {
    dt$imputation_id
  } else {
    1L
  }

  # Create output table
  out <- data.table(patient_id = patient_ids, imputation_id = imputation_ids)

  # Loop through features of interest
  # Want to split their elements in two (0s & 1s) based on a condition
  for (i in seq_len(num_thresholds)) {
    x <- dt[[threshold_names[i]]]

    # Separate conditions according to symbol
    # If condition true then set current element of vector to be 1
    # Else set current element of vector to be 0
    v <- if (threshold_symbols[i] == "<") {
      as.integer(x < threshold_values[i])
    } else if (threshold_symbols[i] == ">") {
      as.integer(x > threshold_values[i])
    } else {
      as.integer(x == threshold_values[i])
    }

    # Handle missing values
    v[is.na(x)] <- NA_integer_

    # Set the vector to be a column of new data frame
    out[[threshold_full[i]]] <- v
  }

  # Remove any rows with missing threshold values
  out <- na.omit(out)
  out
}

# Call Function 1 (read the long-form imputed file)
mi_list <- load_imputed_long_dataset(
  long_path = IMPUTED_LONG_PATH,
  n_imputations = 20L
)

# Loop through each imputed dataset, bin it, and write it
for (imp in seq_along(mi_list)) {
  dt_imp <- mi_list[[imp]]

  # Call Function 4 (create table of binary values)
  binned_imp <- bin_one_dataset(dt_imp)

  # Build the input filename for the current imputation
  out_file <- file.path(BINNED_DIR, sprintf("Imputed_dataset_binned_%02d.csv", imp))
  fwrite(binned_imp, out_file)

  message("Saved: ", out_file)
}