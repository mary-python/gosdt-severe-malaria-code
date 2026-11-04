# MULTIPLE IMPUTATION - MISSINGNESS AT RANDOM (MAR)
# 20 imputed datasets
# 20 iterations per dataset
#
# MICE (Multivariate Imputation by Chained Equations) methods:
# - pmm = predictive mean matching, used for continuous/numeric variables
# - logreg = logistic regression, used for binary variables
# - polyreg = polytomous logistic regression, used for unordered categorical variables
#             with >2 levels
# - polr = proportional odds logistic regression, used for ordered categorical variables

# Load data.table package to read large CSV file
library(data.table)

# Load MICE package to implement the above methods
library(mice)

# Load dplyr package to filter the imputed dataset
library(dplyr)

# Turn off scientific notation
options(scipen = 999)

my_table <- fread("data/Combined_dataset copy.csv")

# Remove patients with gender not M/F
my_table <- my_table[SEX %in% c("M", "F")]

# Remove patients with unidentified outcome
my_table <- my_table[!is.na(Died_28D)]

# Remove patients older than 15 years
my_table <- my_table[Age_category != ">15"]

# Remove studies with fewer than 20 patients
my_table[, study_n := .N, by = STUDY]
my_table <- my_table[study_n >= 20]
my_table[, study_n := NULL]

# Remove irrelevant columns
my_table <- my_table[, c(2, 3, 4, 5, 8, 12, 13, 14, 15, 17, 21, 43, 44, 47, 66, 67, 69)]

# Rename column names with a space
colnames(my_table)[6] <- "Respiratory_distress"
colnames(my_table)[12] <- "Base_Excess_mmol_L"

# Make sure the imputed variables are numeric
my_table[, c("Bicarbonate_mEq_L", "Base_Excess_mmol_L", "Hb", "HCT", "BUN",
             "Creatinine_umol_L", "Glucose_mmol_L", "Lactate", "Temperature_C") :=
           lapply(.SD, as.numeric),
         .SDcols = c("Bicarbonate_mEq_L", "Base_Excess_mmol_L", "Hb", "HCT", "BUN",
                     "Creatinine_umol_L", "Glucose_mmol_L", "Lactate", "Temperature_C")]

# Log transform BUN and Creatinine before imputation
my_table[, `:=`(log_BUN = log1p(BUN), log_Creatinine_umol_L = log1p(Creatinine_umol_L),
                log_Lactate = log1p(Lactate), log_para_ul = log1p(para_ul))]

# Drop the raw versions so MICE uses the log versions instead
my_table[, c("BUN", "Creatinine_umol_L", "Lactate", "para_ul") := NULL]

# Make sure likely categorical variables are factors if they were imported as character
for (v in c("Respiratory_distress", "log_para_ul", "Coma", "Age_category")) {
  if (v %in% names(my_table) && is.character(my_table[[v]])) {
    my_table[, (v) := factor(get(v))]
  }
}

# Variables that are allowed to be imputed (those that contain at least one missing value)
missing_vars <- names(which(colSums(is.na(my_table)) > 0))

# Variables that should not be imputed
# Died_28D is the outcome
# STUDY and SITEID are identifiers / clustering variables
# Age_category is derived / ordinal, so it is left out here
dont_impute <- c("Died_28D", "STUDY", "SITEID", "Age_category")

# Variables to be imputed
impute_vars <- setdiff(missing_vars, dont_impute)

# Function 1: Choose a reasonable MICE method
set_mice_method <- function(x) {
  if (is.ordered(x)) {
    return("polr")
  }
  if (is.factor(x)) {
    if (nlevels(x) == 2) return("logreg")
    return("polyreg")
  }
  return("pmm")
}

# Set up method vector for MICE
method_vec <- make.method(my_table)
method_vec[] <- ""

for (v in impute_vars) {

  # Call Function 1 (choose reasonable MICE method)
  method_vec[v] <- set_mice_method(my_table[[v]])
}

# Set up predictor matrix for MICE
# Start with all variables allowed as predictors
pred_matrix <- make.predictorMatrix(my_table)
pred_matrix[,] <- 1
diag(pred_matrix) <- 0

# Do not use these as predictors
# Died_28D is the outcome
# STUDY and SITEID are identifiers / clustering variables
# Age_category is derived / ordinal, so it is left out here
pred_matrix[, c("Died_28D", "STUDY", "SITEID", "Age_category")] <- 0

# Impute data in the dataset
imputed_data <- mice(my_table, m = 20, maxit = 20, method = method_vec,
                     predictorMatrix = pred_matrix, donors = 10, seed = 123)

# Create folder for plots
plot_dir <- file.path("data", "Imputation_plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Check imputation convergence and chain behaviour via trace plots
# Map variable names to trace plot names
trace_name_map <- c(log_Creatinine_umol_L = "Creatinine [umol/L]",
                    Bicarbonate_mEq_L = "Bicarbonate [mEq/L]", log_BUN = "BUN [mmol/L]",
                    Base_Excess_mmol_L = "Base Excess [mmol/L]",
                    log_Lactate = "Lactate [mmol/L]",
                    HCT = "HCT [L/L]", log_para_ul = "Parasitaemia [/uL]",
                    Respiratory_distress = "Respiratory distress",
                    Glucose_mmol_L = "Glucose [mmol/L]", Hb = "Hb [g/dL]",
                    Temperature_C = "Temperature [C]", Coma = "Coma")

# Back-transform for the log-scale trace plots
trace_back_transform <- c(log_Creatinine_umol_L = "Creatinine", log_BUN = "BUN",
                          log_Lactate = "Lactate", log_para_ul = "Parasitaemia")

# Function 2: Create custom trace plots for a MIDS (Multiply Imputed Data Set)
plot_mids_custom <- function(x, y = NULL, theme = mice.theme(), layout = c(2, 3), type = "l",
                             col = 1:10, lty = 1, var_labels = NULL, ylab = "Value") {

  # Function 3: Create custom strip labels for the trace plots
  strip.combined <- function(which.given, which.panel, factor.levels, ...) {
    lab <- as.character(factor.levels[which.panel[which.given]])

    # Relabel the variable-name strip on the right
    if (!is.null(var_labels) && which.given == 2 && lab %in% names(var_labels)) {
      lab <- unname(var_labels[lab])
    }

    # Draw the left strip (mean / sd)
    if (which.given == 1) {
      lattice::panel.rect(0, 0, 1, 1, col = theme$strip.background$col, border = 1)
      lattice::panel.text(x = 0, y = 0.5, pos = 4, lab = lab)
    }

    # Draw the right strip (variable name)
    if (which.given == 2) {
      lattice::panel.text(x = 1, y = 0.5, pos = 2, lab = lab)
    }
  }

  # Check that the supplied object is a MIDS object
  if (!is.mids(x)) {
    stop("argument 'x' must be a 'mids' object", call. = FALSE)
  }

  # Check that convergence diagnostics are available
  if (is.null(x$chainMean)) {
    stop("no convergence diagnostics found", call. = FALSE)
  }

  # Extract chain means and standard deviations
  mn <- x$chainMean
  sm <- sqrt(x$chainVar)

  # Back-transform selected log-scale traces to original units
  if (length(trace_back_transform) > 0) {
    for (v in intersect(names(trace_back_transform), dimnames(mn)[[1]])) {
      y_mean <- mn[v, , ]
      mn[v, , ] <- expm1(y_mean)
      sm[v, , ] <- exp(y_mean) * sm[v, , ]   # delta-method approximation
    }
  }

  # Select subset of non-missing entries
  obs <- apply(!(is.nan(mn) | is.na(mn)), 1, all)
  varlist <- names(obs)[obs]

  # Create formula if not given in y
  if (missing(y)) {
    formula <- as.formula(paste0(paste0(varlist, collapse = "+"), "~.it|.ms"))
  } else {
    formula <- NULL

    # Plot all variables if y is NULL
    if (is.null(y)) {
      formula <- as.formula(paste0(paste0(varlist, collapse = "+"), "~.it|.ms"))
    }

    # Plot variables supplied by name
    if (is.character(y)) {
      formula <- if (length(y) == 1) {
        as.formula(paste0(y, "~.it|.ms"))
      } else {
        as.formula(paste0(paste0(y, collapse = "+"), "~.it|.ms"))
      }
    }

    # Plot variables supplied by position
    if (is.integer(y) || is.logical(y)) {
      vars <- varlist[y]
      formula <- if (length(vars) == 1) {
        as.formula(paste0(vars, "~.it|.ms"))
      } else {
        as.formula(paste0(paste0(vars, collapse = "+"), "~.it|.ms"))
      }
    }

    # Accept a user-supplied plotting formula
    if (is.null(formula)) {
      formula <- as.formula(y)
    }
  }

  # Extract the number of imputations and iterations
  m <- x$m
  it <- x$iteration

  # Rearrange chain summaries into plotting format
  mn <- matrix(aperm(mn[varlist, , , drop = FALSE], c(2, 3, 1)), nrow = m * it)
  sm <- matrix(aperm(sm[varlist, , , drop = FALSE], c(2, 3, 1)), nrow = m * it)

  # Create plotting dataset for all iterations and imputations
  adm <- expand.grid(seq_len(it), seq_len(m), c("mean", "sd"))
  data <- cbind(adm, rbind(mn, sm))
  colnames(data) <- c(".it", ".m", ".ms", varlist)

  # Remove local object so lattice uses the data column
  .m <- NULL
  rm(.m)

  # Create trace plot
  tp <- lattice::xyplot(x = formula, data = data, groups = .m, type = type, lty = lty,
                        col = col, layout = layout, between = list(x = 0.9, y = 0.4),
                        scales = list(y = list(relation = "free", alternating = FALSE,
                                               tick.number = 3),
                                      x = list(alternating = FALSE)),
                        as.table = TRUE, xlab = "Iteration", ylab = ylab,

                        # Call Function 3 (create custom strip labels)
                        strip = strip.combined, par.strip.text = list(lines = 0.5))

  # Apply the MICE plotting theme
  update(tp, par.settings = theme)
}

# Save trace plot
jpeg(file.path(plot_dir, "trace_plot_1.jpg"), width = 1800, height = 1400, res = 200)

# Call Function 2 (create custom trace plot)
print(plot_mids_custom(imputed_data, c("log_Creatinine_umol_L", "Bicarbonate_mEq_L", "log_BUN"),
                       var_labels = trace_name_map))
dev.off()

# Save trace plot
jpeg(file.path(plot_dir, "trace_plot_2.jpg"), width = 1800, height = 1400, res = 200)

# Call Function 2 (create custom trace plot)
print(plot_mids_custom(imputed_data, c("Base_Excess_mmol_L", "log_Lactate", "HCT"),
                       var_labels = trace_name_map))
dev.off()

# Save trace plot
jpeg(file.path(plot_dir, "trace_plot_3.jpg"), width = 1800, height = 1400, res = 200)

# Call Function 2 (create custom trace plot)
print(plot_mids_custom(imputed_data, c("log_para_ul", "Respiratory_distress", "Glucose_mmol_L"),
                       var_labels = trace_name_map))
dev.off()

# Save trace plot
jpeg(file.path(plot_dir, "trace_plot_4.jpg"), width = 1800, height = 1400, res = 200)

# Call Function 2 (create custom trace plot)
print(plot_mids_custom(imputed_data, c("Hb", "Temperature_C", "Coma"),
                       var_labels = trace_name_map))
dev.off()

# Check distributions of imputed variables
density_vars <- c("log_Creatinine_umol_L", "Bicarbonate_mEq_L", "log_BUN", "Base_Excess_mmol_L",
                  "log_Lactate", "HCT", "log_para_ul", "Respiratory_distress", "Glucose_mmol_L",
                  "Hb", "Temperature_C", "Coma")

# Density plot labels
density_labels <- c(log_Creatinine_umol_L = "Creatinine [umol/L]",
                    Bicarbonate_mEq_L = "Bicarbonate [mEq/L]", log_BUN = "BUN [mmol/L]",
                    Base_Excess_mmol_L = "Base Excess [mmol/L]",
                    log_Lactate = "Lactate [mmol/L]", HCT = "HCT [L/L]",
                    log_para_ul = "Parasitaemia [/uL]",
                    Respiratory_distress = "Respiratory distress",
                    Glucose_mmol_L = "Glucose [mmol/L]", Hb = "Hb [g/dL]",
                    Temperature_C = "Temperature [C]", Coma = "Coma")

# Density plot axis labels for the logged variables
density_ticks <- list(log_Creatinine_umol_L = c(1, 10, 25, 50, 100, 200, 400, 1000),
                      log_BUN = c(0, 1, 2.5, 5, 10, 25, 50, 100),
                      log_Lactate = c(0, 1, 2.5, 5, 10, 20, 40),
                      log_para_ul = c(0, 10, 1000, 100000, 10000000))

for (v in density_vars) {
  jpeg(file.path(plot_dir, paste0("density_", v, ".jpg")),
       width = 1800, height = 1400, res = 200)

  if (v %in% names(density_ticks)) {
    print(densityplot(imputed_data, as.formula(paste("~", v)),
                      scales = list(x = list(cex = 1.8,
                                             at = log1p(density_ticks[[v]]),
                                             labels = format(density_ticks[[v]],
                                                             big.mark = ",",
                                                             scientific = FALSE,
                                                             trim = TRUE)),
                                    y = list(cex = 1.8)),
                      xlab = list(label = density_labels[[v]], cex = 2.2),
                      ylab = list(label = "Density", cex = 2.2)))
  } else {
    print(densityplot(imputed_data, as.formula(paste("~", v)),
                      scales = list(x = list(cex = 1.8), y = list(cex = 1.8)),
                      xlab = list(label = density_labels[[v]], cex = 2.2),
                      ylab = list(label = "Density", cex = 2.2)))
  }

  dev.off()
}

# Create long-format imputed dataset
imputed_table <- complete(imputed_data, action = "long", include = TRUE)
setDT(imputed_table)

# Back-transform to original units
imputed_table[, `:=`(BUN = expm1(log_BUN), Creatinine_umol_L = expm1(log_Creatinine_umol_L),
                     Lactate = expm1(log_Lactate), para_ul = expm1(log_para_ul))]

# Drop the log columns to ensure the saved file stays in original units
imputed_table <- imputed_table %>%
  select(-log_BUN, -log_Creatinine_umol_L, -log_Lactate, -log_para_ul)

# Filter out the original data (.imp == 0)
imputed_table_filtered <- imputed_table %>%
  filter(.imp != 0)

# Save filtered imputed dataset
fwrite(imputed_table_filtered, "data/Imputed_dataset_long.csv")