library(data.table)
library(dplyr)
library(mgcv)

# Set the main folder and the output folders
DATA_FILE <- "data/Combined_dataset copy.csv"
OUT_DIR <- "data"

# Create the output folders if they do not already exist
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Read the dataset
Data <- fread(DATA_FILE)

# Use the death variable that is present in the file
outcome_col <- if ("Died_28D" %in% names(Data)) {
  "Died_28D"
} else {
  "Died within 28 days"
}

# Clean the key columns and make sure they are stored in the right form
Data <- Data %>%
  mutate(Died_28D = as.integer(.data[[outcome_col]] %in% c(1, "1", "Yes", "YES", "yes")),
         Age_category = factor(Age_category), Coma = factor(as.integer(Coma), levels = c(0, 1)),
         STUDY = factor(STUDY), SITEID = factor(SITEID))

# List the markers we want to study
marker_specs <- tibble::tribble(~marker, ~label, ~direction,
                                "Bicarbonate_mEq_L", "Bicarbonate", "lower",
                                "Temperature_C", "Temperature", "lower",
                                "Lactate", "Lactate", "upper",
                                "Glucose_mmol_L", "Glucose", "lower",
                                "Hb", "Hb", "lower")

# Only search for thresholds in the ranges that make sense for each marker
threshold_windows <- list(Bicarbonate_mEq_L = c(5, 25), Temperature_C = c(30, 42),
                          Lactate = c(0, 15), Glucose_mmol_L = c(0, 6), Hb = c(0, 8))

# Fit one model for one marker
fit_marker_gam <- function(marker, data) {
  fml <- as.formula(paste0("Died_28D ~ s(", marker,
                           ", bs = 'ts', by = Age_category, k = 100) + ",
                           "s(STUDY/SITEID, bs = 're', k = 100) + Coma"))

  gam(fml, family = binomial(link = "logit"), data = data, select = TRUE, method = "REML",
      na.action = na.exclude)
}

# Fit all five marker models
models <- setNames(lapply(marker_specs$marker, fit_marker_gam, data = Data),
                   marker_specs$marker)

# Get the study and site levels used in the fitted model
get_ref_levels <- function(model, data) {
  mf <- tryCatch(model.frame(model), error = function(e) NULL)

  if (!is.null(mf) && "STUDY" %in% names(mf) && !is.null(mf$STUDY)) {
    study_ref <- levels(droplevels(mf$STUDY))[1]
  } else {
    study_ref <- levels(data$STUDY)[1]
  }

  if (!is.null(mf) && "SITEID" %in% names(mf) && !is.null(mf$SITEID)) {
    site_ref <- levels(droplevels(mf$SITEID))[1]
  } else {
    site_ref <- levels(data$SITEID)[1]
  }

  list(study_ref = study_ref, site_ref = site_ref)
}

# Make a small prediction dataset for a given marker value
build_newdata <- function(data, marker, x, age_category, coma, study_ref, site_ref) {
  nd <- data.frame(Age_category = factor(age_category, levels = levels(data$Age_category)),
                   Coma = factor(coma, levels = levels(data$Coma)),
                   STUDY = factor(study_ref, levels = levels(data$STUDY)),
                   SITEID = factor(site_ref, levels = levels(data$SITEID)))
  nd <- nd[rep(1, length(x)), , drop = FALSE]
  nd[[marker]] <- x
  nd
}

# Get the allowed search range for each marker
marker_range <- function(marker, data) {
  threshold_windows[[marker]]
}

# Get age weights for coma = 0
age_weights <- function(data, coma = "0") {
  w <- prop.table(table(data$Age_category[data$Coma == coma]))
  w <- w[levels(data$Age_category)]
  w[is.na(w)] <- 0
  as.numeric(w)
}

# Draw the full risk curve for one marker and one coma group
pooled_curve <- function(model, data, marker, coma = "0", n = 5000) {
  refs <- get_ref_levels(model, data)
  x_rng <- marker_range(marker, data)
  x <- seq(x_rng[1], x_rng[2], length.out = n)
  age_levels <- levels(data$Age_category)
  w <- age_weights(data, coma)

  risk_mat <- sapply(age_levels, function(a) {
    nd <- build_newdata(data, marker, x, a, coma, refs$study_ref, refs$site_ref)
    as.numeric(predict(model, newdata = nd, type = "response"))
  })

  if (is.null(dim(risk_mat))) {
    risk_mat <- matrix(risk_mat, ncol = 1)
  }

  pooled_risk <- as.numeric(risk_mat %*% w / sum(w))

  data.frame(marker = marker, x = x, Coma = coma, risk = pooled_risk)
}

# Find where the curve crosses 5%
find_crossings <- function(x, p, target = 0.05) {
  s <- p - target
  roots <- numeric(0)

  # If any point lands exactly on 5%, keep it
  exact_hits <- which(abs(s) < 1e-8)
  if (length(exact_hits) > 0) {
    roots <- c(roots, x[exact_hits])
  }

  # Find places where the curve moves from above 5% to below 5%, or the other way around
  idx <- which(s[-length(s)] * s[-1] < 0)
  if (length(idx) > 0) {
    roots <- c(roots, mapply(function(i)
                 x[i] + (target - p[i]) * (x[i + 1] - x[i]) / (p[i + 1] - p[i]), idx))
  }

  sort(unique(round(roots, 8)))
}

# Find the 5% threshold from the pooled curve
find_5pct_threshold <- function(model, data, marker, coma = "0", target = 0.05) {
  curve <- pooled_curve(model, data, marker, coma = coma, n = 5000)
  roots <- find_crossings(curve$x, curve$risk, target = target)

  if (length(roots) > 0) {
    threshold <- roots[1]
  } else {
    idx <- which.min(abs(curve$risk - target))
    threshold <- curve$x[idx]
  }

  data.frame(marker = marker, label = marker_specs$label[marker_specs$marker == marker],
             Coma = coma, threshold_5pct = threshold)
}

# Build the pooled threshold table for coma = 0
threshold_summary <- bind_rows(lapply(marker_specs$marker, function(m) {
  find_5pct_threshold(models[[m]], Data, m, coma = "0")})) %>%
  select(marker, label, Coma, threshold_5pct) %>%
  arrange(marker)

# Keep the same object name for a simple handoff to the second file
summary_coma0 <- threshold_summary

# Save the tables
write.csv(summary_coma0, file.path(OUT_DIR, "Thresholds_coma0.csv"), row.names = FALSE)

# Save the dataset, model list, and results for later use
saveRDS(list(data = Data, marker_specs = marker_specs, models = models,
             threshold_summary = threshold_summary, summary_coma0 = summary_coma0),
        file.path(OUT_DIR, "Thresholds_for_mortality.rds"))

summary_coma0