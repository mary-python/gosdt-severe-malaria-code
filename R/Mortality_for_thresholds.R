library(dplyr)
library(tidyr)
library(ggplot2)
library(mgcv)

# Set the project folder and output folder
OUT_DIR <- "data"

# Load the saved models and pooled thresholds from the first script
obj <- readRDS(file.path(OUT_DIR, "Thresholds_for_mortality.rds"))

# Extract the data, marker information, fitted models, and pooled thresholds
Data <- obj$data
marker_specs <- obj$marker_specs
models <- obj$models

# List the thresholds that we want to evaluate
threshold_grid <- tibble::tribble(~marker, ~label, ~value,
                                  "Bicarbonate_mEq_L", "Bicarbonate", 10,
                                  "Bicarbonate_mEq_L", "Bicarbonate", 12,
                                  "Bicarbonate_mEq_L", "Bicarbonate", 15,
                                  "Temperature_C", "Temperature", 35,
                                  "Temperature_C", "Temperature", 37.5,
                                  "Temperature_C", "Temperature", 40,
                                  "Lactate", "Lactate", 5,
                                  "Lactate", "Lactate", 7,
                                  "Lactate", "Lactate", 10,
                                  "Glucose_mmol_L", "Glucose", 2.2,
                                  "Glucose_mmol_L", "Glucose", 3.0,
                                  "Glucose_mmol_L", "Glucose", 3.5,
                                  "Glucose_mmol_L", "Glucose", 4.0,
                                  "Hb", "Hb", 4.0,
                                  "Hb", "Hb", 5.0,
                                  "Hb", "Hb", 6.0)

# Get the study and site levels used in the fitted model
get_ref_levels <- function(model, data) {
  mf <- tryCatch(model.frame(model), error = function(e) NULL)

  study_ref <- if (!is.null(mf) && "STUDY" %in% names(mf) && !is.null(mf$STUDY)) {
    levels(droplevels(mf$STUDY))[1]
  } else {
    levels(data$STUDY)[1]
  }

  site_ref <- if (!is.null(mf) && "SITEID" %in% names(mf) && !is.null(mf$SITEID)) {
    levels(droplevels(mf$SITEID))[1]
  } else {
    levels(data$SITEID)[1]
  }

  list(study_ref = study_ref, site_ref = site_ref)
}

# Create a small dataset for making predictions
build_newdata <- function(data, marker, value, age_category, coma, study_ref, site_ref) {
  nd <- data.frame(Age_category = factor(age_category, levels = levels(data$Age_category)),
                   Coma = factor(coma, levels = levels(data$Coma)),
                   STUDY = factor(study_ref, levels = levels(data$STUDY)),
                   SITEID = factor(site_ref, levels = levels(data$SITEID)))

  nd[[marker]] <- value
  nd
}

# Calculate the predicted mortality risk for one threshold
risk_at_threshold <- function(model, data, marker, value, age_category, coma) {
  refs <- get_ref_levels(model, data)
  nd <- build_newdata(data, marker, value, age_category, coma, refs$study_ref, refs$site_ref)
  as.numeric(predict(model, newdata = nd, type = "response"))
}

# Calculate the pooled mortality risk across age groups for coma = 0
pooled_risk_at_threshold <- function(model, data, marker, value, coma = "0") {

  age_levels <- levels(data$Age_category)

  w <- prop.table(table(data$Age_category[data$Coma == coma]))
  w <- w[age_levels]
  w[is.na(w)] <- 0

  risks <- sapply(age_levels, function(a) {
    risk_at_threshold(model, data, marker, value, a, coma)
  })

  as.numeric(sum(risks * w) / sum(w))
}

# Calculate pooled mortality risk for every threshold using coma = 0
result_table <- bind_rows(lapply(seq_len(nrow(threshold_grid)), function(i) {
  m <- threshold_grid$marker[i]
  v <- threshold_grid$value[i]
  lab <- threshold_grid$label[i]

  data.frame(marker = m,
             label = lab,
             threshold_value = v,
             Coma = 0,
             mortality_risk_pct =
               round(pooled_risk_at_threshold(models[[m]], Data, m, v, coma = "0") * 100, 1))
})) %>%
  arrange(label, threshold_value)

# Save the results table
write.csv(result_table, file.path(OUT_DIR, "Mortality_for_thresholds.csv"), row.names = FALSE)

# Display the final table
result_table