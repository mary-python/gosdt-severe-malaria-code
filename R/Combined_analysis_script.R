# This script loads a dataset and deep GOSDT tree (in the form of a table).
# It sends the data through the tree one patient at a time,
# and then follows a simple classification algorithm to prune the tree.
# A diagnostic performance table and RDS bundle are exported.
# The GOSDT tree plotting script uses this bundle to visualise the original and pruned trees.

library(data.table)
options(scipen = 999)

# Set up config (user preferences)
CFG <- list(dataset_path = "data/Combined_dataset copy.csv",
            splits_path = "Python/Combined_dataset_splits.csv",
            dp_csv = "outputs/Combined_diagnostic_performance.csv",
            result_rds = "R/Combined_tree_analysis_result.rds",
            thresh = 0.05, min_n = 5)

# ---- Function 1: Load the data and tree splits ----
load_data <- function(cfg) {
  dataset <- fread(cfg$dataset_path)
  split_info <- read.csv(cfg$splits_path)
  list(dataset = dataset, split_info = split_info)
}

# ---- Function 2: Separate the tree objects into nodes and leaves ----
separate_tree <- function(split_info) {
  list(node_df = subset(split_info, NodeType == "node"),
       leaf_df = subset(split_info, NodeType == "leaf"))
}

# ---- Function 3: Prepare the analysis data ----
prepare_data <- function(dataset, node_df, cfg) {

  # Remove patients with gender not M/F
  dataset <- dataset[SEX == "M" | SEX == "F"]

  # Remove patients with unidentified outcome
  dataset <- dataset[!is.na(Died_28D)]

  # Remove patients older than 15 years
  dataset <- dataset[Age_category != ">15"]

  # Remove studies with fewer than 20 patients
  dataset[, study_n := .N, by = STUDY]
  dataset <- dataset[study_n >= 20]
  dataset[, study_n := NULL]

  # Change respiratory distress from TRUE/FALSE to 1/0
  dataset$`Respiratory distress` <- as.integer(dataset$`Respiratory distress`)

  # This is the complete-case cohort used for pruning:
  # keep patients complete on all predictor columns in the original dataset
  cc_cohort <- na.omit(dataset[, c(3, 12, 13, 14, 15, 17, 21, 44, 66, 67, 69)])

  list(dataset = dataset, cc_cohort = cc_cohort)
}

# ---- Function 4: Create comparison function ----
create_comparison <- function() {
  function(dataset, predictor, symbol, value, lr) {

    if (!predictor %in% names(dataset)) {
      stop(paste0("Predictor not found in dataset: ", predictor))
    }

    if (!lr %in% c("L", "R")) {
      stop(paste0("Branch must be 'L' or 'R', not: ", lr))
    }

    x <- dataset[[predictor]]
    value_num <- suppressWarnings(as.numeric(value))

    if (is.na(value_num)) {
      stop(paste0("Non-numeric split value for ", predictor, ": ", value))
    }

    if (symbol == ">") {
      return(if (lr == "L") x > value_num else x <= value_num)
    } else if (symbol == "<") {
      return(if (lr == "L") x < value_num else x >= value_num)
    } else if (symbol == "<=") {
      return(if (lr == "L") x <= value_num else x > value_num)
    } else if (symbol == ">=") {
      return(if (lr == "L") x >= value_num else x < value_num)
    } else if (symbol %in% c("=", "==")) {
      return(if (lr == "L") x == value_num else x != value_num)
    }

    stop(paste0("Unrecognised comparison symbol: ", symbol))
  }
}

# ---- Function 5: Send patients down the tree ----
send_patients_down_tree <- function(node_df, leaf_df, comparison_function) {

  # Apply function to one patient (stored as a data frame)
  function(patient_df) {

    # Create an empty path (at root of the tree)
    current <- ""

    # Safety limit to avoid infinite loop
    for (k in 1:100) {

      # Look up the current node
      node <- node_df[node_df$Path == current, , drop = FALSE]
      if (nrow(node) == 0) break

      # Decide if the path should go left or right using the comparison function
      go_left <- comparison_function(patient_df, node$Predictor[1], node$Symbol[1],
                                     node$Value[1], "L")
      if (is.na(go_left[1])) break

      # Add an "L" or "R" to the path
      current <- paste0(current, if (go_left[1]) "L" else "R")

      # Once a leaf is found, return the path of the patient (a series of "L"s and "R"s)
      if (current %in% leaf_df$Path) return(current)
    }

    # If the loop finishes without successfully finding a leaf, return NA
    NA_character_
  }
}

# ---- Function 6: Assign patients to leaves in the tree ----
assign_patients_to_leaves <- function(my_subset, send_patients) {

  # Send each patient down the tree and record their path to a leaf (links to Function 5)
  path_to_leaf <- vapply(seq_len(nrow(my_subset)),
                      function(i) send_patients(my_subset[i, , drop = FALSE]), character(1))

  # Keep only the patients that successfully reach a leaf
  valid_idx <- !is.na(path_to_leaf)

  # Group patients by the leaf they land in
  leaf_members <- split(which(valid_idx), path_to_leaf[valid_idx])

  # Associate each leaf with the probability of death of its patients
  leaf_probs <- as.numeric(sapply(leaf_members, function(idx) mean(my_subset$Died_28D[idx])))
  names(leaf_probs) <- names(leaf_members)

  # Return the patients in each leaf and their probability of death
  list(path_to_leaf = path_to_leaf, leaf_members = leaf_members, leaf_probs = leaf_probs)
}

# ---- Function 7: Prune the tree ----
# Starting at the deepest node, prune it if:
# - both terminal leaves are above the mortality-risk threshold
# - both terminal leaves are below the mortality-risk threshold
# - at least one terminal leaf contains fewer than min_n patients
# Otherwise, keep the node.
prune_tree <- function(my_subset, leaf_probs, leaf_members, thresh, min_n) {
  repeat {

    # Extract all paths (consisting of "L" and "R")
    paths <- names(leaf_probs)

    # Look at the nodes one level above each of the leaves
    parents <- unique(substr(paths, 1, nchar(paths) - 1))
    parents <- parents[parents != ""]

    # Examine the nodes in reverse order (deepest first)
    parents <- parents[order(nchar(parents), decreasing = TRUE)]

    changed <- FALSE

    for (p in parents) {

      # Identify the left and right leaves of the current node
      left <- paste0(p, "L")
      right <- paste0(p, "R")

      left_exists <- left %in% names(leaf_probs)
      right_exists <- right %in% names(leaf_probs)

      # Only merge if BOTH children are current terminal leaves
      if (!(left_exists && right_exists)) {
        next
      }

      # Fetch the number of patients and probability of death in both left and right leaves
      n1 <- length(leaf_members[[left]])
      n2 <- length(leaf_members[[right]])
      s1 <- leaf_probs[[left]]
      s2 <- leaf_probs[[right]]

      should_merge <- FALSE

        # If either leaf has fewer than min_n patients assigned to it, remove the node
        if (n1 < min_n || n2 < min_n) {
          should_merge <- TRUE

        # If both leaves are on the same side of mortality-risk threshold, remove the node
        } else if ((s1 >= thresh && s2 >= thresh) ||
                   (s1 < thresh && s2 < thresh)) {
          should_merge <- TRUE
        }

        # If the node is removed
        if (should_merge) {

          # Combine the patients from both leaves
          new_members <- c(leaf_members[[left]], leaf_members[[right]])

          # Put the patients and combined probability of death in place of the node
          leaf_probs[[p]] <- as.numeric(mean(my_subset$Died_28D[new_members]))
          leaf_members[[p]] <- new_members

          # Remove the leaves from their old positions
          leaf_probs <- leaf_probs[!names(leaf_probs) %in% c(left, right)]
          leaf_members <- leaf_members[!names(leaf_members) %in% c(left, right)]

          changed <- TRUE
        }
      }

    # If there are no more changes, stop pruning
    if (!changed) break
  }

  # Return list of patients and probability of death assigned to each leaf
  list(probs = leaf_probs, members = leaf_members)
}

# ---- Function 8: Return all node paths in the pruned tree ----
return_node_paths <- function(paths) {

  # Always include the root
  out <- c("")

  # For each path to a leaf
  for (p in paths) {

    # Skip invalid or trivial paths
    if (is.na(p) || !nzchar(p) || nchar(p) <= 1) next

    # Remove the final character from each valid path
    out <- c(out, vapply(1:(nchar(p) - 1), function(k) substr(p, 1, k), character(1)))
  }

  # Return all node paths without duplicates
  unique(out)
}

# ---- Function 9: Send patients down the pruned tree ----
send_patients_down_pruned_tree <- function(node_df, pruned_paths_to_leaf, pruned_node_paths,
                                           comparison_function) {

  # Store the nodes that exist after pruning
  pruned_node_df <- node_df[node_df$Path %in% pruned_node_paths, , drop = FALSE]

  # Apply function to one patient (stored as a data frame)
  function(patient_df) {

    # Start at the root of the tree
    current <- ""

    # Safety limit to avoid infinite loop
    for (k in 1:100) {

      # Look up the current node
      node <- pruned_node_df[pruned_node_df$Path == current, , drop = FALSE]
      if (nrow(node) == 0) break

      # Decide if the path should go left or right using the comparison function
      go_left <- comparison_function(patient_df, node$Predictor[1], node$Symbol[1],
                                     node$Value[1], "L")
      if (is.na(go_left[1])) return(NA_character_)

      # Add this left/right step to the path
      next_path <- paste0(current, if (go_left[1]) "L" else "R")

      # If the current path matches a valid path to a leaf, return it
      if (next_path %in% pruned_paths_to_leaf) return(next_path)

      # If the current path is not present in the tree, then return NA
      if (!(next_path %in% pruned_node_paths)) return(NA_character_)

      # Otherwise, keep traversing
      current <- next_path
    }

    # If the loop finishes without successfully finding a leaf, return NA
    NA_character_
  }
}

# ---- Function 10: Classify patients using original/pruned tree ----
# Either: 1 = high-risk / died
#     Or: 0 = low-risk / survived
# Prediction rule: 1 if leaf death probability >= 0.05, else 0
classify_with_tree <- function(dat, leaf_probs_used, route_function, thresh = 0.05) {

  # Loop through each patient row
  paths <- vapply(seq_len(nrow(dat)),

                  # Determine which leaf the patient lands in
                  function(i) route_function(dat[i, , drop = FALSE]), character(1))

  # Remove invalid rows
  ok <- !is.na(paths) & !is.na(dat$Died_28D)

  # Return an empty table if nothing can be evaluated
  if (!any(ok)) {
    return(data.table(path = character(0), score = numeric(0), pred = integer(0),
                      truth = integer(0)))
  }

  # Get the probability of death for the patients in current leaf
  scores <- as.numeric(unlist(leaf_probs_used[paths[ok]]))

  # Convert score (0.00-1.00) to prediction (0 or 1)
  # Set prediction to 0 if score less than thresh
  # Set prediction to 1 if score greater or equal to thresh
  pred <- as.integer(scores >= thresh)
  truth <- as.integer(dat$Died_28D[ok])

  # Return path, score, and prediction
  data.table(path = paths[ok], score = scores, pred = pred, truth = truth)
}

# ---- Function 11: Find the Patient ID column if one exists ----
find_id_col <- function(dat) {

  # Store common names used for patient ID columns
  candidates <- c("PATIENT_ID", "Patient_ID", "PatientID", "ID", "subject_id", "SUBJECT_ID")

  # Identify candidate names present in the dataset
  hit <- candidates[candidates %in% names(dat)]

  # Return the first matching ID column, or NULL if none are found
  if (length(hit) == 0) NULL else hit[1]
}

# ---- Function 12: Create tables to double-check diagnostic performance of trees ----
classify_check <- function(dat, leaf_probs_used, route_function, thresh = 0.05,
                                      id_col = NULL) {

  # Loop through each patient row
  paths <- vapply(seq_len(nrow(dat)),

                  # Determine which leaf the patient lands in
                  function(i) route_function(dat[i, , drop = FALSE]), character(1))

  # Remove invalid rows
  ok <- !is.na(paths) & !is.na(dat$Died_28D)

  # Return an empty table if nothing can be evaluated
  if (!any(ok)) {
    return(data.table(row_in_subset = integer(0), patient_id = character(0),
                      path = character(0), leaf_death_prob = numeric(0),
                      leaf_pred_class = integer(0), truth = integer(0),
                      correct = integer(0)))
  }

  # Store row indices of valid patients
  idx <- which(ok)

  # Get the probability of death for the patients in current leaf
  scores <- as.numeric(unlist(leaf_probs_used[paths[ok]]))

  # Convert score (0.00-1.00) to prediction (0 or 1)
  # Set prediction to 0 if score less than thresh
  # Set prediction to 1 if score greater than thresh
  pred <- as.integer(scores >= thresh)
  truth <- as.integer(dat$Died_28D[ok])

  # Fetch Patient ID from relevant column if available
  # If not available then use row number
  patient_id <- if (!is.null(id_col) && id_col %in% names(dat)) {
    as.character(dat[[id_col]][ok])
  } else {
    as.character(idx)
  }

  # Return patient-level table for checking diagnostic performance
  data.table(row_in_subset = idx, patient_id = patient_id, path = paths[ok],
             leaf_death_prob = scores, leaf_pred_class = pred, truth = truth,
             correct = as.integer(pred == truth))
}

# ---- Function 13: Build a summary table of diagnostic performance ----
build_summary <- function(pred_df, tree, cohort) {

  # Store number of rows in current classification data frame
  n <- nrow(pred_df)

  # Compute true positives, false positives, true negatives, and false negatives
  tp <- sum(pred_df$pred == 1 & pred_df$truth == 1, na.rm = TRUE)
  fp <- sum(pred_df$pred == 1 & pred_df$truth == 0, na.rm = TRUE)
  tn <- sum(pred_df$pred == 0 & pred_df$truth == 0, na.rm = TRUE)
  fn <- sum(pred_df$pred == 0 & pred_df$truth == 1, na.rm = TRUE)

  # Compute sensitivity, specificity, and accuracy
  sensitivity <- if ((tp + fn) > 0) 100 * tp / (tp + fn) else NA_real_
  specificity <- if ((tn + fp) > 0) 100 * tn / (tn + fp) else NA_real_
  accuracy <- if (n > 0) 100 * (tp + tn) / n else NA_real_

  # Return diagnostic performance table
  data.table(tree = tree, cohort = cohort, n = n, TP = tp, FP = fp, TN = tn, FN = fn,
             sensitivity = round(sensitivity, 1), specificity = round(specificity, 1),
             accuracy = round(accuracy, 1))
}

# ---- Function 14: Classify patients using severe malaria criteria ----
classify_with_criteria <- function(dat) {

  # Convert indicator variables to TRUE/FALSE values
  present <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    !is.na(x) & x > 0
  }

  # Apply each severe malaria criterion and record whether it is present
  coma <- if ("Coma" %in% names(dat)) {
    present(dat$Coma)
  } else rep(FALSE, nrow(dat))

  resp <- if ("Respiratory distress" %in% names(dat)) {
    present(dat$`Respiratory distress`)
  } else rep(FALSE, nrow(dat))

  bun <- if ("BUN" %in% names(dat)) {
    suppressWarnings(as.numeric(dat$BUN)) > 20
  } else rep(FALSE, nrow(dat))

  bicarb <- if ("Bicarbonate_mEq_L" %in% names(dat)) {
    suppressWarnings(as.numeric(dat$Bicarbonate_mEq_L)) < 15
  } else rep(FALSE, nrow(dat))

  lactate <- if ("Lactate" %in% names(dat)) {
    suppressWarnings(as.numeric(dat$Lactate)) > 5
  } else rep(FALSE, nrow(dat))

  paras <- if ("para_ul" %in% names(dat)) {
    suppressWarnings(as.numeric(dat$para_ul)) > 100000
  } else rep(FALSE, nrow(dat))

  glucose <- if ("Glucose_mmol_L" %in% names(dat)) {
    suppressWarnings(as.numeric(dat$Glucose_mmol_L)) < 2.2
  } else rep(FALSE, nrow(dat))

  hb <- if ("Hb" %in% names(dat)) {
    suppressWarnings(as.numeric(dat$Hb)) < 5
  } else rep(FALSE, nrow(dat))

  # Classify patients as having severe malaria if any criterion is present
  pred <- as.integer(coma | resp | bun | bicarb | lactate | paras | glucose | hb)
  truth <- as.integer(dat$Died_28D)
  ok <- !is.na(truth)

  # Use patient_id if available, otherwise use row number
  pid <- if ("patient_id" %in% names(dat)) dat$patient_id[ok] else which(ok)

  # Return patient ID, prediction, and true outcome
  data.table(patient_id = pid, pred = pred[ok], truth = truth[ok])
}

# ---- Function 15: Build a one-row diagnostic performance table ----
build_comparison_summary <- function(pred_df) {

  # Compute true positives, false positives, true negatives, and false negatives
  n <- nrow(pred_df)
  tp <- sum(pred_df$pred == 1 & pred_df$truth == 1, na.rm = TRUE)
  fp <- sum(pred_df$pred == 1 & pred_df$truth == 0, na.rm = TRUE)
  tn <- sum(pred_df$pred == 0 & pred_df$truth == 0, na.rm = TRUE)
  fn <- sum(pred_df$pred == 0 & pred_df$truth == 1, na.rm = TRUE)

  # Compute sensitivity, specificity, and accuracy
  sensitivity <- if ((tp + fn) > 0) 100 * tp / (tp + fn) else NA_real_
  specificity <- if ((tn + fp) > 0) 100 * tn / (tn + fp) else NA_real_
  accuracy <- if (n > 0) 100 * (tp + tn) / n else NA_real_

  # Return a one-row diagnostic performance table
  data.table(n = n, TP = tp, FP = fp, TN = tn, FN = fn,
             sensitivity = round(sensitivity, 1),
             specificity = round(specificity, 1),
             accuracy = round(accuracy, 1))
}

# ======================
# RUN THE FULL ANALYSIS
# ======================
run_analysis <- function(cfg) {

  # Call Function 1 (load data)
  inputs <- load_data(cfg)

  # Call Function 2 (separate tree objects)
  tree_frames <- separate_tree(inputs$split_info)

  # Call Function 3 (prepare analysis data)
  prep <- prepare_data(inputs$dataset, tree_frames$node_df, cfg)

  # Prepare objects for next functions
  node_df <- tree_frames$node_df
  leaf_df <- tree_frames$leaf_df
  dataset <- prep$dataset
  cc_cohort <- prep$cc_cohort

  # Call Function 4 (create comparison function)
  comparison_function <- create_comparison()

  # Call Function 5 (send patients down tree)
  send_patients <- send_patients_down_tree(node_df, leaf_df, comparison_function)

  # Call Function 6 (assign patients to leaves)
  original <- assign_patients_to_leaves(cc_cohort, send_patients)
  orig_leaf_probs <- original$leaf_probs
  orig_leaf_members <- original$leaf_members

  # Call Function 7 (prune tree)
  pruned_tree <- prune_tree(my_subset = cc_cohort, leaf_probs = original$leaf_probs,
                            leaf_members = original$leaf_members, thresh = cfg$thresh,
                            min_n = cfg$min_n)

  # Extract paths from pruned tree that end in a leaf
  pruned_paths_to_leaf <- names(pruned_tree$probs)

  # Call Function 8 (return node paths)
  pruned_node_paths <- return_node_paths(pruned_paths_to_leaf)

  # Call Function 9 (send patients down pruned tree)
  send_patients_pruned <- send_patients_down_pruned_tree(node_df, pruned_paths_to_leaf,
                                                         pruned_node_paths,
                                                         comparison_function)

  # Call Function 10 (classify complete-case cohort using original tree)
  pred_cc_orig <- classify_with_tree(cc_cohort, orig_leaf_probs,
                                     route_function = send_patients, thresh = cfg$thresh)

  # Call Function 10 (classify complete-case cohort using pruned tree)
  pred_cc_pruned <- classify_with_tree(cc_cohort, pruned_tree$probs,
                                       route_function = send_patients_pruned,
                                       thresh = cfg$thresh)

  # Call Function 11 (find patient ID column)
  id_col <- find_id_col(dataset)

  # Call Function 12 (create table to double-check classification of complete-case cohort
  #                   using original tree)
  check_tables <- list(original_cc_cohort = classify_check(cc_cohort, orig_leaf_probs,
                                                           send_patients,
                                                           thresh = cfg$thresh,
                                                           id_col = id_col),

                       # Call Function 12 (create table to double-check classification of
                       #                   complete-case cohort using pruned tree)
                       pruned_cc_cohort = classify_check(cc_cohort, pruned_tree$probs,
                                                         send_patients_pruned,
                                                         thresh = cfg$thresh,
                                                         id_col = id_col))

  # Create output folder
  check_dir <- file.path("outputs", "Combined_performance_check")
  dir.create(check_dir, recursive = TRUE, showWarnings = FALSE)

  # Export one CSV per patient-level table
  for (nm in names(check_tables)) {
    fwrite(check_tables[[nm]], file = file.path(check_dir, paste0(nm, ".csv")))
  }

  # Call Function 13 (build diagnostic performance table)
  diagnostic_performance <- rbindlist(list(build_summary(pred_cc_orig, "original",
                                                         "complete-case cohort"),
                                           build_summary(pred_cc_pruned, "pruned",
                                                         "complete-case cohort")))

  # Export diagnostic performance table
  write.csv(diagnostic_performance, file = cfg$dp_csv)

  # Create a complete-case dataset containing the severe malaria criteria variables
  cols_criteria <- unique(c("patient_id", "Coma", "Respiratory distress", "BUN",
                            "Bicarbonate_mEq_L", "Lactate", "para_ul",
                            "Glucose_mmol_L", "Hb", "Died_28D"))
  cols_criteria <- cols_criteria[cols_criteria %in% names(dataset)]
  criteria_subset <- na.omit(dataset[, ..cols_criteria])

  # Call Function 14 (classify patients using severe malaria criteria)
  pred_criteria <- classify_with_criteria(criteria_subset)

  # Call Function 15 (build a one-row diagnostic performance table)
  comparison_summary <- build_comparison_summary(pred_criteria)

  fwrite(comparison_summary,
         file = file.path("outputs", "Combined_diagnostic_perf_comparison.csv"))

  # Export an RDS bundle containing objects needed by the plotting script
  saveRDS(list(node_df = node_df, leaf_df = leaf_df, orig_leaf_probs = orig_leaf_probs,
               orig_leaf_members = orig_leaf_members, pruned_tree = pruned_tree,
               thresh = cfg$thresh, pruned_paths_to_leaf = pruned_paths_to_leaf,
               pruned_node_paths = pruned_node_paths), file = cfg$result_rds)

  # Return analysis objects without printing them to the console
  invisible(list(diagnostic_performance = diagnostic_performance,
                 check_tables = check_tables,
                 pruned_tree = pruned_tree, cc_cohort = cc_cohort,
                 original_leaf_probs = orig_leaf_probs,
                 original_leaf_members = orig_leaf_members))
}

# Generate result
result <- run_analysis(CFG)