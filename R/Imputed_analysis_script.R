# This script loads a file containing 20 imputations of a dataset,
# and a deep GOSDT tree (in the form of a table) for each imputation.
# It sends each imputation through its tree one patient at a time,
# and then follows a simple classification algorithm to prune the tree.
# A diagnostic performance table and RDS bundle are exported for each imputation.
# The GOSDT tree plotting script uses each bundle to visualise the original and pruned trees
# for each imputation.

library(data.table)
options(scipen = 999)

# State paths of input files
IMPUTED_LONG_PATH <- "data/Imputed_dataset_long.csv"
SPLITS_DIR <- file.path("Python", "Imputed_dataset_splits")

# State paths of output folders
OUT_DIR <- "outputs"
DIAG_DIR <- file.path(OUT_DIR, "Imputed_diagnostic_performance")
RDS_DIR <- file.path("R", "Imputed_tree_analysis_result")

# If any output folder does not exist then create it
for (d in c(OUT_DIR, DIAG_DIR, RDS_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# State analysis settings
CFG <- list(thresh = 0.05, min_n = 424)

# ---- Function 1: Load the imputed dataset ----
# The CSV is expected to contain the 20 imputations stacked one after another
load_imputed_long_dataset <- function(long_path, n_imputations = 20L) {
  mi_dataset <- fread(file = long_path)

  # Check whether the dataset already has an explicit imputation_id column
  if ("imputation_id" %in% names(mi_dataset)) {
    mi_dataset[, imputation_id := as.integer(imputation_id)]
    setorder(mi_dataset, imputation_id)

    # Ensure each imputation has a patient_id column
    # Create one if needed
    if (!"patient_id" %in% names(mi_dataset)) {
      mi_dataset[, patient_id := seq_len(.N), by = imputation_id]
    } else {
      mi_dataset[, patient_id := as.integer(patient_id)]
    }

  # If no imputation_id column exists, infer the imputation structure from row order
  } else {
    n_rows <- nrow(mi_dataset)

    # Check that the number of rows is compatible with 20 imputations
    if (n_rows %% n_imputations != 0) {
      stop(sprintf(paste0("Expected the imputed file to contain a multiple of %d rows,",
                          " but found %d."), n_imputations, n_rows))
    }

    # Compute the number of patients per imputation
    n_patients <- as.integer(n_rows / n_imputations)

    # Create imputation_id and patient_id columns
    mi_dataset[, imputation_id := rep(seq_len(n_imputations), each = n_patients)]
    mi_dataset[, patient_id := rep(seq_len(n_patients), times = n_imputations)]
  }

  # Split the stacked dataset into a list of 20 imputation-specific data.tables
  split(mi_dataset, by = "imputation_id", keep.by = FALSE)
}

# ---- Function 2: Separate the tree object into nodes and leaves ----
separate_tree <- function(split_info) {
  list(node_df = subset(split_info, NodeType == "node"),
       leaf_df = subset(split_info, NodeType == "leaf"))
}

# ---- Function 3: Prepare the imputed dataset ----
prepare_imputed_dataset <- function(dataset) {
  dataset <- copy(dataset)

  # Convert Respiratory distress from logical/character form to integer if present
  if ("Respiratory distress" %in% names(dataset)) {
    dataset[, `Respiratory distress` := as.integer(`Respiratory distress`)]
  }

  # Convert Died_28D to integer if present and verify that no values fail conversion
  if ("Died_28D" %in% names(dataset)) {
    raw_outcome <- trimws(as.character(dataset$Died_28D))
    dataset[, Died_28D := as.integer(raw_outcome)]

    # Stop if any outcome values could not be converted to integer
    if (anyNA(dataset$Died_28D)) {
      bad_vals <- unique(raw_outcome[is.na(dataset$Died_28D)])
      stop(sprintf("Died_28D contains values that could not be converted to integer: %s",
                   paste(bad_vals, collapse = ", ")))
    }
  }

  dataset
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
assign_patients_to_leaves <- function(my_subset, send_patients, leaf_df) {

  # Send each patient down the tree and record their path to a leaf (links to Function 5)
  path_to_leaf <- vapply(seq_len(nrow(my_subset)),
                         function(i) send_patients(my_subset[i, , drop = FALSE]),
                         character(1))

  # Store all leaf paths present in the tree
  all_leaf_paths <- as.character(leaf_df$Path)

  # Create an empty patient-membership list for every leaf
  leaf_members <- setNames(vector("list", length(all_leaf_paths)), all_leaf_paths)

  # Create an empty probability-of-death vector for every leaf
  leaf_probs <- setNames(rep(NA_real_, length(all_leaf_paths)), all_leaf_paths)

  # Keep only the patients that successfully reach a leaf
  valid_idx <- !is.na(path_to_leaf)

  # Identify routed patients and group them by leaf path
  routed_idx <- which(valid_idx)
  routed_paths <- path_to_leaf[valid_idx]
  split_members <- split(routed_idx, routed_paths)

  # Fill in patient memberships and probabilities for leaves that contain patients
  for (p in names(split_members)) {
    leaf_members[[p]] <- split_members[[p]]
    leaf_probs[[p]] <- mean(my_subset$Died_28D[split_members[[p]]], na.rm = TRUE)
  }

  # Return the patients in each leaf and their probability of death
  list(path_to_leaf = path_to_leaf, leaf_members = leaf_members, leaf_probs = leaf_probs)
}

# ---- Function 7: Prune the tree ----
# Starting at the deepest node, prune it if:
# - both terminal leaves are above the mortality-risk threshold
# - both terminal leaves are below txhe mortality-risk threshold
# - at least one terminal leaf contains fewer than min_n patients
# Otherwise, keep the node.
prune_tree <- function(my_subset, leaf_probs, leaf_members, thresh, min_n) {

  # Function 7a: Safely compute the probability of death of the patients in the input
  compute_prob_death <- function(idx) {
    if (length(idx) == 0) return(NA_real_)
    mean(my_subset$Died_28D[idx], na.rm = TRUE)
  }

  # Store the current leaves and their patient memberships
  current_probs <- leaf_probs
  current_members <- leaf_members

  # Function 7b: Remove any remaining terminal leaves with fewer than min_n patients
  remove_small_leaves <- function() {

    leaf_sizes <- lengths(current_members)
    small_leaf_paths <- names(leaf_sizes)[leaf_sizes < min_n]

    if (length(small_leaf_paths) == 0) return(FALSE)

    current_probs <<- current_probs[!names(current_probs) %in% small_leaf_paths]
    current_members <<- current_members[!names(current_members) %in% small_leaf_paths]

    TRUE
  }

  repeat {

    # Extract all current terminal leaf paths
    paths <- names(current_probs)
    if (length(paths) == 0) break

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

      left_is_leaf <- left %in% names(current_probs)
      right_is_leaf <- right %in% names(current_probs)

      # Only merge if BOTH children are current terminal leaves
      if (!(left_is_leaf && right_is_leaf)) {
        next
      }

      # Fetch the number of patients and probability of death in both left and right leaves
      n1 <- length(current_members[[left]])
      n2 <- length(current_members[[right]])
      s1 <- current_probs[[left]]
      s2 <- current_probs[[right]]

      should_merge <- FALSE

      # If either leaf has fewer than min_n patients assigned to it, remove the node
      if (n1 < min_n || n2 < min_n) {
        should_merge <- TRUE

      # If both leaves are on the same side of mortality-risk threshold, remove the node
      } else if (!is.na(s1) && !is.na(s2) &&
                 ((s1 >= thresh && s2 >= thresh) ||
                  (s1 < thresh && s2 < thresh))) {
        should_merge <- TRUE
      }

      # If the node is removed
      if (should_merge) {

        # Combine the patients from both leaves
        new_members <- c(current_members[[left]], current_members[[right]])

        # Call Function 7a (compute probability of death of patients)
        current_probs[[p]] <- compute_prob_death(new_members)

        # Update the patient membership of the current leaf
        current_members[[p]] <- new_members

        # Remove the leaves from their old positions
        current_probs <- current_probs[!names(current_probs) %in% c(left, right)]
        current_members <- current_members[!names(current_members) %in% c(left, right)]

        changed <- TRUE
      }
    }

    # Call Function 7b (look for terminal leaves with fewer than min_n patients)
    if (remove_small_leaves()) {
      changed <- TRUE
    }

    # If there are no more changes, stop pruning
    if (!changed) break
  }

  # Return list of patients and probability of death assigned to each leaf
  list(probs = current_probs, members = current_members)
}

# ---- Function 8: Compact the pruned tree structure ----
# Remove decision nodes that have only one surviving child, and move the
# surviving branch up to the nearest remaining decision point.
compact_pruned_tree <- function(node_df, leaf_members, min_n) {

  # Standardise the input node table
  node_df <- as.data.table(node_df)
  node_df[, Path := as.character(Path)]
  node_df[, Predictor := as.character(Predictor)]
  node_df[, Symbol := as.character(Symbol)]
  node_df[, Value := as.character(Value)]

  # Keep only surviving leaves with valid names
  keep <- !is.na(names(leaf_members))
  leaf_members <- leaf_members[keep]

  # Keep only leaves with at least min_n patients
  leaf_sizes <- lengths(leaf_members)
  leaf_paths <- names(leaf_sizes)[leaf_sizes >= min_n]

  # Store the compacted node table and leaf information as they are built
  node_rows <- list()
  compact_leaf_paths <- character(0)

  # Check whether any surviving leaf lies below the requested prefix
  has_surviving_leaf <- function(prefix) {
    any(startsWith(leaf_paths, prefix))
  }

  # Add one node row to the compacted table
  add_node <- function(path, row) {
    node_rows[[length(node_rows) + 1]] <<- data.table(
      Path = path,
      Predictor = row$Predictor[1],
      Symbol = row$Symbol[1],
      Value = row$Value[1]
    )
  }

  # Add one leaf row to the compacted table
  add_leaf <- function(path) {
    compact_leaf_paths <<- c(compact_leaf_paths, path)
  }

  # Walk the original tree and rebuild the compacted one
  recurse <- function(original_path, compact_path) {

    left_orig <- paste0(original_path, "L")
    right_orig <- paste0(original_path, "R")
    left_has <- has_surviving_leaf(left_orig)
    right_has <- has_surviving_leaf(right_orig)

    # Stop if no surviving leaf remains below this point
    if (!left_has && !right_has) {
      if (original_path %in% leaf_paths) {
        add_leaf(compact_path)
      }
      return(invisible(NULL))
    }

    # If both children remain, keep this node and follow both branches
    if (left_has && right_has) {
      row <- node_df[Path == original_path]
      if (nrow(row) > 0) {
        add_node(compact_path, row)
      }

      recurse(left_orig, paste0(compact_path, "L"))
      recurse(right_orig, paste0(compact_path, "R"))
      return(invisible(NULL))
    }

    # If only one side remains, remove the current decision and keep walking
    if (left_has) {
      recurse(left_orig, compact_path)
    } else if (right_has) {
      recurse(right_orig, compact_path)
    }

    invisible(NULL)
  }

  recurse("", "")

  # Build the compacted outputs
  pruned_node_df <- if (length(node_rows) > 0) {
    rbindlist(node_rows, fill = TRUE)
  } else {
    data.table(Path = character(0), Predictor = character(0),
               Symbol = character(0), Value = character(0))
  }

  if (nrow(pruned_node_df) > 0) {
    setorder(pruned_node_df, Path)
  }

  pruned_leaf_df <- data.table(Path = compact_leaf_paths)
  if (nrow(pruned_leaf_df) > 0) {
    setorder(pruned_leaf_df, Path)
  }

  list(node_df = pruned_node_df, leaf_df = pruned_leaf_df)
}

# ---- Function 9: Build a stable signature for the pruned tree ----
make_pruned_tree_signature <- function(node_df, pruned_node_paths) {

  # Keep only the nodes that remain after pruning and extract their split rules
  sig_df <- node_df[node_df$Path %in% pruned_node_paths,
                    .(Path, Predictor, Symbol, Value)]

  # Order nodes by path so identical trees produce identical signatures
  setorder(sig_df, Path)

  # Concatenate node path, predictor, operator, and threshold into a unique signature
  paste(apply(sig_df, 1, function(x) {paste(x["Path"], x["Predictor"], x["Symbol"],
                                            x["Value"], sep = ":")}), collapse = "|")
}
# ---- Function 10: Send patients down the pruned tree ----
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

# ---- Function 11: Classify patients using original/pruned tree ----
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
    return(data.table(patient_id = integer(0), path = character(0), score = numeric(0),
                      pred = integer(0), truth = integer(0)))
  }

  # Get the probability of death for the patients in current leaf
  scores <- as.numeric(leaf_probs_used[paths[ok]])

  # Convert score (0.00-1.00) to prediction (0 or 1)
  # Set prediction to 0 if score less than thresh
  # Set prediction to 1 if score greater or equal to thresh
  pred <- as.integer(scores >= thresh)
  truth <- as.integer(dat$Died_28D[ok])

  # Use patient_id if available, otherwise use row number
  pid <- if ("patient_id" %in% names(dat)) dat$patient_id[ok] else which(ok)

  # Return path, score, and prediction
  data.table(patient_id = pid, path = paths[ok], score = scores, pred = pred, truth = truth)
}

# ---- Function 12: Build a summary table of diagnostic performance ----
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

# ---- Function 13: Classify patients using severe malaria criteria ----
classify_with_criteria <- function(dat) {

  # Convert indicator variables to TRUE/FALSE values
  present <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    !is.na(x) & x > 0
  }

  # Apply each severe malaria criterion and record whether it is present
  if ("Coma" %in% names(dat)) {
    coma <- present(dat$Coma)
  } else {
    coma <- rep(FALSE, nrow(dat))
  }

  if ("Respiratory distress" %in% names(dat)) {
    resp <- present(dat$`Respiratory distress`)
  } else {
    resp <- rep(FALSE, nrow(dat))
  }

  if ("BUN" %in% names(dat)) {
    bun <- suppressWarnings(as.numeric(dat$BUN)) > 20
  } else {
    bun <- rep(FALSE, nrow(dat))
  }

  if ("Bicarbonate_mEq_L" %in% names(dat)) {
    bicarb <- suppressWarnings(as.numeric(dat$Bicarbonate_mEq_L)) < 15
  } else {
    bicarb <- rep(FALSE, nrow(dat))
  }

  if ("Lactate" %in% names(dat)) {
    lactate <- suppressWarnings(as.numeric(dat$Lactate)) > 5
  } else {
    lactate <- rep(FALSE, nrow(dat))
  }

  if ("para_ul" %in% names(dat)) {
    paras <- suppressWarnings(as.numeric(dat$para_ul)) > 100000
  } else {
    paras <- rep(FALSE, nrow(dat))
  }

  if ("Glucose_mmol_L" %in% names(dat)) {
    glucose <- suppressWarnings(as.numeric(dat$Glucose_mmol_L)) < 2.2
  } else {
    glucose <- rep(FALSE, nrow(dat))
  }

  if ("Hb" %in% names(dat)) {
    hb <- suppressWarnings(as.numeric(dat$Hb)) < 5
  } else {
    hb <- rep(FALSE, nrow(dat))
  }

  # Classify patients as having severe malaria if any criterion is present
  pred <- as.integer(coma | resp | bun | bicarb | lactate | paras | glucose | hb)
  truth <- as.integer(dat$Died_28D)
  ok <- !is.na(truth)

  # Use patient_id if available, otherwise use row number
  pid <- if ("patient_id" %in% names(dat)) dat$patient_id[ok] else which(ok)

  # Return patient ID, prediction, and true outcome
  data.table(patient_id = pid, pred = pred[ok], truth = truth[ok])
}

# ---- Function 14: Build an average diagnostic table for a group of trees ----
build_group_diagnostic_table <- function(diag_files) {

  # Read and combine all diagnostic performance tables in the group
  all_diag <- rbindlist(lapply(diag_files, fread), fill = TRUE)

  # Average diagnostic performance separately for original and pruned trees
  out <- all_diag[, .(n = round(mean(n, na.rm = TRUE)), TP = round(mean(TP, na.rm = TRUE)),
                      FP = round(mean(FP, na.rm = TRUE)), TN = round(mean(TN, na.rm = TRUE)),
                      FN = round(mean(FN, na.rm = TRUE)),
                      sensitivity = mean(sensitivity, na.rm = TRUE),
                      specificity = mean(specificity, na.rm = TRUE),
                      accuracy = mean(accuracy, na.rm = TRUE)), by = .(tree, cohort)]

  # Round performance measures to one decimal place
  out[, sensitivity := round(sensitivity, 1)]
  out[, specificity := round(specificity, 1)]
  out[, accuracy := round(accuracy, 1)]

  out[]
}

# ---- Function 15: Build an averaged pruned-tree bundle for one repeated group ----
build_group_pruned_bundle <- function(rds_files) {

  # Read all RDS bundles belonging to the current repeated tree group
  bundles <- lapply(rds_files, readRDS)

  # Use the first bundle as the structural template
  template <- bundles[[1]]

  # Extract the leaf and node paths from the template bundle
  leaf_paths <- as.character(template$pruned_leaf_df$Path)
  node_paths <- as.character(template$pruned_node_df$Path)

  # Average the probability of death across bundles for each leaf
  avg_leaf_probs <- setNames(round(vapply(leaf_paths, function(path) {
    mean(vapply(bundles, function(b) {
      as.numeric(b$pruned_leaf_probs[[path]])
    }, numeric(1)), na.rm = TRUE)
  }, numeric(1)), 3), leaf_paths)

  # Average the number of patients in each leaf across bundles
  avg_leaf_n <- setNames(as.integer(round(vapply(leaf_paths, function(path) {
    mean(vapply(bundles, function(b) {
      length(b$pruned_leaf_members[[path]])
    }, integer(1)), na.rm = TRUE)
  }, numeric(1)))), leaf_paths)

  # Build the averaged leaf table
  pruned_leaf_df <- data.table(Path = leaf_paths, n = as.integer(avg_leaf_n[leaf_paths]),
                               p = as.numeric(avg_leaf_probs[leaf_paths]))

  # Create dummy patient membership vectors so leaf sizes remain available for plotting
  pruned_leaf_members <- setNames(lapply(leaf_paths, function(path) {
    seq_len(max(1L, as.integer(avg_leaf_n[[path]])))
  }), leaf_paths)

  # Start from the template node table and standardise key columns
  pruned_node_df <- copy(template$pruned_node_df)
  pruned_node_df[, Path := as.character(Path)]
  pruned_node_df[, Predictor := as.character(Predictor)]
  pruned_node_df[, Symbol := as.character(Symbol)]
  pruned_node_df[, Value := as.character(Value)]

  # Return the averaged tree bundle
  list(pruned_node_df = pruned_node_df, pruned_leaf_df = pruned_leaf_df,
       pruned_leaf_probs = avg_leaf_probs, pruned_leaf_members = pruned_leaf_members,
       thresh = template$thresh)
}

# ---- Function 16: Build a one-row diagnostic performance table ----
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
  data.table(n = n, TP = tp, FP = fp, TN = tn, FN = fn, sensitivity = round(sensitivity, 1),
             specificity = round(specificity, 1), accuracy = round(accuracy, 1))
}

# ---- Function 17: Run the analysis for one imputation ----
run_one_analysis <- function(dataset, split_file, imp, cfg) {
  split_info <- fread(file = split_file)

  # Call Function 2 (separate tree objects)
  tree_frames <- separate_tree(split_info)

  # Prepare objects for next functions
  node_df <- tree_frames$node_df
  leaf_df <- tree_frames$leaf_df

  # Call Function 3 (prepare imputed dataset)
  dataset <- prepare_imputed_dataset(dataset)

  # Prepare objects for next functions
  tree_features <- unique(node_df$Predictor)
  tree_features <- tree_features[tree_features %in% names(dataset)]
  cols_full_tree <- unique(c("patient_id", tree_features, "Died_28D"))
  cols_full_tree <- cols_full_tree[cols_full_tree %in% names(dataset)]
  my_subset <- na.omit(dataset[, ..cols_full_tree])

  if (nrow(my_subset) != nrow(dataset)) {
    stop(paste0("Complete-case filtering removed patients in imputation ",
                imp, ": ", nrow(dataset) - nrow(my_subset), " rows removed."))
  }

  # Call Function 4 (create comparison function)
  comparison_function <- create_comparison()

  # Call Function 5 (send patients down tree)
  send_patients <- send_patients_down_tree(node_df, leaf_df, comparison_function)

  # Call Function 6 (assign patients to leaves)
  original <- assign_patients_to_leaves(my_subset, send_patients, leaf_df)
  orig_leaf_probs <- original$leaf_probs
  orig_leaf_members <- original$leaf_members

  # Step 0: original tree, with all leaves still present
  step0_layout <- list(node_df = node_df, leaf_df = leaf_df, leaf_probs = orig_leaf_probs,
                       leaf_members = orig_leaf_members)

  # Build the step1 tree structure by removing small leaves and compacting the tree
  step1_structure <- compact_pruned_tree(node_df, original$leaf_members, cfg$min_n)

  # Route every patient through the step1 tree structure so removed leaves are absorbed
  step1_send <- send_patients_down_pruned_tree(node_df = step1_structure$node_df,
    pruned_paths_to_leaf = step1_structure$leaf_df$Path,
    pruned_node_paths = step1_structure$node_df$Path,
    comparison_function = comparison_function)

  step1_routed <- assign_patients_to_leaves(my_subset, step1_send, step1_structure$leaf_df)

  # Keep one step1 object that contains the tree plus the routed memberships
  step1_layout <- list(node_df = step1_structure$node_df, leaf_df = step1_structure$leaf_df,
                       leaf_probs = step1_routed$leaf_probs,
                       leaf_members = step1_routed$leaf_members)

  # Call Function 7 (prune tree) using the step1 tree, not the raw original tree
  pruned_tree <- prune_tree(my_subset = my_subset, leaf_probs = step1_layout$leaf_probs,
                            leaf_members = step1_layout$leaf_members, thresh = cfg$thresh,
                            min_n = cfg$min_n)

  # Build the pruned tree structure from the pruned memberships
  pruned_structure <- compact_pruned_tree(step1_layout$node_df, pruned_tree$members,
                                          cfg$min_n)

  # Call Function 10 (send patients down pruned tree)
  pruned_send <- send_patients_down_pruned_tree(node_df = pruned_structure$node_df,
    pruned_paths_to_leaf = pruned_structure$leaf_df$Path,
    pruned_node_paths = pruned_structure$node_df$Path,
    comparison_function = comparison_function)

  pruned_routed <- assign_patients_to_leaves(my_subset, pruned_send,
                                             pruned_structure$leaf_df)

  # Keep one pruned object that contains the tree plus the routed memberships
  pruned_layout <- list(node_df = pruned_structure$node_df,
                        leaf_df = pruned_structure$leaf_df,
                        leaf_probs = pruned_routed$leaf_probs,
                        leaf_members = pruned_routed$leaf_members)

  # Call Function 9 (make signature for pruned tree)
  pruned_signature <- make_pruned_tree_signature(pruned_layout$node_df,
                                                 pruned_layout$node_df$Path)

  # Call Function 11 (classify imputed cohort using original tree)
  pred_full_orig <- classify_with_tree(my_subset, orig_leaf_probs,
                                       route_function = send_patients, thresh = cfg$thresh)

  # Call Function 11 (classify imputed cohort using pruned tree)
  pred_full_pruned <- classify_with_tree(my_subset, pruned_tree$probs,
                                         route_function = pruned_send,
                                         thresh = cfg$thresh)

  # State output file names for this imputation
  diag_file <- file.path(DIAG_DIR, sprintf("Imputed_diagnostic_performance_%02d.csv", imp))
  rds_file <- file.path(RDS_DIR, sprintf("Imputed_tree_analysis_result_%02d.rds", imp))

  # Call Function 12 (build diagnostic performance table)
  diagnostic_performance <- rbindlist(list(build_summary(pred_full_orig, "original",
                                                         "imputed cohort"),
                                           build_summary(pred_full_pruned, "pruned",
                                                         "imputed cohort")))

  # Export diagnostic performance table
  fwrite(as.data.table(diagnostic_performance), diag_file)

  # Export an RDS bundle containing objects needed by the plotting script
  saveRDS(list(step0_node_df = step0_layout$node_df, step0_leaf_df = step0_layout$leaf_df,
               step0_leaf_probs = step0_layout$leaf_probs,
               step0_leaf_members = step0_layout$leaf_members,
               step1_node_df = step1_layout$node_df, step1_leaf_df = step1_layout$leaf_df,
               step1_leaf_probs = step1_layout$leaf_probs,
               step1_leaf_members = step1_layout$leaf_members,
               pruned_node_df = pruned_layout$node_df,
               pruned_leaf_df = pruned_layout$leaf_df,
               pruned_leaf_probs = pruned_layout$leaf_probs,
               pruned_leaf_members = pruned_layout$leaf_members, thresh = cfg$thresh),
          file = rds_file)

  # Return analysis objects without printing them to the console
  invisible(list(diagnostic_performance = as.data.table(diagnostic_performance),
                 pred_full_orig = as.data.table(pred_full_orig),
                 pred_full_pruned = as.data.table(pred_full_pruned),
                 pruned_signature = pruned_signature, diag_file = diag_file,
                 rds_file = rds_file))
}

# ======================
# RUN THE FULL ANALYSIS
# ======================

# Call Function 1 (load imputed dataset)
mi_list <- load_imputed_long_dataset(long_path = IMPUTED_LONG_PATH, n_imputations = 20L)

analysis_results <- vector("list", length(mi_list))

# Loop over the 20 imputations and run the full analysis for each one
for (imp in seq_along(mi_list)) {
  dataset <- mi_list[[imp]]
  split_file <- file.path(SPLITS_DIR, sprintf("Imputed_dataset_splits_%02d.csv", imp))

  # Call Function 17 (run the analysis for one imputation)
  analysis_results[[imp]] <- run_one_analysis(dataset, split_file, imp, CFG)
  message("Finished imputation ", imp)
}

# Create a table linking each imputation to its diagnostic file and pruned-tree signature
sig_dt <- rbindlist(lapply(seq_along(analysis_results), function(i) {
  data.table(imp = i, diag_file = analysis_results[[i]]$diag_file,
             rds_file = analysis_results[[i]]$rds_file,
             signature = analysis_results[[i]]$pruned_signature)
}))

# Assign a group number to each unique tree signature and keep groups with >1 tree
sig_dt[, group_id := match(signature, unique(signature))]
group_sizes <- sig_dt[, .N, by = group_id][N > 1]

# Loop over repeated tree structures and export an averaged diagnostic table and pruned tree
for (g in seq_len(nrow(group_sizes))) {
  gid <- group_sizes$group_id[g]

  diag_files <- sig_dt[group_id == gid, diag_file]
  rds_files <- sig_dt[group_id == gid, rds_file]

  # Call Function 14 (build an average diagnostic table for a group of trees)
  group_table <- build_group_diagnostic_table(diag_files)
  fwrite(group_table, file.path(DIAG_DIR, sprintf("Imputed_diagnostic_perf_group_%d.csv",
                                                  g)))

  # Call Function 15 (build an averaged pruned-tree bundle for one repeated group)
  group_bundle <- build_group_pruned_bundle(rds_files)
  saveRDS(group_bundle, file.path(RDS_DIR,
                                  sprintf("Imputed_tree_analysis_result_group_%d.rds", g)))
}

# Run severe-malaria criteria table on imputation 1 only
dataset1 <- mi_list[[1]]
dataset1 <- prepare_imputed_dataset(dataset1)

# Create a complete-case dataset containing the severe malaria criteria variables
cols_criteria <- unique(c("patient_id", "Coma", "Respiratory distress", "BUN",
                          "Bicarbonate_mEq_L", "Lactate", "para_ul", "Glucose_mmol_L", "Hb",
                          "Died_28D"))
cols_criteria <- cols_criteria[cols_criteria %in% names(dataset1)]
criteria_subset <- na.omit(dataset1[, ..cols_criteria])

# Call Function 13 (classify patients using severe malaria criteria)
pred_criteria <- classify_with_criteria(criteria_subset)

# Call Function 16 (build a one-row diagnostic performance table)
comparison_summary <- build_comparison_summary(pred_criteria)

fwrite(comparison_summary, file.path(DIAG_DIR, "Imputed_diagnostic_perf_comparison.csv"))

message("Finished imputed analysis.")