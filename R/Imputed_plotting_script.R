# This script loads an RDS bundle for every imputation (20 imputations in total),
# containing objects from the GOSDT tree analysis script:
# - tables storing the paths, predictors, symbols, values, and node/leaf types in the tree,
# - maps storing probabilities and patient memberships for leaves,
# - and the probability threshold used to colour leaf nodes.
# It plots each tree before pruning (step 0) and after pruning.

library(DiagrammeR)
library(DiagrammeRsvg)
library(data.table)
options(scipen = 999)

# State paths of output folders
RDS_DIR <- file.path("R", "Imputed_tree_analysis_result")
OUT_DIR <- "outputs"
PRUNED_DIR <- file.path(OUT_DIR, "Imputed_pruned_tree")
STEP0_DIR <- file.path(OUT_DIR, "Imputed_step0_tree")

# If any output folder does not exist then create it
for (d in c(OUT_DIR, PRUNED_DIR, STEP0_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# ---- Function 1: Clean text so it can be written safely inside Graphviz labels ----
clean_text <- function(x) {

  # Convert the input to text
  x <- as.character(x)

  # Replace missing values with empty text
  x[is.na(x)] <- ""

  # Replace each backslash with a double backslash
  x <- gsub("\\\\", "\\\\\\\\", x)

  # Replace quotation marks with escaped quotation marks
  x <- gsub('"', '\\"', x)

  # Replace line breaks with spaces
  x <- gsub("\n", " ", x)
  x
}

# ---- Function 2: Safely convert a single value to numeric ----
convert_num <- function(x, default = NA_real_) {

  # If the input is empty then return the default value
  if (is.null(x) || length(x) == 0) return(default)

  # Convert the first input value to numeric
  x <- as.numeric(x)[1]

  # If conversion fails then return the default value
  # Otherwise return the numeric value
  if (is.na(x)) default else x
}

# ---- Function 3: Safely read one named value from a list or named vector ----
read_named_value <- function(x, key, default = NA_real_) {

  # If object is empty or requested name does not exist, then return default value
  if (is.null(x) || length(x) == 0 || is.null(names(x)) || !key %in% names(x)) {
    return(default)
  }

  # Extract the value associated with the requested name
  val <- x[[key]]

  # If the extracted value is empty then return the default value
  if (is.null(val) || length(val) == 0) return(default)

  # Call Function 2 (safely convert value to numeric)
  convert_num(val, default = default)
}

# ---- Function 4: Safely read the length of a named list entry ----
read_named_length <- function(x, key, default = 0L) {
  if (is.null(x) || length(x) == 0 || is.null(names(x)) || !key %in% names(x)) {
    return(default)
  }
  val <- x[[key]]
  if (is.null(val) || length(val) == 0) return(default)

  # Return the number of elements in the extracted value
  as.integer(length(val))
}

# ---- Function 5: Standardise comparison symbols for display ----
normalise_symbol <- function(symbol) {

  # Convert the symbol to a character string and keep only the first value
  symbol <- as.character(symbol)[1]
  if (is.na(symbol) || !nzchar(symbol)) return("")

  # Replace double equals with a single equals sign
  if (symbol == "==") return("=")
  symbol
}

# ---- Function 6: Format a threshold value ----
format_value <- function(v) {

  # Convert the value to a character string
  v <- as.character(v)[1]

  # If the value is missing or empty, return an empty string
  if (is.na(v) || !nzchar(v)) return("")

  # Convert value to numeric
  num <- suppressWarnings(as.numeric(v))
  if (!is.na(num)) {

    # If the number is a whole number, show it without decimal places
    if (abs(num - round(num)) < 1e-8) {
      return(as.character(as.integer(round(num))))
    }

    # Return value in plain text without scientific notation
    return(format(num, scientific = FALSE, trim = TRUE))
  }

  v
}

# ---- Function 7: Replace predictor name with display name ----
display_predictor_name <- function(predictor) {
  predictor <- as.character(predictor)[1]
  if (is.na(predictor) || !nzchar(predictor)) return("")

  if (predictor == "Glucose_mmol_L") return("Glucose")
  if (predictor == "para_ul") return("Parasitaemia")
  if (predictor == "Bicarbonate_mEq_L") return("Bicarbonate")
  if (predictor == "Temperature_C") return("Temperature")
  if (predictor == "BUN") return("BUN")
  if (predictor == "Lactate") return("Lactate")
  if (predictor == "Hb") return("Hb")
  if (predictor == "Respiratory_distress") return("Respiratory distress")

  predictor
}

# ---- Function 8: Fetch the display unit for a predictor ----
fetch_unit <- function(predictor) {

  # Call Function 7 (replace predictor name with display name)
  predictor <- display_predictor_name(predictor)

  switch(predictor, "BUN" = "mmol/L", "Lactate" = "mmol/L", "Glucose" = "mmol/L",
         "Bicarbonate" = "mEq/L", "Temperature" = "°C", "Parasitaemia" = "parasites/µL",
         "Hb" = "g/dL", "")
}

# ---- Function 9: Attach the unit to a split label ----
attach_unit <- function(label, predictor) {

  # Call Function 8 (fetch unit of predictor)
  unit <- fetch_unit(predictor)

  if (!nzchar(unit)) {
    return(label)
  }
  paste0(label, " ", unit)
}

# Function 10: Create a readable split label from predictor, symbol, and threshold value
format_split_label <- function(predictor, symbol, value) {

  # Call Function 7 (replace predictor name with display name)
  predictor <- display_predictor_name(predictor)
  if (is.na(predictor) || !nzchar(predictor)) return("Root")

  # Call Function 5 (normalise symbol)
  symbol <- normalise_symbol(symbol)

  # Call Function 6 (format threshold value)
  value <- format_value(value)

  # If either part is missing, show only the predictor name
  if (!nzchar(symbol) || !nzchar(value)) {
    return(predictor)
  }

  # Combine the predictor, symbol, and threshold value into one label
  label <- sprintf("%s %s %s", predictor, symbol, value)

  # Call Function 9 (attach unit)
  attach_unit(label, predictor)
}

# ---- Function 11: Associate each node with a path ----
node_id <- function(path) {
  path <- as.character(path)

  # If path is empty then return root identifier
  if (identical(path, "")) {
    "root"

    # Replace non-alphanumeric characters with underscores
  } else {
    paste0("n_", gsub("[^A-Za-z0-9]", "_", path))
  }
}

# ---- Function 12: Associate each leaf with a path ----
leaf_id <- function(path) {
  paste0("leaf_", gsub("[^A-Za-z0-9]", "_", as.character(path)))
}

# ---- Function 13: Colour a leaf based on whether the probability of death is above/below
# the threshold (high-risk leaves are red and low-risk leaves are green) ----
prob_fill <- function(p, thresh) {

  # Call Function 2 (safely convert probability to numeric)
  p <- convert_num(p, default = NA_real_)
  if (is.na(p)) return("#cccccc")
  if (p >= thresh) "#d73027" else "#2ca25f"
}

# ---- Function 14: Choose black or white text based on the leaf colour ----
text_colour_for_fill <- function(fill) {

  # If no valid fill colour is available, use black text
  if (is.null(fill) || length(fill) == 0 || all(is.na(fill)) || !nzchar(fill[1])) {
    return("#000000")
  }

  # Convert the fill colour to RGB and calculate its brightness
  rgb <- grDevices::col2rgb(fill[1]) / 255
  lum <- 0.299 * rgb[1] + 0.587 * rgb[2] + 0.114 * rgb[3]

  # If the fill colour is in the darker half of the spectrum, choose white text
  # Otherwise choose black text
  if (lum < 0.55) "#FFFFFF" else "#000000"
}

# ---- Function 15: Return the node information given a path ----
return_node_info <- function(node_df, path) {
  node_df <- as.data.table(node_df)
  node_df[, Path := as.character(Path)]
  node_df[, Predictor := as.character(Predictor)]
  node_df[, Symbol := as.character(Symbol)]
  node_df[, Value := as.character(Value)]

  # Find the row that matches the requested path
  nd <- node_df[Path == as.character(path), , drop = FALSE]
  if (nrow(nd) == 0) return(data.table())

  data.table(predictor = as.character(nd$Predictor[1]),
             symbol = normalise_symbol(nd$Symbol[1]),
             value = as.character(nd$Value[1]))
}

# ---- Function 16: Return the leaf statistics (n and probability) for a path ----
leaf_stats <- function(leaf_probs, leaf_members, path) {

  # Call Function 4 (safely read number of patients in leaf)
  list(n = read_named_length(leaf_members, path, default = 0L),

       # Call Function 3 (safely read leaf probability)
       p = read_named_value(leaf_probs, path, default = NA_real_))
}

# ---- Function 17: Export a DiagrammeR graph as a PDF via SVG ----
export_tree_pdf <- function(viz, file) {
  if (!requireNamespace("rsvg", quietly = TRUE)) {
    stop(paste0("Package 'rsvg' is required for PDF export.\n",
                "Install it with install.packages('rsvg')."))
  }

  # Convert the Graphviz figure to SVG format
  svg_txt <- DiagrammeRsvg::export_svg(viz)
  tmp_svg <- tempfile(fileext = ".svg")
  writeLines(svg_txt, tmp_svg)
  on.exit(unlink(tmp_svg), add = TRUE)

  # Save the SVG as a PDF file
  rsvg::rsvg_pdf(tmp_svg, file)
}

# ---- Function 18: Draw a tree exactly as stored in the relevant node and leaf tables ----
draw_tree_graph <- function(node_df, leaf_paths, leaf_probs, leaf_members, thresh,
                            root_gap = 0, base_nodesep = 0.20, base_ranksep = 0.50,
                            t_label_distance_by_depth = c("0" = 3, "1" = 3, "2" = 3,
                                                          "3" = 3),
                            f_label_distance_by_depth = c("0" = 3, "1" = 3, "2" = 3,
                                                          "3" = 3),
                            t_label_angle_by_depth = c("0" = 12, "1" = 14, "2" = 14,
                                                       "3" = 14),
                            f_label_angle_by_depth = c("0" = -12, "1" = -14, "2" = -14,
                                                       "3" = -14)) {

  # Convert the node table into a data table and standardise key columns
  node_df <- as.data.table(node_df)
  node_df[, Path := as.character(Path)]
  node_df[, Predictor := as.character(Predictor)]
  node_df[, Symbol := as.character(Symbol)]
  node_df[, Value := as.character(Value)]

  # Remove missing values from the leaf paths
  leaf_paths <- as.character(leaf_paths)
  leaf_paths <- unique(leaf_paths[!is.na(leaf_paths)])

  # Stop if there are no leaf paths left to draw
  if (length(leaf_paths) == 0) {
    stop("No leaves remain after filtering.")
  }

  # Identify every node path that must remain visible because it lies on a path
  # from the root to one of the remaining leaves
  visible_node_paths <- c("")
  for (p in leaf_paths) {
    if (nchar(p) > 1) {
      visible_node_paths <- c(visible_node_paths,
                              vapply(1:(nchar(p) - 1), function(k) substr(p, 1, k),
                                     character(1)))
    }
  }

  # Keep only these node paths
  visible_node_paths <- unique(visible_node_paths)
  node_df <- node_df[node_df$Path %in% visible_node_paths, , drop = FALSE]

  # Function 18a: Look up layout values by depth
  get_depth_value <- function(x, depth, default = NA_real_) {
    key <- as.character(depth)
    if (!is.null(names(x)) && key %in% names(x)) {
      val <- suppressWarnings(as.numeric(x[[key]]))
      if (!is.na(val)) return(val)
    }
    if (!is.null(names(x)) && length(names(x)) > 0) {
      last_key <- names(x)[length(names(x))]
      val <- suppressWarnings(as.numeric(x[[last_key]]))
      if (!is.na(val)) return(val)
    }
    default
  }

  # Function 18b: Add a node to the tree drawing
  add_node <- function(id, label, shape = "box", fill = "#f2f2f2", fontcolor = "#000000") {

    # Call Function 1 (clean label text)
    sprintf('%s [label = "%s", shape = %s, style = "filled", fillcolor = "%s",
            fontcolor = "%s"];', id, clean_text(label), shape, fill, fontcolor)
  }

  # Function 18c: Add a directed edge with a T/F label
  add_edge <- function(from, to, lbl, depth) {
    label_dist <- if (lbl == "F") {
      get_depth_value(f_label_distance_by_depth, depth, default = 1.05)
    } else {
      get_depth_value(t_label_distance_by_depth, depth, default = 1)
    }

    # Use a different label position depending on whether the edge is T or F
    label_angle <- if (lbl == "F") {
      get_depth_value(f_label_angle_by_depth, depth, default = -20)
    } else {
      get_depth_value(t_label_angle_by_depth, depth, default = 20)
    }

    # Graphviz code for this edge
    sprintf('%s -> %s [taillabel = "%s", labeldistance = %.2f, labelangle = %d,\n
            labelfloat = true, fontname = Helvetica, fontsize = 10, fontcolor = "#666666",\n
            color = "#666666", penwidth = 1.1];', from, to, lbl, label_dist, label_angle)
  }

  # Function 18d: Add an invisible node to keep branches aligned
  add_dummy_node <- function(id) {
    sprintf('%s [label = "", shape = point, width = 0.01, style = invis];', id)
  }

  # Function 18e: Add an invisible node so the first split can be spaced more neatly
  add_first_split_gap <- function(id, width) {
    sprintf('%s [label = "", shape = box, style = invis, fixedsize = true, width = %.2f,\n
            height = 0.01];', id, width)
  }

  # Function 18f: Draw a path in the tree
  draw_path <- function(path, is_root = FALSE, depth = 0L) {

    # Call Function 15 (return node information given path)
    this_term <- return_node_info(node_df, path)

    # Try extending the path to the left and right child nodes
    left_path <- paste0(path, "L")
    right_path <- paste0(path, "R")

    # Check for existence
    left_exists <- left_path %in% visible_node_paths || left_path %in% leaf_paths
    right_exists <- right_path %in% visible_node_paths || right_path %in% leaf_paths
    n_kids <- as.integer(left_exists) + as.integer(right_exists)

    # If current node is a leaf
    if (n_kids == 0) {
      if (!(path %in% leaf_paths)) return(NULL)

      # Call Function 16 (return leaf statistics)
      s <- leaf_stats(leaf_probs, leaf_members, path)

      # Call Function 13 (colour leaf based on probability)
      fill <- prob_fill(s$p, thresh)

      # Call Function 14 (choose text colour based on colour of leaf)
      font <- text_colour_for_fill(fill)

      # Build the text shown inside the leaf
      lab <- if (is.na(s$p)) {
        sprintf("n = %d\\np = NA", s$n)
      } else {
        sprintf("n = %d\\np = %.3f", s$n, s$p)
      }

      # Call Function 12 (create a path identifier for the leaf)
      leaf_dot_id <- leaf_id(path)

      # Graphviz code for this leaf
      dot <<- c(dot, sprintf('%s [label = "%s", shape = triangle, style = "filled",\n
                             fillcolor = "%s", fontcolor = "%s", fixedsize = true,\n
                             width = 1.35, height = 1.05, margin = 0.02];',
                             leaf_dot_id, lab, fill, font))

      return(leaf_dot_id)
    }

    # Call Function 11 (create a path identifier for a non-leaf node)
    this_id <- if (is_root) "root" else node_id(path)

    # Create the display label for the node
    this_lab <- if (nrow(this_term) == 0) {
      "Root"
    } else {

      # Call Function 10 (create readable split label)
      format_split_label(this_term$predictor[1], this_term$symbol[1], this_term$value[1])
    }

    # Add the current node label to the tree drawing
    dot <<- c(dot, add_node(this_id, this_lab, shape = "box", fill = "#f2f2f2"))

    left_id <- NULL
    right_id <- NULL

    # Call Function 18f recursively (draw the left child subtree)
    if (left_exists) {
      left_id <- draw_path(left_path, is_root = FALSE, depth = depth + 1L)
    }

    # Call Function 18f recursively (draw the right child subtree)
    if (right_exists) {
      right_id <- draw_path(right_path, is_root = FALSE, depth = depth + 1L)
    }

    # If only one child exists, add an invisible dummy node so the tree stays balanced
    if (!is.null(left_id) && is.null(right_id)) {
      dummy_id <- paste0(this_id, "_dummyR")
      dot <<- c(dot, add_dummy_node(dummy_id),
                sprintf('{ rank = same; %s; %s; }', left_id, dummy_id),
                sprintf('%s -> %s [style = invis, weight = 100];', left_id, dummy_id))
    }

    if (is.null(left_id) && !is.null(right_id)) {
      dummy_id <- paste0(this_id, "_dummyL")
      dot <<- c(dot, add_dummy_node(dummy_id),
                sprintf('{ rank = same; %s; %s; }', dummy_id, right_id),
                sprintf('%s -> %s [style = invis, weight = 100];', dummy_id, right_id))
    }

    # Add extra spacing after the root if requested
    if (depth == 0L && root_gap > 0 && !is.null(left_id) && !is.null(right_id)) {
      gap_id <- paste0(this_id, "_rootgap")

      dot <<- c(dot, add_first_split_gap(gap_id, root_gap),
                sprintf('{ rank=same; %s; %s; %s; }', left_id, gap_id, right_id),
                sprintf('%s -> %s [style = invis, weight = 1000, minlen = 0];',
                        left_id, gap_id),
                sprintf('%s -> %s [style = invis, weight = 1000, minlen = 0];',
                        gap_id, right_id))
    }

    # Add the true edges from this node to its left child
    if (!is.null(left_id)) {
      dot <<- c(dot, add_edge(this_id, left_id, "T", depth))
    }

    # Add the false edge from this node to its right child
    if (!is.null(right_id)) {
      dot <<- c(dot, add_edge(this_id, right_id, "F", depth))
    }

    this_id
  }

  # Start building the Graphviz document
  dot <- c("digraph G {", sprintf("graph [rankdir = TB, ordering = out,
                                  splines = polyline,\n newrank = true, nodesep = %.2f,
                                  ranksep = %.2f, pad = 0.35,\n margin = 0.05];",
                                  base_nodesep, base_ranksep),
           "node [fontname = Helvetica, fontsize = 11, margin = 0.08];",
           "edge [arrowhead = none, color = \"#666666\", penwidth = 1.1,\n
           fontname = Helvetica, fontsize = 10];")

  # Draw the tree starting from the root
  draw_path("", is_root = TRUE, depth = 0L)

  # Finish the Graphviz document
  dot <- c(dot, "}")
  grViz(paste(dot, collapse = "\n"))
}

# Root spacing values used to improve the layout
pruned_root_gap <- 0
step0_root_gap <- 3

# --- Function 19: Render one imputation ---
render_one_result <- function(result_rds, step0_tree_pdf, pruned_tree_pdf) {

  # Load the saved objects from the analysis script
  obj <- readRDS(result_rds)
  step0_node_df <- as.data.table(obj$step0_node_df)
  step0_leaf_df <- as.data.table(obj$step0_leaf_df)
  step0_leaf_probs <- obj$step0_leaf_probs
  step0_leaf_members <- obj$step0_leaf_members

  pruned_node_df <- as.data.table(obj$pruned_node_df)
  pruned_leaf_df <- as.data.table(obj$pruned_leaf_df)
  pruned_leaf_probs <- obj$pruned_leaf_probs
  pruned_leaf_members <- obj$pruned_leaf_members

  THRESH <- obj$thresh

  # Call Function 18 (draw the step-0 tree)
  viz0 <- draw_tree_graph(node_df = step0_node_df,
                          leaf_paths = step0_leaf_df$Path,
                          leaf_probs = step0_leaf_probs,
                          leaf_members = step0_leaf_members,
                          thresh = THRESH,
                          root_gap = step0_root_gap)

  # Call Function 17 (export tree pdf)
  export_tree_pdf(viz0, step0_tree_pdf)

  # Call Function 18 (draw the pruned tree)
  viz <- draw_tree_graph(node_df = pruned_node_df,
                         leaf_paths = pruned_leaf_df$Path,
                         leaf_probs = pruned_leaf_probs,
                         leaf_members = pruned_leaf_members,
                         thresh = THRESH,
                         root_gap = pruned_root_gap)

  # Call Function 17 (export tree pdf)
  export_tree_pdf(viz, pruned_tree_pdf)
}

# ---- Function 20: Render one averaged pruned group tree ----
render_group_result <- function(result_rds, pruned_tree_pdf) {

  obj <- readRDS(result_rds)

  pruned_node_df <- as.data.table(obj$pruned_node_df)
  pruned_leaf_df <- as.data.table(obj$pruned_leaf_df)
  pruned_leaf_probs <- obj$pruned_leaf_probs
  pruned_leaf_members <- obj$pruned_leaf_members
  THRESH <- obj$thresh

  viz <- draw_tree_graph(node_df = pruned_node_df, leaf_paths = pruned_leaf_df$Path,
                         leaf_probs = pruned_leaf_probs, leaf_members = pruned_leaf_members,
                         thresh = THRESH, root_gap = pruned_root_gap)

  export_tree_pdf(viz, pruned_tree_pdf)
}

# Run the main code for each of 20 imputations
for (imp in 1:20) {

  # State paths of output files
  result_rds <- file.path(RDS_DIR, sprintf("Imputed_tree_analysis_result_%02d.rds", imp))
  step0_tree_pdf <- file.path(STEP0_DIR, sprintf("Imputed_step0_tree_%02d.pdf", imp))
  pruned_tree_pdf <- file.path(PRUNED_DIR, sprintf("Imputed_pruned_tree_%02d.pdf", imp))

  # Call Function 19 (render one imputation)
  render_one_result(result_rds, step0_tree_pdf, pruned_tree_pdf)

  message("Saved plots for imputation ", imp)
}


# Render one averaged pruned tree for each repeated tree group
group_files <- list.files(RDS_DIR,
                          pattern = "^Imputed_tree_analysis_result_group_[0-9]+\\.rds$",
                          full.names = TRUE)

for (result_rds in sort(group_files)) {

  group_num <- sub("^.*group_([0-9]+)\\.rds$", "\\1", basename(result_rds))

  pruned_tree_pdf <- file.path(PRUNED_DIR, sprintf("Imputed_pruned_tree_group_%s.pdf",
                                                   group_num))

  # Call Function 20 (render one averaged pruned group tree)
  render_group_result(result_rds, pruned_tree_pdf)

  message("Saved group plot ", group_num)
}