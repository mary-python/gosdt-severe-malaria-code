# CREATE A TABLE CONTAINING A BINARY COLUMN FOR EACH SPLIT VALUE

# Load data.table package to read large CSV file
library(data.table)

# Turn off scientific notation
options(scipen = 999)

my_table <- fread("data/Combined_dataset copy.csv")

# Remove patients with gender not M/F
my_table <- my_table[SEX == "M" | SEX == "F"]

# Remove patients with unidentified outcome
my_table <- my_table[!is.na(Died_28D)]

# Remove patients older than 15 years
my_table <- my_table[Age_category != ">15"]

# Remove studies with fewer than 20 patients
my_table[, study_n := .N, by = STUDY]
my_table <- my_table[study_n >= 20]
my_table[, study_n := NULL]

# Reduce to threshold columns and remove NAs
my_table <- na.omit(my_table[, c(3, 12, 13, 14, 15, 17, 21, 44, 66, 67, 69)])

# Change respiratory distress from TRUE/FALSE to 1/0
my_table$`Respiratory distress` <- as.integer(my_table$`Respiratory distress`)

# Provide threshold names, symbols, and values
reps_of_each <- c(3, 1, 3, 3, 1, 3, 4, 3, 4, 1)
threshold_names <- rep(c("Hb", "Respiratory distress", "para_ul", "Temperature_C", "Coma",
                         "Lactate", "Glucose_mmol_L", "Bicarbonate_mEq_L", "BUN",
                         "Died_28D"),
                       times = reps_of_each)
threshold_symbols <- rep(c("<", "=", ">", "<", "=", ">", "<", "<", ">", "="),
                         times = reps_of_each)
threshold_values <- c(4, 5, 6, 1, 100000, 250000, 500000, 35, 37.5, 40, 1, 5, 7, 10,
                      2.2, 3, 3.5, 4, 10, 12, 15, 8, 10, 15, 20, 1)
num_thresholds <- length(threshold_names)

# Write full names of threshold columns
threshold_full <- 1:num_thresholds

for (i in 1:num_thresholds) {
  threshold_full[i] <- sprintf("%s:%s%s", threshold_names[i], threshold_symbols[i],
                               threshold_values[i])
}

# Create new data frame to store columns of binary data
Combined_dataset_binned <- data.frame(matrix(NA, nrow = nrow(my_table), ncol = num_thresholds))
names(Combined_dataset_binned) <- threshold_full

# Create function which increments itself
inc <- function(x) {
  eval.parent(substitute(x <- x + 1))
}

# Fill data frame with columns of 1s & 0s to describe presence/absence of each feature
for (i in 1:num_thresholds) {
  my_condition <- 1:nrow(my_table)
  j <- 0

  # Loop through elements of a feature of interest
  # Want to split them in two (0s & 1s) based on a condition
  for (x in my_table[[threshold_names[i]]]) {
    inc(j)

    # Handle missing values
    if (is.na(x)) {
      my_condition[j] = NA

    # Separate conditions according to symbol
    } else if (threshold_symbols[i] == "<") {

      # If condition true then set current element of vector to be 1
      if (x < threshold_values[i]) {
        my_condition[j] = 1

      # Else set current element of vector to be 0
      } else {
        my_condition[j] = 0
      }
    } else if (threshold_symbols[i] == ">") {
      if (x > threshold_values[i]) {
        my_condition[j] = 1
      } else {
        my_condition[j] = 0
      }
    } else {
      if (x == threshold_values[i]) {
        my_condition[j] = 1
      } else {
        my_condition[j] = 0
      }
    }
  }

  # Set the vector to be a column of new data frame
  Combined_dataset_binned[[threshold_full[i]]] <- my_condition
  Combined_dataset_binned
}

# State location where data frame will be saved
my_output <- "Python/Combined_dataset_binned.csv"

# Save data frame as CSV without row names
write.csv(Combined_dataset_binned, file = my_output, row.names = FALSE)