import pandas as pd
from gosdt import GOSDTClassifier
import warnings
warnings.filterwarnings("ignore")
from sklearn.model_selection import train_test_split
from sklearn.metrics import (accuracy_score, balanced_accuracy_score,
                             confusion_matrix, precision_score, recall_score,
                             f1_score)

# Hyperparameters for the model
REGULARIZATION = 1e-3
SIMILAR_SUPPORT = False
DEPTH_BUDGET = 4
TIME_LIMIT = 60
VERBOSE = True
ALLOW_SMALL_REG = True
BALANCE = True

# Validation settings
TEST_SIZE = 0.20
RANDOM_STATE = 123

N_IMPUTATIONS = 20
OUTCOME = "Died_28D:=1"

# Main function: Validate the GOSDT algorithm
def gosdt_validation(input_csv, train_ids, test_ids, imp):

    # Read one imputed dataset
    data = pd.read_csv(input_csv, sep = ",")

    # Apply exactly the same patient-level training/test split
    # to every imputation
    train_data = data[data["patient_id"].isin(train_ids)].copy()
    test_data = data[data["patient_id"].isin(test_ids)].copy()

    # patient_id and imputation_id are metadata and must never be supplied to GOSDT
    X_train = train_data.drop(columns = ["patient_id", "imputation_id", OUTCOME])
    y_train = train_data[OUTCOME]
    X_test = test_data.drop(columns = ["patient_id", "imputation_id", OUTCOME])
    y_test = test_data[OUTCOME]

    print(f"\nImputation {imp:02d}")
    print(f"Training rows: {len(X_train)}")
    print(f"Test rows: {len(X_test)}")
    print(f"Candidate predictors: {X_train.shape[1]}")

    # Initialise the GOSDT model
    gosdt_model = GOSDTClassifier(regularization = REGULARIZATION,
                                  similar_support = SIMILAR_SUPPORT,
                                  time_limit = TIME_LIMIT,
                                  depth_budget = DEPTH_BUDGET, verbose = VERBOSE,
                                  allow_small_reg = ALLOW_SMALL_REG, balance = BALANCE)

    # Train the GOSDT model on the training dataset only
    gosdt_model.fit(X_train, y_train)

    # Extract the first tree from the fitted GOSDT model
    tree0 = gosdt_model.trees_[0].tree

    # Function 1: Extract dictionary from node
    def extract_dict(node):

        # If node is empty then return an empty dictionary
        if node is None:
            return {}

        # If node is already a dictionary then return it
        if isinstance(node, dict):
            return node

        # If node is an object with a dictionary attribute then return this dictionary
        if hasattr(node, "__dict__"):
            return vars(node)

        return {}

    # Function 2: Read one value from a node
    def read_node_value(node, *names, default = None):

        # Call Function 1 (extract dictionary from node)
        data = extract_dict(node)

        # Try to match input name with each dictionary name
        for name in names:

            # If a matching name is found then return its value
            if name in data:
                return data[name]

            # If input name exists as attribute of node then return its value
            if hasattr(node, name):
                return getattr(node, name)

        # If no matching value is found, return the default value
        return default

    # Function 3: Check whether this node is a leaf
    def is_leaf(node):

        # Call Function 2 (get prediction and children of node)
        has_prediction = read_node_value(node, "prediction") is not None
        has_left_child = read_node_value(node, "left_child") is not None
        has_right_child = read_node_value(node, "right_child") is not None

        # The node is a leaf if it has a prediction and no children
        return has_prediction or (not has_left_child and not has_right_child)

    # Function 4: Predict the outcome of one patient using the fitted tree
    def predict_one(node, row):

        # Start at the root of the tree
        current_node = node

        # Continue until a leaf is reached
        while not is_leaf(current_node):

            # Call Function 2 (get feature and children of node)
            feature = int(read_node_value(current_node, "feature"))
            left_child = read_node_value(current_node, "left_child")
            right_child = read_node_value(current_node, "right_child")

            # Get the patient's value for this binary feature
            feature_value = row.iloc[feature]

            # Follow the appropriate branch
            # This direction reproduces the loss reported internally by GOSDT
            if feature_value <= 0.5:
                current_node = right_child
            else:
                current_node = left_child

        # Return the prediction at the leaf
        return int(read_node_value(current_node, "prediction"))

    # Function 5: Predict the outcomes of all patients in a dataset
    def predict_dataset(node, features):

        # Apply Function 4 to every patient
        predictions = [predict_one(node, row) for _, row in features.iterrows()]

        return predictions

    # Function 6: Calculate model performance
    def evaluate_model(node, features, target, dataset_name):

        # Call Function 5 (predict outcomes)
        predictions = predict_dataset(node, features)

        # Calculate the confusion matrix
        tn, fp, fn, tp = confusion_matrix(target, predictions,
                                          labels = [0, 1]).ravel()

        # Calculate validation statistics
        accuracy = accuracy_score(target, predictions)
        balanced_accuracy = balanced_accuracy_score(target, predictions)
        sensitivity = recall_score(target, predictions, pos_label = 1,
                                   zero_division = 0)
        specificity = recall_score(target, predictions, pos_label = 0,
                                   zero_division = 0)
        precision = precision_score(target, predictions, pos_label = 1,
                                    zero_division = 0)
        f1 = f1_score(target, predictions, pos_label = 1,
                      zero_division = 0)

        # Return one row containing the validation results
        return pd.DataFrame({"Imputation": [imp], "Dataset": [dataset_name],
                             "N": [len(target)], "Accuracy": [accuracy],
                             "BalancedAccuracy": [balanced_accuracy],
                             "Sensitivity": [sensitivity],
                             "Specificity": [specificity],
                             "Precision": [precision], "F1": [f1],
                             "TN": [tn], "FP": [fp], "FN": [fn], "TP": [tp]})

    # Call Function 6 for the training dataset
    training_results = evaluate_model(tree0, X_train, y_train, "Training")

    # Call Function 6 for the held-out test dataset
    test_results = evaluate_model(tree0, X_test, y_test, "Test")

    # Combine the training and test results
    validation_results = pd.concat([training_results, test_results],
                                   ignore_index = True)

    # Print the validation results
    print("\nValidation results:")
    print(validation_results.to_string(index = False))

    return validation_results

# Create the same common patient split as Imputed_gosdt.py
reference_file = "Python/Imputed_dataset_binned/Imputed_dataset_binned_01.csv"
reference_data = pd.read_csv(reference_file)

train_ids, test_ids = train_test_split(reference_data["patient_id"],
                                       test_size = TEST_SIZE,
                                       random_state = RANDOM_STATE,
                                       stratify = reference_data[OUTCOME])
train_ids = set(train_ids)
test_ids = set(test_ids)

print(f"Training patients: {len(train_ids)}")
print(f"Test patients: {len(test_ids)}")

# Validate all 20 imputations
all_results = []

for imp in range(1, N_IMPUTATIONS + 1):
    input_csv = f"Python/Imputed_dataset_binned/Imputed_dataset_binned_{imp:02d}.csv"
    results = gosdt_validation(input_csv, train_ids, test_ids, imp)
    all_results.append(results)

# Combine all validation results
validation_all = pd.concat(all_results, ignore_index = True)

# Save all 20 imputations
validation_file = "Python/Imputed_dataset_validation.csv"
validation_all.to_csv(validation_file, index = False)
print(f"\nSaved '{validation_file}'")

# Summarise validation performance across imputations
metric_columns = ["Accuracy", "BalancedAccuracy", "Sensitivity", "Specificity",
                  "Precision", "F1"]

validation_summary = (validation_all.groupby("Dataset")[metric_columns].agg(["mean",
                                                                             "std",
                                                                             "min",
                                                                             "max"]))
summary_file = "Python/Imputed_dataset_validation_summary.csv"
validation_summary.to_csv(summary_file)

print("\nValidation summary across imputations:")
print(validation_summary)
print(f"\nSaved '{summary_file}'")
