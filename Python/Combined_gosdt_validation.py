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

# Main function: Validate the GOSDT algorithm
def gosdt_validation(input_csv, output_csv):

    # Read the data from the CSV file
    data = pd.read_csv(input_csv, sep = ",")

    # Separate the input columns from the output column
    features, target = data.iloc[:, :-1], data.iloc[:, -1]

    # Features have already been binarised at prespecified thresholds
    features_with_thresholds = features.copy()

    # Split the dataset into the same training and test datasets used in
    # Combined_gosdt.py
    X_train, X_test, y_train, y_test = train_test_split(features_with_thresholds,
                                                        target,
                                                        test_size = TEST_SIZE,
                                                        random_state = RANDOM_STATE,
                                                        stratify = target)

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
        return pd.DataFrame({"Dataset": [dataset_name], "N": [len(target)],
                             "Accuracy": [accuracy],
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

    # Save the validation results as a CSV
    validation_results.to_csv(output_csv, index = False)
    print(f"Saved '{output_csv}'")

# Run the main function
gosdt_validation("Combined_dataset_binned.csv", "Combined_dataset_validation.csv")