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
VERBOSE = False
ALLOW_SMALL_REG = True
BALANCE = True

# Validation settings
TEST_SIZE = 0.20
RANDOM_STATES = range(123, 143)

N_IMPUTATIONS = 20
OUTCOME = "Died_28D:=1"

# Main function: Run repeated validation of the GOSDT algorithm
def gosdt_repeated_validation():

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

        # Call Function 4 for every patient
        predictions = [predict_one(node, row) for _, row in features.iterrows()]

        return predictions

    # Function 6: Calculate model performance
    def evaluate_model(node, features, target, dataset_name, repeat, random_state, imp):

        # Call Function 5 (predict outcomes)
        predictions = predict_dataset(node, features)

        # Calculate the confusion matrix
        tn, fp, fn, tp = confusion_matrix(target, predictions,
                                           labels = [0, 1]).ravel()

        # Return one row containing the validation results
        return pd.DataFrame({"Repeat": [repeat], "RandomState": [random_state],
                             "Imputation": [imp],
                             "Dataset": [dataset_name], "N": [len(target)],
                             "Accuracy": [accuracy_score(target, predictions)],
                             "BalancedAccuracy":
                                 [balanced_accuracy_score(target, predictions)],
                             "Sensitivity":
                                 [recall_score(target, predictions, pos_label = 1,
                                               zero_division = 0)],
                             "Specificity":
                                 [recall_score(target, predictions, pos_label = 0,
                                               zero_division = 0)],
                             "Precision":
                                 [precision_score(target, predictions, pos_label = 1,
                                                  zero_division = 0)],
                             "F1":
                                 [f1_score(target, predictions, pos_label = 1,
                                           zero_division = 0)],
                             "TN": [tn], "FP": [fp], "FN": [fn], "TP": [tp]})

    # Read imputation 1 to define the common patient population
    reference_file = ("Python/Imputed_dataset_binned/Imputed_dataset_binned_01.csv")
    reference_data = pd.read_csv(reference_file)
    reference_ids = reference_data["patient_id"]
    reference_target = reference_data[OUTCOME]

    # Store results from every repeat and imputation
    all_results = []

    # Repeat the stratified 80:20 patient-level validation
    for repeat, random_state in enumerate(RANDOM_STATES, start = 1):

        print(f"\nRunning repetition {repeat} of {len(RANDOM_STATES)} "
              f"(random state = {random_state})")

        # Create one patient split for this repetition
        # Exactly the same patients are used for all 20 imputations
        train_ids, test_ids = train_test_split(reference_ids,
                                               test_size = TEST_SIZE,
                                               random_state = random_state,
                                               stratify = reference_target)
        train_ids = set(train_ids)
        test_ids = set(test_ids)

        for imp in range(1, N_IMPUTATIONS + 1):
            print(f"  Imputation {imp:02d}/{N_IMPUTATIONS}")

            input_file = f"Python/Imputed_dataset_binned/Imputed_dataset_binned_{imp:02d}.csv"
            data = pd.read_csv(input_file)

            # Apply the same patient split to this imputation
            train_data = data[data["patient_id"].isin(train_ids)].copy()
            test_data = data[data["patient_id"].isin(test_ids)].copy()

            # Remove metadata and outcome from predictors
            X_train = train_data.drop(columns = ["patient_id",
                                                 "imputation_id",
                                                 OUTCOME])
            y_train = train_data[OUTCOME]

            X_test = test_data.drop(columns = ["patient_id",
                                               "imputation_id",
                                               OUTCOME])
            y_test = test_data[OUTCOME]

            # Initialise the balanced GOSDT model
            gosdt_model = GOSDTClassifier(regularization = REGULARIZATION,
                                          similar_support = SIMILAR_SUPPORT,
                                          time_limit = TIME_LIMIT,
                                          depth_budget = DEPTH_BUDGET,
                                          verbose = VERBOSE,
                                          allow_small_reg = ALLOW_SMALL_REG,
                                          balance = BALANCE)

            # Train on the training patients only
            gosdt_model.fit(X_train, y_train)

            tree0 = (gosdt_model.trees_[0].tree)

            # Evaluate training performance
            training_results = evaluate_model(tree0, X_train, y_train, "Training",
                                              repeat, random_state, imp)

            # Evaluate held-out performance
            test_results = evaluate_model(tree0, X_test, y_test, "Test",
                                          repeat, random_state, imp)

            all_results.append(training_results)
            all_results.append(test_results)

    # Combine all 800 validation results
    validation_results = pd.concat(all_results, ignore_index = True)

    # Save all individual results
    output_csv = "Python/Imputed_dataset_repeated_validation.csv"
    validation_results.to_csv(output_csv, index = False)
    print(f"\nSaved '{output_csv}'")

    # Keep held-out test results
    test_results = validation_results[validation_results["Dataset"] == "Test"].copy()

    metrics = ["Accuracy", "BalancedAccuracy", "Sensitivity", "Specificity",
               "Precision", "F1"]

    # Average across the 20 imputations within each train/test partition
    repeat_results = (test_results.groupby(["Repeat", "RandomState"],
                                           as_index = False)[metrics].mean())

    repeat_csv = "Python/Imputed_dataset_repeated_validation_by_repeat.csv"
    repeat_results.to_csv(repeat_csv, index = False)
    print(f"Saved '{repeat_csv}'")

    # Summarise variation across the 20 patient splits
    summary_results = (repeat_results[metrics].agg(["mean", "std", "median",
                                                    "min", "max"]).T.reset_index())

    summary_results.columns = ["Metric", "Mean", "SD", "Median", "Minimum", "Maximum"]

    summary_results["Q25"] = [repeat_results[metric].quantile(0.25)
                              for metric in metrics]

    summary_results["Q75"] = [repeat_results[metric].quantile(0.75)
                              for metric in metrics]

    summary_results = summary_results[["Metric", "Mean", "SD", "Median",
                                       "Q25", "Q75", "Minimum", "Maximum"]]

    summary_csv = "Python/Imputed_dataset_repeated_validation_summary.csv"
    summary_results.to_csv(summary_csv, index = False)

    print("\nRepeated held-out validation summary:")
    print(summary_results.to_string(index = False))
    print(f"\nSaved '{summary_csv}'")

# Run repeated MI validation
gosdt_repeated_validation()
