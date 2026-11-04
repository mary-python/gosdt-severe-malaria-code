import pandas as pd
from gosdt import GOSDTClassifier
import warnings
warnings.filterwarnings("ignore")
import re
import os
from sklearn.model_selection import train_test_split

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

# Main function: Run the GOSDT algorithm
def gosdt_function(input_csv, output_csv, train_ids, test_ids):

    # Read one imputed dataset
    data = pd.read_csv(input_csv, sep = ",")

    # Apply exactly the same patient-level training/test split
    # to every imputation
    train_data = data[data["patient_id"].isin(train_ids)].copy()
    test_data = data[data["patient_id"].isin(test_ids)].copy()

    # patient_id and imputation_id are metadata and must never
    # be supplied to GOSDT
    X_train = train_data.drop(columns = ["patient_id", "imputation_id", "Died_28D:=1"])
    y_train = train_data["Died_28D:=1"]
    X_test = test_data.drop(columns = ["patient_id", "imputation_id", "Died_28D:=1"])
    y_test = test_data["Died_28D:=1"]

    print(f"Training rows: {len(X_train)}")
    print(f"Test rows: {len(X_test)}")
    print(f"Candidate predictors: {X_train.shape[1]}")

    # Initialise the balanced GOSDT model
    gosdt_model = GOSDTClassifier(regularization = REGULARIZATION,
                                  similar_support = SIMILAR_SUPPORT,
                                  time_limit = TIME_LIMIT,
                                  depth_budget = DEPTH_BUDGET,
                                  verbose = VERBOSE,
                                  allow_small_reg = ALLOW_SMALL_REG,
                                  balance = BALANCE)

    # Train the GOSDT model on the training dataset only
    gosdt_model.fit(X_train, y_train)

    # Function 1: Extract dictionary from node
    def extract_dict(node):

        if node is None:
            return {}

        if isinstance(node, dict):
            return node

        if hasattr(node, "__dict__"):
            return vars(node)

        return {}

    # Function 2: Read one value from a node
    def read_node_value(node, *names, default = None):
        node_data = extract_dict(node)

        for name in names:
            if name in node_data:
                return node_data[name]

            if hasattr(node, name):
                return getattr(node, name)

        return default

    # Function 3: Check whether this node is a leaf
    def is_leaf(node):

        has_prediction = (read_node_value(node, "prediction") is not None)
        has_left_child = (read_node_value(node, "left_child") is not None)
        has_right_child = (read_node_value(node, "right_child") is not None)

        return (has_prediction or (not has_left_child and not has_right_child))

    # Function 4: Split a feature label into predictor name,
    # comparison symbol and threshold value
    def extract_split_details(label):
        label = str(label).strip()
        label = re.sub(r"\s*<=\s*0\.5\s*$", "", label)

        match = re.search(r"(.+?)(?:[:\s]*)(<=|>=|<|>|==|=)\s*([-\d\.]+)$", label)

        if not match:
            return label.strip(), "", ""

        predictor = match.group(1).strip()
        symbol = match.group(2).strip()
        value = match.group(3).strip()

        return predictor, symbol, value

    # Function 5: Get the left and right children
    def get_children(node):

        return (read_node_value(node, "left_child", default = None),
                read_node_value(node, "right_child", default = None))

    # Function 6: Convert tree to table
    def tree_to_table_rows(node, path = ""):
        rows = []

        if node is None:
            return rows

        if is_leaf(node):

            rows.append({"Path": path, "NodeType": "leaf", "Predictor": "",
                         "Symbol": "", "Value": "",
                         "Prediction": int(read_node_value(node, "prediction"))})

            return rows

        feat = read_node_value(node, "feature")

        if feat is not None:
            colname = str(X_train.columns[int(feat)])
        else:
            colname = ""

        predictor, symbol, value = (extract_split_details(colname))

        rows.append({"Path": path, "NodeType": "node", "Predictor": predictor,
                     "Symbol": symbol, "Value": value, "Prediction": ""})

        left_child, right_child = (get_children(node))

        rows.extend(tree_to_table_rows(left_child, path + "L"))
        rows.extend(tree_to_table_rows(right_child, path + "R"))

        return rows

    # Extract the first fitted tree
    tree0 = gosdt_model.trees_[0].tree

    rows = tree_to_table_rows(tree0)

    # Save tree structure
    data_out = pd.DataFrame(rows,
                            columns = ["Path", "NodeType", "Predictor",
                                       "Symbol", "Value", "Prediction"])

    data_out.to_csv(output_csv, index = False)

    print(f"Saved '{output_csv}'")


# Create one common patient-level train/test split
reference_file = "Python/Imputed_dataset_binned/Imputed_dataset_binned_01.csv"

reference_data = pd.read_csv(reference_file)

train_ids, test_ids = train_test_split(reference_data["patient_id"],
                                       test_size = TEST_SIZE,
                                       random_state = RANDOM_STATE,
                                       stratify = reference_data["Died_28D:=1"])

train_ids = set(train_ids)
test_ids = set(test_ids)

print(f"Training patients: {len(train_ids)}")
print(f"Test patients: {len(test_ids)}")

# Create output directory
directory = "Python/Imputed_dataset_splits"
os.makedirs(directory, exist_ok = True)

# Fit one independent GOSDT tree per imputation
for imp in range(1, N_IMPUTATIONS + 1):
    print(f"\nImputation {imp:02d}/{N_IMPUTATIONS}\n")

    input_csv = f"Python/Imputed_dataset_binned/Imputed_dataset_binned_{imp:02d}.csv"
    output_csv = f"Python/Imputed_dataset_splits/Imputed_dataset_splits_{imp:02d}.csv"

    gosdt_function(input_csv, output_csv, train_ids, test_ids)
