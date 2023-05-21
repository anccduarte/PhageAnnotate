# -*- coding: utf-8 -*-

import joblib
import pandas as pd
from pathlib import Path
from sklearn.ensemble import HistGradientBoostingClassifier as HGBR
from sklearn.metrics import accuracy_score, precision_score, recall_score
from sklearn.model_selection import train_test_split

"""
AS OF NOW, THE FINAL MODEL IS ONLY TRAINED ON 75% OF THE DATA. THE REMAINING
25% IS USED FOR TESTING. IN A MORE ADVANCED PHASE, BUILD MODEL ON 75% OF THE
DATA, TEST THE MODEL ON THE REMAINING 25%, AND RE-BUILD THE MODEL USING ALL
DATA AVALILABLE (???)
"""

def build_model_and_display_stats(file_name: str) -> None:
    """
    Builds an HGBR model fed on the data present in 'database/<file_name>.csv',
    saves it to a .csv file and prints some useful statisitics (accuracy,
    precision and recall).
    
    Parameters
    ----------
    file_name: str
        The name of the .csv file containing the data
    """
    # read data
    data = pd.read_csv(f"../database/{file_name}.csv")
    # split features from label
    X, y = data.iloc[:, 1:-1], data.iloc[:, -1]
    # train-test split
    X_train, X_test, y_train, y_test = train_test_split(X,
                                                        y,
                                                        test_size=0.25,
                                                        random_state=42)
    # build model at save it
    model = HGBR().fit(X_train, y_train)
    joblib.dump(model, f"../models/{file_name}.joblib")
    # print statistics (accuracy, precision and recall)
    y_pred = model.predict(X_test)
    accuracy = accuracy_score(y_test, y_pred)
    precision = precision_score(y_test, y_pred, average="macro")
    recall = recall_score(y_test, y_pred, average="macro")
    print(f"Metrics on testing data ('{file_name}'):\n"
          f"- {accuracy = :.2%}\n- {precision = :.2%}\n- {recall = :.2%}")
    
    
if __name__ == "__main__":
    
    Path("../models").mkdir(exist_ok=True)
    
    datasets = ("all",
                "dna_modification",
                "dna_replication",
                "lysis",
                "lysogeny_repressor",
                "packaging",
                "structural",
                "other")
    
    print("---")
    for fname in datasets:
        build_model_and_display_stats(file_name=fname)
        print("---")
        