# -*- coding: utf-8 -*-

import joblib
import numpy as np
import os
import pandas as pd
from pathlib import Path
from sklearn.base import BaseEstimator
from sklearn.ensemble import RandomForestClassifier as RFC
from sklearn.feature_selection import SelectFromModel
from sklearn.metrics import accuracy_score, precision_score, recall_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler

# create directory (if not already created)
MODELS = "../models"
Path(MODELS).mkdir(exist_ok=True)

def get_scaler(dataset: str, X_train: np.ndarray) -> MinMaxScaler:
    """
    If the process of feature scaling was not previously performed for the
    dataset <dataset>, the scaler is computed. Otherwise, the scaler is loaded
    from a .joblib file. Returns the scaler object.
    
    Parameters
    ----------
    dataset: str
        The name of the dataset
    X_train: np.ndarray
        The scaled training data
    """
    scaler_name = f"{MODELS}/{dataset}_scaler.joblib"
    if not os.path.exists(scaler_name):
        scaler = MinMaxScaler().fit(X_train)
        # only saved for use in this module (saves some time)
        joblib.dump(scaler, scaler_name)
    else:
        scaler = joblib.load(scaler_name)
    return scaler

def get_support_vector(dataset: str,
                       X_train: np.ndarray,
                       y_train: np.ndarray) -> np.ndarray:
    """
    If the process of feature selection was not previously performed for the
    dataset <dataset>, the features are selected and the respective support
    vector is computed. Otherwise, the support vector is loaded from a .joblib
    file. Returns the support vector.
    
    Parameters
    ----------
    dataset: str
        The name of the dataset
    X_train: np.ndarray
        The scaled training data
    y_train: np.ndarray
        The label vector (response variable)
    """
    support_name = f"{MODELS}/{dataset}_support.joblib"
    if not os.path.exists(support_name):
        selector = SelectFromModel(estimator=RFC()).fit(X_train, y_train)
        support = selector.get_support()
        # only saved for use in this module (saves a lot of time)
        joblib.dump(support, support_name)
    else:
        support = joblib.load(support_name)
    return support

def build_model_and_display_stats(dataset: str, model: BaseEstimator) -> None:
    """
    Builds an ML model fed on the data present in 'database/<file_name>.csv',
    saves it to a .joblib file (TEMPORARY) and prints some useful statisitics
    on the testing data (accuracy, precision and recall).
    
    Parameters
    ----------
    dataset: str
        The name of the .csv file containing the data
    """
    # ---
    # NOTE: THIS WILL EXCLUSIVELY BE USED TO TRAIN A MODEL WITH 75% OF THE DATA
    # AND TEST ON THE REMAINING 25% (HENCE, THE MODELS, SCALERS AND SEELCTORS
    # WILL NOT BE SAVED HERE...). THE ONLY PURPOSE OF THIS FUNCTION IS TO TRAIN
    # AN ML MODEL, TEST IT AND DISPLAY THE RESULTS. THIS INCLUDES: SPLIT THE
    # DATA INTO TRAIN AND TEST, SCALE THE DATA BASED ON THE "CHARACTERISTICS"
    # OF THE TRAINING DATA, SELECT FEATURES (ALSO, BASED ON THE TRAINING DATA),
    # TRAIN THE MODEL, TEST THE MODEL, AND DISPLAY THE RESULTS OF TESTING. ONLY
    # THE BEST MODELS (ML ALGORITHM + DATASET) WILL BE FULLY TRAINED IN A
    # SEPARATE FUNCTION / CLASS (ONLY THEN THE MODELS, SCALERS AND SEELCTORS
    # WILL BE SAVED). AS OF NOW, THE MODELS AND SCALERS ARE SAVED HERE FOR
    # PRACTICALITY REASONS.
    # ---
    # read data
    print(f"{dataset.upper()} -> {model().__class__.__name__}")
    print("Reading and splitting data...")
    data = pd.read_csv(f"../database/{dataset}.csv")
    # ---
    # split features (explanatory) and label (response)
    X, y = data.iloc[:, 1:-1], data.iloc[:, -1]
    # ---
    # train-test split
    X_train, X_test, y_train, y_test = train_test_split(X,
                                                        y,
                                                        test_size=0.25,
                                                        random_state=42)
    # ---
    # scale data (X_train and X_test) according to X_train
    print("Scaling data...")
    scaler = get_scaler(dataset, X_train)
    X_train, X_test = scaler.transform(X_train), scaler.transform(X_test)
    # ---
    # select features (based on the Random Forest Classifier...)
    print("Selecting features...")
    support = get_support_vector(dataset, X_train, y_train)
    X_train, X_test = X_train[:, support], X_test[:, support]
    # ---
    # build model
    print("Building model on 75% of the data...")
    estimator = model().fit(X_train, y_train)
    # ---
    # TEMPORARY -> save model (only for practicality reasons)
    joblib.dump(estimator, f"{MODELS}/{dataset}.joblib")
    # ---
    # print statistics (accuracy, precision and recall)
    print("Testing model on 25% of the data...")
    y_pred = estimator.predict(X_test)
    accuracy = accuracy_score(y_test, y_pred)
    precision = precision_score(y_test, y_pred, average="macro")
    recall = recall_score(y_test, y_pred, average="macro")
    print(f"Metrics on testing data ({dataset!r}):\n"
          f"- {accuracy = :.2%}\n- {precision = :.2%}\n- {recall = :.2%}")
    
    
if __name__ == "__main__":
    
    from sklearn.ensemble import HistGradientBoostingClassifier as HGBR
    # ---
    # from sklearn.naive_bayes import CategoricalNB as NaiveBayes
    # from sklearn.neighbors import KNeighborsClassifier as KNC
    # from sklearn.neural_network import MLPClassifier as MLPC
    # from sklearn.svm import SVC
    # from sklearn.tree import DecisionTreeClassifier as DTC
    # ---
    
    Path("../models").mkdir(exist_ok=True)
    
    # tuple of datasets that will feed the ML models
    datasets = ("all",
                "dna_modification",
                "dna_replication",
                "lysis",
                "lysogeny_repressor",
                "packaging",
                "structural",
                "other")
    
    # tuple of ML algorithms used to build the models
    algorithms = (HGBR,) # (HGBR, NaiveBayes, KNC, DTC, RF, SVC, MLPC)
    
    print("---")
    for dataset in datasets:
        for algo in algorithms:
            build_model_and_display_stats(dataset=dataset, model=algo)
            print("---")
        