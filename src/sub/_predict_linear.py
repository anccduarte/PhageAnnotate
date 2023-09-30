# -*- coding: utf-8 -*-

import joblib
import numpy as np
import pandas as pd
import sys
sys.path.append("..")
from Bio import SeqIO
from ml_dataset import MLDataset

def predict_function(path_to_file: str,
                     models_dir: str,
                     threshold: float) -> tuple:
    """
    Predicts the function associated to each sequence record present in
    <path_to_file>. A prediction is considered "valid" only if the class
    predicted by the model in <models_dir> has an associated probability
    equal to or larger than <threshold>. The result is saved to a .csv
    file in the parent directory of <models_dir>.
    
    Parameters
    ----------
    path_to_file: str
        The path to the file containing the sequences
    models_dir: str
        The relative path to the directory storing the model, scaler and
        support vector used in the predictive task
    threshold: float
        The threshold considered when establishing whether a prediction
        may be regarded as "valid"
    """
    # ---
    def _get_data(path_to_file: str) -> pd.DataFrame:
        # initilialize MLDataset object and call "build_dataset"
        data = MLDataset(file=path_to_file,
                         prot_name="unknown",
                         ttable=TTABLE,
                         icodons=ICODONS).build_dataset()
        # retrieve descriptions and feature vectors
        descriptions = data["Description"].tolist()
        X_to_pred = data.iloc[:, :-2]
        # return descriptions and feature vectors
        return descriptions, X_to_pred
    # ---
    def _get_model(models_dir: str) -> tuple:
        # load model, scaler and support vector
        model = joblib.load(f"{models_dir}/model-all.joblib")
        scaler = joblib.load(f"{models_dir}/scaler-all.joblib")
        support = joblib.load(f"{models_dir}/support-all.joblib")
        # return tuple containing the model, scaler and support vector
        return model, scaler, support
    # ---
    def _predict() -> tuple:
        # initialize list of predictions
        predictions = []
        # get descriptions and feature vectors
        descriptions, X_to_pred = _get_data(path_to_file=path_to_file)
        # load model, scaler and support vector
        model, scaler, support = _get_model(models_dir=models_dir)
        # main loop -> predict
        for _, row in X_to_pred.iterrows():
            # transform feature vector
            vec = np.array(row).reshape(1,-1)
            X_vec = scaler.transform(vec)[:, support]
            # predict and update list of predictions
            vld = np.amax(model.predict_proba(X_vec)) >= threshold
            func = model.predict(X_vec)[0] if vld else "other-function"
            predictions.append(func)
        # return list of predictions
        return descriptions, predictions
    # ---
    def save_to_csv() -> None:
        # get descriptions and predictions
        descriptions, predictions = _predict()
        # save predictions to .csv file
        data = pd.DataFrame(data={"Description": descriptions,
                                  "Function": predictions})
        data.to_csv(f"{RESULTS_DIR}/{NAME}.csv")
    # ---
    return save_to_csv()


if __name__ == "__main__":
    
    """
    Important note: For reasons of simplicity, renamed the files in
    "naive/models" to not include the suffix "-tol100". This padronizes
    the names of the files in "naive/models" and "linear/models".
    """
    
    import utils
    from pathlib import Path
    
    # read arguments (sys.argv)
    args = utils.get_args(("-path_to_file",),
                          ("-models_dir",),
                          ("-threshold",))
    path_to_file = args.path_to_file
    models_dir = args.models_dir
    threshold = float(args.threshold)
            
    # initialize TTABLE and ICODONS
    TTABLE = "11"
    ICODONS = ("TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG")
    
    # get name of the file
    NAME = path_to_file.split("/")[-1].split(".")[0]
    
    # create "results" directory (if not already created)
    RESULTS_DIR = "/".join(models_dir.split("/")[:-1]) + "/results"
    Path(RESULTS_DIR).mkdir(exist_ok=True)
    
    # predict and save results to .csv file
    predict_function(path_to_file, models_dir, threshold)
