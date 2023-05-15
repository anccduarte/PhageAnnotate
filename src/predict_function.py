# -*- coding: utf-8 -*-

import joblib
import numpy as np
import pandas as pd
from ml_dataset import MLDataset
from pathlib import Path

class PredictFunction:
    
    """
    Allows to predict the function of the proteins coded by the DNA sequences
    contained in the .fasta file pointed to by <path>.
    """
    
    def __init__(self, path: str) -> None:
        """
        Initializes an instance of PredictFunction.
        
        Parameters
        ----------
        path: str
            The path to the .fasta file
            
        Attributes
        ----------
        _models: list
            A list of ML models
        """
        self.path = path
        self._models = self._load_models()
        
    def __repr__(self) -> str:
        """
        Returns the string representation of the object.
        """
        class_ = self.__class__.__name__
        return f"{class_}({self.path!r})"
        
    @staticmethod
    def _load_models() -> list:
        """
        Loads and returns the HGBR models in 'models'.
        """
        models = ("all",
                  "lysis",
                  "modification_replication",
                  "packaging",
                  "structural")
        out = []
        for model in models:
            out.append(joblib.load(f"../models/{model}.joblib"))
        return out

    def _get_dataset(self) -> pd.DataFrame:
        """
        Constructs and returns a featurized dataset from the sequences contained in
        the file pointed to by <path>.
        """
        data = MLDataset(file=self.path, protein_name="unknown")
        X = data.build_dataset().iloc[:, :-1]
        return X
    
    def _get_descriptions(self) -> list:
        """
        Parses the .fasta file pointed to by <path> and saves the descriptions of the
        sequences in a list.
        """
        descrips = []
        with open(self.path) as handle:
            lines = handle.readlines()
        for line in lines:
            if line.startswith(">"):
                descrip = line[2:-1]
                descrips.append(descrip)
        return descrips

    def _predict(self) -> tuple:
        """
        Predicts the functional class and function associated to each feature vector
        in <X> (each representing a DNA sequence). Returns a tuple of two lists
        containing the predictions.
        """
        # get featurized dataset
        X = self._get_dataset()
        # get models
        MODEL_ALL, MODEL_LYSIS, MODEL_MOD_REP, MODEL_PACK, MODEL_STRUCT = self._models
        # predict functional class
        preds_func_class = MODEL_ALL.predict(X)
        # predict function
        preds_func = []
        # iterate through rows of <X>
        for i, vec in X.iterrows():
            # get appropriate model
            func_class = preds_func_class[i]
            if func_class == "lysis": model = MODEL_LYSIS
            elif func_class == "modification-replication": model = MODEL_MOD_REP
            elif func_class == "packaging": model = MODEL_PACK
            elif func_class == "structural": model = MODEL_STRUCT
            # predict and save
            vec = np.array(vec).reshape(1,-1)
            preds_func.append(model.predict(vec)[0])
        # return predictions
        return preds_func_class, preds_func

    def save_results(self, name: str) -> None:
        """
        Saves results (predictions) to a .csv file.

        Parameters
        ----------
        name: str
            The name to be given to the .csv file
        """
        # get descriptions in the .fasta file
        descrips = self._get_descriptions()
        # get predictions
        func_class, func = self._predict()
        # construct dataframe
        df = pd.DataFrame(data={"Description": descrips,
                                "Functional-Class": func_class,
                                "Function": func})
        # save to .csv file
        df.to_csv(f"../results/{name}.csv")
