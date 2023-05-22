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
    
    def __init__(self, path: str, ttable: str, icodons: tuple) -> None:
        """
        Initializes an instance of PredictFunction.
        
        Parameters
        ----------
        path: str
            The path to the .fasta file
        ttable: str
            The identifier of the translation table to be used
        icodons: tuple
            A tuple of possible initiation codons given the table <ttable>
            
        Attributes
        ----------
        _models: list
            A list of ML models
        """
        self.path = path
        self.ttable = ttable
        self.icodons = icodons
        self._models = PredictFunction._load_models()
        
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
                 "dna_modification",
                 "dna_replication",
                 "lysis",
                 "lysogeny_repressor",
                 "packaging",
                 "structural",
                 "other")
        out = []
        for model in models:
            out.append(joblib.load(f"../models/{model}.joblib"))
        return out

    def _get_dataset(self) -> pd.DataFrame:
        """
        Constructs and returns a featurized dataset from the sequences contained in
        the file pointed to by <path>.
        """
        data = MLDataset(file=self.path,
                         prot_name="unknown",
                         ttable=self.ttable,
                         icodons=self.icodons).build_dataset()
        return data.iloc[:, :-1]
    
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
        # load scaler used to scale features in "all.csv"
        scaler_all = joblib.load("../models/all_scaler.joblib")
        # get models
        ALL, MOD, REP, LYSIS, LYS_REP, PACK, STRUCT, OTHER = self._models
        # predict functional class
        preds_func_class = ALL.predict(scaler_all.transform(X))
        # initialize empty list to store predicted functions
        preds_func = []
        # iterate through rows of <X> and predict function
        for i, vec in X.iterrows():
            # get appropriate model and load scaler
            func_class = preds_func_class[i]
            if func_class == "dna-modification":
                model = MOD
                scaler = joblib.load("../models/dna_modification_scaler.joblib")
            elif func_class == "dna-replication":
                model = REP
                scaler = joblib.load("../models/dna_replication_scaler.joblib")
            elif func_class == "lysis":
                model = LYSIS
                scaler = joblib.load("../models/lysis_scaler.joblib")
            elif func_class == "lysogeny-repressor":
                model = LYS_REP
                scaler = joblib.load("../models/lysogeny_repressor_scaler.joblib")
            elif func_class == "packaging":
                model = PACK
                scaler = joblib.load("../models/packaging_scaler.joblib")
            elif func_class == "structural":
                model = STRUCT
                scaler = joblib.load("../models/structural_scaler.joblib")
            elif func_class == "other":
                model = OTHER
                scaler = joblib.load("../models/other_scaler.joblib")
            # predict and save
            vec = np.array(vec).reshape(1,-1)
            preds_func.append(model.predict(scaler.transform(vec))[0])
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
