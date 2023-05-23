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
    
    # initialize class variable (all instances of PredictFunction share it)
    DATASETS = ("all", "dna-modification", "dna-replication", "lysis",
                "lysogeny-repressor", "packaging", "structural", "other")
    
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
        _models: dict
            A dictionary of ML models
        _scalers: dict
            A dictionary of scaler objects
        _support_vecs: dict
            A dictionary of support vectors (feature selection)
        """
        # parameters
        self.path = path
        self.ttable = ttable
        self.icodons = icodons
        # attributes
        DATASETS = PredictFunction.DATASETS
        models, scalers, support_vecs = PredictFunction._load()
        self._models = {ds: mo for (ds, mo) in zip(DATASETS, models)}
        self._scalers = {ds: sc for (ds, sc) in zip(DATASETS, scalers)}
        self._support_vecs = {ds: sp for (ds, sp) in zip(DATASETS, support_vecs)}
        
    def __repr__(self) -> str:
        """
        Returns the string representation of the object.
        """
        class_ = self.__class__.__name__
        return f"{class_}({self.path!r})"
        
    @staticmethod
    def _load() -> tuple:
        """
        Loads and returns the models, scalers and support vectors corresponding to
        each dataset in "DATASETS".
        """
        models, scalers, support_vecs = [], [], []
        for ds in PredictFunction.DATASETS:
            name_ds = "_".join(ds.split("-"))
            models.append(joblib.load(f"../models/{name_ds}.joblib"))
            scalers.append(joblib.load(f"../models/{name_ds}_scaler.joblib"))
            support_vecs.append(joblib.load(f"../models/{name_ds}_support.joblib"))
        return models, scalers, support_vecs

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
        in <X_pred> (each representing a DNA sequence). Returns a tuple of two lists
        containing the predictions.
        """
        # get featurized dataset
        X_pred = self._get_dataset()
        # load model, scaler and support vector for "all.csv"
        model_all = self._models["all"]
        scaler_all = self._scalers["all"]
        support_all = self._support_vecs["all"]
        # predict functional class
        X_all = scaler_all.transform(X_pred)[:, support_all]
        preds_func_class = model_all.predict(X_all)
        # iterate through rows of <X>, predict function and save prediction
        preds_func = []
        for i, vec in X_pred.iterrows():
            # get/load appropriate model, scaler and feature suppport
            func_class = preds_func_class[i]
            model = self._models[func_class]
            scaler = self._scalers[func_class]
            support = self._support_vecs[func_class]
            # predict and save
            vec = np.array(vec).reshape(1,-1)
            X_vec = scaler.transform(vec)[:, support]
            preds_func.append(model.predict(X_vec)[0])
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
