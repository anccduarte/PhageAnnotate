# -*- coding: utf-8 -*-

import joblib
import numpy as np
import pandas as pd
from ml_dataset import MLDataset

# REMOVE in later stage of development
# (used in "_predict_func_class" and "_predict_function")
from tqdm import tqdm

class PredictFunction:
    
    """
    Allows to predict the function of the proteins coded by the DNA sequences
    contained in the .fasta file pointed to by <path>.
    """
    
    # initialize class variable (all instances of PredictFunction share it)
    DATASETS = ("all", "dna-modification", "dna-replication", "lysis",
                "lysogeny-repressor", "packaging", "structural")
    
    def __init__(self,
                 path: str,
                 models_dir: str,
                 thresholds: list,
                 ttable: str,
                 icodons: tuple) -> None:
        """
        Initializes an instance of PredictFunction.
        
        Parameters
        ----------
        path: str
            The path to the .fasta file
        models_dir: str
            The name of the directory containing the models
        thresholds: list
            A list of two thresholds: the first relative to the functional class, the
            second relative to the specific function
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
        self.models_dir = models_dir
        self.thresholds = thresholds
        self.ttable = ttable
        self.icodons = icodons
        # attributes
        DATASETS = PredictFunction.DATASETS
        models, scalers, support_vecs = self._load() #PredictFunction
        self._models = {ds: mo for (ds, mo) in zip(DATASETS, models)}
        self._scalers = {ds: sc for (ds, sc) in zip(DATASETS, scalers)}
        self._support_vecs = {ds: sp for (ds, sp) in zip(DATASETS, support_vecs)}
        
    def __repr__(self) -> str:
        """
        Returns the string representation of the object.
        """
        class_ = self.__class__.__name__
        return f"{class_}({self.path!r})"
        
    def _load(self) -> tuple:
        """
        Loads and returns the models, scalers and support vectors corresponding to
        each dataset in "DATASETS".
        """
        models, scalers, support_vecs = [], [], []
        for ds in PredictFunction.DATASETS:
            # get model, scaler and support for <ds>
            model = joblib.load(f"../{self.models_dir}/model-{name_ds}.joblib")
            scaler = joblib.load(f"../{self.models_dir}/scaler-{name_ds}.joblib")
            support = joblib.load(f"../{self.models_dir}/support-{name_ds}.joblib")
            # append the latter to the respective lists
            models.append(model)
            scalers.append(scaler)
            support_vecs.append(support)
        return models, scalers, support_vecs
    
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
        
    def _predict_all(self,
                     vec: pd.Series,
                     model: str,
                     default: str,
                     threshold: float) -> str:
        """
        Labels the feature vector <vec> using the model <model>. If the vector of
        probabilities "prob_vec" is not valid (max(prob_vec) < threshold), the label
        of <vec> is set to <default>.
        
        Parameters
        ----------
        vec: pd.Series
            The feature vector
        model: str
            The model to be used when making the prediction on <vec>
        default: str
            The default label to be attributed to <vec> if the probability vector
            generated by predicting the class of <vec> is not valid
        threshold: float
            The probability threshold for the model to be able to confidently assign
            a label to a particular data instance
        """
        # load appropriate model, scaler and feature suppport
        estimator = self._models[model]
        scaler = self._scalers[model]
        support = self._support_vecs[model]
        # scale and select features of the feature vector
        vec = np.array(vec).reshape(1,-1)
        X_vec = scaler.transform(vec)[:, support]
        # evaluate "prob_vec" -> if not valid, function is set to <default>
        #                     -> otherwise, predict function using <estimator>
        prob_vec = estimator.predict_proba(X_vec)
        if np.amax(prob_vec) < threshold:
            function = default
        else:
            function = estimator.predict(X_vec)[0]
        # return the prediction
        return function
    
    def _predict_func_class(self,
                            X_pred: pd.DataFrame,
                            threshold: float) -> list:
        """
        Predicts the functional class associated to each feature vector present in
        the pd.DataFrame <X_pred>. Returns a list object containing the predictions.
        
        Parameters
        ----------
        X_pred: pd.DataFrame
            The dataframe containing the feature vectors
        threshold: float
            The probability threshold for the model to be able to confidently assign
            a label to a particular data instance
        """
        func_class = []
        for _, vec in tqdm(X_pred.iterrows()):
            pred = self._predict_all(vec, "all", "other", threshold)
            func_class.append(pred)
        return func_class
    
    def _predict_function(self,
                          X_pred: pd.DataFrame,
                          func_class: list,
                          threshold: float) -> list:
        """
        Predicts the function associated to each feature vector present in the
        pd.DataFrame <X_pred> based on the respective functional class given by the
        list object <func_class>. Returns a list containing the predictions.
        
        Parameters
        ----------
        X_pred: pd.DataFrame
            The dataframe containing the feature vectors
        func_class: list
            A list containing functional class predictions
        threshold: float
            The probability threshold for the model to be able to confidently assign
            a label to a particular data instance
        """
        preds_func = []
        for i, vec in tqdm(X_pred.iterrows()):
            class_ = func_class[i]
            if class_ == "other":
                func = ""
            else:
                func = self._predict_all(vec, class_, "other-function", threshold)
            preds_func.append(func)
        return preds_func
        
    def _predict(self) -> tuple:
        """
        Predicts the functional class and function associated to each feature vector
        in <X_pred> (each representing a DNA sequence). Returns a tuple of two lists
        containing the predictions.
        """
        # get featurized dataset and probability thresholds
        X_pred = self._get_dataset()
        thresh1, thresh2 = self.thresholds
        # predict functional class
        func_class = self._predict_func_class(X_pred, thresh1)
        # predict function based on the functional class of each protein
        preds_func = self._predict_function(X_pred, func_class, thresh2)
        # return predictions
        return func_class, preds_func

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
