# -*- coding: utf-8 -*-

# basic imports
import collections
import joblib
import logging
import numpy as np
import os
import pandas as pd
import sys
import warnings
# sklearn estimators
from sklearn.ensemble import HistGradientBoostingClassifier as HGBC
from sklearn.ensemble import RandomForestClassifier as RFC
from sklearn.naive_bayes import GaussianNB as GNB
from sklearn.neighbors import KNeighborsClassifier as KNC
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier as DTC
# MLPClassifier wrapper (works fine with BayesSearchCV -> temporary?)
from utils import MLPCWrapper as MLPC
# other sklearn
from sklearn.base import BaseEstimator
from sklearn.feature_selection import SelectFromModel
from sklearn.metrics import accuracy_score, f1_score
from sklearn.metrics import precision_score, recall_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler
# skopt
from skopt import BayesSearchCV
from skopt.space import Categorical, Integer, Real

# ignore sklearn warnings
warnings.filterwarnings("ignore")

class MLModel:
    
    """
    Constructs a machine learning (ML) model based on an ML algorithm and data provided
    by the user. The validation and test sizes are also adjustable.
    """
    
    ALGOS = {"gradient-boosting": HGBC,
             "random-forest": RFC,
             "naive-bayes": GNB,
             "k-nearest-neighbors": KNC,
             "neural-network": MLPC,
             "support-vector-machine": SVC,
             "decision-tree": DTC}
    
    HYPER = {"gradient-boosting": {"learning_rate": Real(0.05, 0.2), # [0.05, 0.1, 0.2]
                                   "max_iter": Integer(100, 300), # [100, 150, 200]
                                   "max_leaf_nodes": Integer(15, 63), # [15, 31, 63]
                                   "min_samples_leaf": Integer(10, 30)}, # [10, 20, 30]
             # ---
             "random-forest": {"n_estimators": Integer(100, 300), # [100, 150, 200]
                               "criterion": Categorical(["gini", "entropy"]),
                               "min_samples_split": Integer(2, 5), # [2, 3, 4]
                               "min_samples_leaf": Integer(1, 3), # [1, 2]
                               "max_features": Categorical(["sqrt", "log2"])},
             # ---
             "naive-bayes": {}, # no optimization
             # ---
             "k-nearest-neighbors": {"n_neighbors": Integer(3, 8),
                                     "algorithm": Categorical(["ball_tree", # stochastic
                                                               "kd_tree"]), # stochastic
                                     "p": Integer(1, 2)}, # p -> minkowski power
             # ---
             "neural-network": {"num_layers": Integer(1, 6),
                                "activation": Categorical(["relu", "tanh"]),
                                "solver": Categorical(["adam", "sgd"]),
                                "learning_rate": Categorical(["constant", "invscaling"]),
                                "learning_rate_init": Real(0.001, 0.01), # [0.001, 0.005]
                                "max_iter": Integer(200, 2000)}, # [500, 1000, 2000]
             # ---
             "support-vector-machine": {"C": Real(0.1, 10), # [0.1, 1, 10]
                                        "kernel": Categorical(["linear", "poly", "rbf"]),
                                        "degree": Integer(2, 6), # [2, 3, 4]
                                        "decision_function_shape": Categorical(["ovo",
                                                                                "ovr"]),
                                        "probability": Categorical([True])},
             # ---
             "decision-tree": {"criterion": Categorical(["gini", "entropy"]),
                               "min_samples_split": Integer(2, 5), # [2, 3, 4]
                               "min_samples_leaf": Integer(1, 3)}} # [1, 2]
    
    def __init__(self,
                 data_path: str,
                 models_dir: str,
                 algorithm: str,
                 test_size: float,
                 random_state: int = 0,
                 init: bool = True,
                 final_model: bool = False) -> None:
        """
        Initializes an instance of MLModel.
        
        Parameters
        ----------
        data_path: str
            The path to the .csv file containing the data feeding the model
        models_dir: str
            The name of the directory where the .joblib files storing information
            on the constructed model are to be saved (only relevant if <final_model>
            is set to True)
        algorithm: str
            The name of the ML algorithm to be used
        test_size: float
            The proportion of the data to be used for testing (only relevant if
            <init> is set to False - if <init> is set to True, the train-test split
            is carried out via loading the respective sets from memory)
        random_state: int (default=0)
            Controls the split of dataset into training, testing and validation
        init: bool (default=True)
            Whether the model being built belongs to the first set of models (note
            that the second set of models are fed on datasets which encompass
            predictions made by the first set of models)
        final_model: bool (default=False)
            Whether the model to be trained is definitive. If set to True, the
            scaler, the support vector (used for feature selection) and the model
            itself are saved in .joblib files for further use (make predictions on
            new data). Moreover, a .txt file storing the indices of the rows that
            comprise the testing set is created if <init> is also set to True
            
        Attributes
        ----------
        _dataset: pd.DataFrame
            The dataset containing the data that will feed the model
        _name_set: str
            The name of the dataset (extension stripped off)
        _estimator: BaseEstimator
            The sklearn estimator to be used when constructing the ML model
        _logger: logging.Logger
            A logging.Logger instance
        """
        # parameters
        self.data_path = data_path
        self.models_dir = models_dir
        self.algorithm = algorithm
        self.test_size = test_size
        self.random_state = random_state
        self.init = init
        self.final_model = final_model
        # attributes
        self._dataset = MLModel._read_csv(data_path)
        self._name_set = data_path.split("/")[-1].split(".")[0].replace("_","-")
        self._estimator = MLModel.ALGOS[algorithm]
        self._logger = self._get_logger()
    
    @staticmethod
    def _read_csv(data_path: str) -> pd.DataFrame:
        """
        Reads and returns the dataset pointed to by <data_path>. The contents of all
        columns except for the last two are casted to np.float32, reducing RAM usage
        by approximately 50%. Returns the newly created dataset.
        
        Parameters
        ----------
        data_path: str
            The path to the .csv file containing the data feeding the model
        """
        # construct mapping for data types
        map_ = collections.defaultdict(lambda: np.float32)
        map_["Description"] = object
        map_["Function"] = object
        # read .csv (with types specified by <map_>) and return the resulting dataset
        return pd.read_csv(data_path, dtype=map_)
    
    def _get_logger(self):
        """
        Creates a logging.Logger composed of two distinct handlers: a stream handler
        and a file handler. This means that every log message has as its destination
        both the standard output stream (prints messages to screen) and a file.
        """
        # initialize logger instance and set level
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.INFO)
        # remove handlers from logger (if existent)
        logger.handlers.clear()
        # add stream handler to the logger
        logger.addHandler(logging.StreamHandler(sys.stdout))
        # add file handler to the logger
        nfile = f"{self.models_dir}/info-{self._name_set}-{self.algorithm}.log"
        logger.addHandler(logging.FileHandler(nfile))
        # return the logger instance
        return logger
    
    def _log(self, message: str) -> None:
        """
        Wrapper for "self._logger". Logs with level "info".
        
        Parameters
        ----------
        message: str
            The message to be logged
        """
        self._logger.info(message)
    
    def __repr__(self) -> str:
        """
        Returns the string representation of the object.
        """
        c = self.__class__.__name__
        d, m, a, t = self.data_path, self.models_dir, self.algorithm, self.test_size
        r, i, f = self.random_state, self.init, self.final_model
        res = f"{c}({d!r}, {m!r}, {a!r}, {t!r}, {r!r}, {i!r}, {f!r})"
        return res
            
    @staticmethod
    def _save_test_data_info(x_test: pd.DataFrame,
                             dl_test: pd.DataFrame,
                             file_name: str) -> None:
        """
        Saves a .txt file containing the labels and descriptions associated to the
        rows assigned to testing data, and a .csv file containing the testing data
        itself. It is only called when <self.final_model> and <self.init> are set to
        True. This allows for testing the second set of models (the ones fed on an
        expanded dataset encompassing predictions made by the first set of models) on
        the same data used to test the first set of models.
        
        Parameters
        ----------
        x_test: pd.DataFrame
            The feature vectors of the testing portion of the data
        dl_test: pd.DataFrame
            The labels vector of the testing portion of the data, and the descriptions
            associated to those labels
        file_name: str
            The name to be given to the files storing information on testing data
        """
        # save sequence labels and descriptions to a .txt file
        with open(file_name+".txt", "w") as fout:
            for _, row in dl_test.iterrows():
                label, description = row["Function"], row["Description"]
                fout.write(f"- {label}: {description}\n")
        # save test set to a .csv file (to be loaded in "_train_test_split_cs")
        test_data = x_test.copy(deep=True)
        test_data["Function"] = dl_test.iloc[:, -1]
        test_data.to_csv(file_name+".csv")
    
    def _train_test_split_init(self) -> tuple:
        """
        Executes the division of the dataset into training and testing. Returns a
        tuple containing 2 pd.DataFrame objects (corresponding to the training and
        testing feature vectors) and 2 pd.Series objects (corresponding to the
        training and testing label vectors). It executes the train-test split by
        calling sklearn's function "train_test_split". The method is employed when
        <self.init> is set to True.
        """
        # log state of the process (train-test split)
        self._log("Performing 'init' train-test split on the dataset...")
        # split the dataset into features and labels
        df_feats = self._dataset.iloc[:, 1:-2] # remove "Description" and "Function"
        df_other = self._dataset.iloc[:, -2:] # keep "Description" and "Function"
        # split the data into training and testing sets
        x_trn, x_tst, dl_trn, dl_tst = train_test_split(df_feats,
                                                        df_other,
                                                        stratify=df_other.iloc[:, -1],
                                                        test_size=self.test_size,
                                                        random_state=self.random_state)
        # if <self.final_model> and <self.init> are set to True, save test set info
        if self.final_model:
            fname = self.models_dir + "/test-data-" + self._name_set
            MLModel._save_test_data_info(x_tst, dl_tst, fname)
        # return tuple of pandas DataFrames/Series
        # (.iloc[:, -1] to discard the column "Description" from dl_trn and dl_tst)
        return x_trn, x_tst, dl_trn.iloc[:, -1], dl_tst.iloc[:, -1]
        
    def _train_test_split_cs(self) -> tuple:
        """
        Executes the division of the dataset into training and testing. Returns a
        tuple containing 2 pd.DataFrame objects (corresponding to the training and
        testing feature vectors) and 2 pd.Series objects (corresponding to the
        training and testing label vectors). It executes the train-test split by
        loading "init" test and trainign sets. The method is employed when <self.init>
        is set to False.
        """
        # log state of the process (train-test split)
        self._log("Performing 'cs' train-test split on the dataset...")
        # read test and train sets from .csv files
        # (assumes that the training set is already processed, that is, that redundant
        # sequences in the training data relative to the sequences in the testing data
        # were previously removed using the software CD-HIT -> cd-hit-est-2d)
        train_set = self._dataset
        test_set = pd.read_csv(f"../models/test-data-{self._name_set}.csv")
        # split features and labels (note that "Description" is only present in the
        # training data (*) -> see "_train_test_split_init")
        x_trn, y_trn = train_set.iloc[:, 1:-2], train_set.iloc[:, -1] # (*) 1:-2
        x_tst, y_tst = test_set.iloc[:, 1:-1], test_set.iloc[:, -1] # (*) 1:-1
        # return tuple of pandas DataFrames/Series
        return x_trn, x_tst, y_trn, y_tst
        
    def _scale_data(self, x_train: np.ndarray, x_test: np.ndarray) -> tuple:
        """
        Returns a scaled version of the training and testing data. To scale the data,
        the training set is used to fit a min-max scaler, and the resulting object is
        used to transform the original training and testing sets.
        
        Parameters
        ----------
        x_train: np.ndarray
            The feature vectors of the training portion of the data
        x_test: np.ndarray
            The feature vectors of the testing portion of the data
        """
        # log state of the process (scale data)
        self._log("Scaling training and testing data...")
        # initialize MinMaxScaler and fit it to <x_train>
        scaler = MinMaxScaler().fit(x_train)
        # save scaler to a .joblib file if <final_model> is set to True
        if self.final_model:
            scaler_name = self.models_dir + "/scaler-" + self._name_set
            joblib.dump(scaler, scaler_name+".joblib")
        # return <x_train> and <x_test> scaled
        return scaler.transform(x_train), scaler.transform(x_test)
    
    def _select_features(self,
                         x_train: np.ndarray,
                         y_train: np.ndarray,
                         x_test: np.ndarray) -> tuple:
        """
        Returns shrunken versions of the training and testing sets by selecting the
        more important features of the dataset. This is achieved by creating an
        instance of sklearn's SelectFromModel and fitting it on a fitted sklearn's
        RandomForestClassifier (RFC) estimator (the estimator is fitted on the
        training data, <x_train>). The selector is specifically fitted on a RFC
        estimator, since this particular estimator initializes a mandatory attribute
        "feature_importances_" which is necessary for fitting SelectFromModel. The
        selector works by keeping the features whose absolute importance value is
        greater than or equal to a specified threshold (the mean of the importances
        of all features) and discarding the remaining features.
        
        Parameters
        ----------
        x_train: np.ndarray
            The feature vectors of training portion of the data
        y_train: np.ndarray
            The label vector of the training portion of the data
        x_test: np.ndarray
            The feature vectors of the testing portion of the data
        """
        # log state of the process (select features)
        self._log("Selecting features...")
        # compute and save support vector for feature selection if non existent;
        # otherwise, load it from the appropriate .joblib file
        support_name = self.models_dir + "/support-" + self._name_set
        if os.path.exists(support_name+".joblib"):
            support_vector = joblib.load(support_name+".joblib")
        else:
            selector = SelectFromModel(estimator=RFC(n_estimators=200,
                                                     criterion="entropy"))
            selector.fit(x_train, y_train)
            support_vector = selector.get_support()
            joblib.dump(support_vector, support_name+".joblib")
        # display number of selected features
        self._log(f"- num_selected = {len(support_vector[support_vector==True])}")
        # return shrunken versions of <x_train> and <x_test>
        return x_train[:, support_vector], x_test[:, support_vector]
    
    def _optimize_hyperparameters(self,
                                  x_train: np.ndarray,
                                  y_train: np.ndarray) -> dict:
        """
        Computes the best combination of hyperparameters, for a given estimator and
        training set, by fitting a BayesSearchCV object. Returns a dictionary object
        whose keys are the names of the hyperparameters and whose values are the best
        value found for the respective hyperparameter.
        
        Parameters
        ----------
        x_train: np.ndarray
            The feature vectors of training portion of the data
        y_train: np.ndarray
            The label vector of the training portion of the data
        """
        # log state of the process (optimize hyperparameters)
        self._log("Optimizing hyperparameters...")
        # select hyperparameter grid based on <self.algorithm>
        search_spaces = MLModel.HYPER[self.algorithm] # param_grid
        # skip optimization if <search_spaces> is empty (avoids unnecessary fit)
        if not search_spaces: # param_grid
            best_params = {}
        # otherwise, determine optimal combination of hyperparameters
        else:
            # the search's exhaustiveness and robustness depends on <self.final_model>
            num_iter = 32 if self.final_model else 16
            cv_folds = 5 if self.final_model else 3
            # (cv: "For integer/None inputs, if the estimator is a classifier and y is
            # either binary or multiclass, StratifiedKFold is used. In all other cases,
            # KFold is used. These splitters are instantiated with shuffle=False so the
            # splits will be the same across calls.")
            opt = BayesSearchCV(estimator=self._estimator(), # GridSearchCV
                                search_spaces=search_spaces, # param_grid
                                n_iter=num_iter,
                                scoring="f1_macro",
                                cv=cv_folds,
                                refit=False).fit(x_train, y_train)
            # assign <opt.best_params_> to <best_params> to match the first case
            # (cast to dict -> <opt.best_params_> is an OrderedDict)
            best_params = dict(opt.best_params_)
        # display and return best combination of hyperparameters
        self._log(f"- {best_params = }")
        return best_params
    
    def _test_model(self,
                    model: BaseEstimator,
                    x_train: np.ndarray,
                    x_test: np.ndarray,
                    y_test: np.ndarray) -> None:
        """
        Tests the model <model> on the testing data <x_test> and <y_test>. Computes
        the accuracy score, precision score, recall score and f1-score of the model.
        
        Parameters
        ----------
        model: BaseEstimator
            An estimator fitted on training data
        x_train: np.ndarray
            The feature vectors of training portion of the data (in this context, it
            is only needed for determining the number of examples in the training set)
        x_test: np.ndarray
            The feature vectors of testing portion of the data
        y_test: np.ndarray
            The label vector of the testing portion of the data
        """
        # get percentage of data used for testing
        rows_trn, rows_tst = x_train.shape[0], x_test.shape[0]
        perc_test = rows_tst / (rows_trn + rows_tst)
        # log state of the process (test model)
        self._log(f"\nTesting {self.algorithm.upper()} model ({self._name_set})")
        self._log(f"('x_test' corresponds to {perc_test:.2%} of the data)\n---")
        # compute predictions
        y_pred = model.predict(x_test)
        # compute metrics and display them on screen
        accuracy = accuracy_score(y_test, y_pred)
        precision = precision_score(y_test, y_pred, average="macro")
        recall = recall_score(y_test, y_pred, average="macro")
        f1_scr = f1_score(y_test, y_pred, average="macro")
        self._log(f"Metrics on testing data:\n{accuracy = :.2%}\n{precision = :.2%}"
                  f"\n{recall = :.2%}\nf1-score = {f1_scr:.2%}\n\n")
    
    def build_model(self) -> None:
        """
        Builds an ML model given the algorithm and the data provided by the user.
        To do so, splits the dataset in training and testing sets, scales them using
        a MinMaxScaler, selects the most important features by fitting a sklearn's
        SelectFromModel instance on a RandomForestClassifier (fitted on the training
        data) and optimizes the hyperparameters of the estimator chosen by the user.
        The model is then fitted on the training data using these hyperparameters.
        Finally, the model is tested on the testing data. If <self.final_model> is
        set to True, the model is refitted to the whole dataset (train + test) and
        saved for further use.
        """
        # log state of the process (build model)
        algo = self.algorithm.upper()
        self._log(f"\nBuilding {algo} model ({self._name_set})\n---")
        # train-test split (<self.init> dictates which method is called)
        if self.init:
            x_train, x_test, y_train, y_test = self._train_test_split_init()
        else:
            x_train, x_test, y_train, y_test = self._train_test_split_cs()
        # scale <x_train> and <x_test>
        x_train, x_test = self._scale_data(x_train, x_test)
        # select most important features
        x_train, x_test = self._select_features(x_train, y_train, x_test)
        # optimize hyperparameters and fit the model based on them
        params = self._optimize_hyperparameters(x_train, y_train)
        # fit model on the best combination of hyperparamters
        self._log("Fitting estimator on best combination of hyperparameters...")
        model = self._estimator(**params).fit(x_train, y_train)
        self._log("DONE")
        # test the fitted estimator (model) on the testing data
        self._test_model(model, x_train, x_test, y_test)
        # save model to a .joblib file if <final_model>
        if self.final_model:
            model_name = self.models_dir + "/model-" + self._name_set
            joblib.dump(model, model_name+".joblib")
