# -*- coding: utf-8 -*-

# basic imports
import collections
import joblib
import numpy as np
import os
import pandas as pd
import warnings
# sklearn estimators
from sklearn.ensemble import HistGradientBoostingClassifier as HGBC
from sklearn.ensemble import RandomForestClassifier as RFC
from sklearn.naive_bayes import GaussianNB as GNB
from sklearn.neighbors import KNeighborsClassifier as KNC
from sklearn.neural_network import MLPClassifier as MLPC
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier as DTC
# other sklearn
from sklearn.base import BaseEstimator
from sklearn.feature_selection import SelectFromModel
from sklearn.metrics import accuracy_score, f1_score
from sklearn.metrics import precision_score, recall_score
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV # StratifiedKFold <----
from sklearn.preprocessing import MinMaxScaler

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
    
    HYPER = {"gradient-boosting": {0: {"max_iter": [100, 150, 200]},
                                   1: {"learning_rate": [0.05, 0.1, 0.2],
                                       "max_iter": [100, 150, 200],
                                       "max_leaf_nodes": [15, 31, 63],
                                       "min_samples_leaf": [10, 20, 30]}},
             # ---
             "random-forest": {0: {}, # no optimization
                               1: {"n_estimators": [100, 150, 200],
                                   "criterion": ["gini", "entropy"],
                                   "min_samples_split": [2, 3, 4],
                                   "min_samples_leaf": [1, 2],
                                   "max_features": ["sqrt", "log2"]}},
             # ---
             "naive-bayes": {0: {}, # no optimization
                             1: {}}, # no optimization
             # ---
             "k-nearest-neighbors": {0: {}, # no optimization
                                     1: {"n_neighbors": [5, 8, 10],
                                         "algorithm": ["ball_tree", "kd_tree"], # stoch
                                         "p": [1, 2]}}, # p -> minkowski power
             # ---
             "neural-network": {0: {"max_iter": [500, 1000, 2000]},
                                1: {"hidden_layer_sizes": [(100,), # 1 hidden layer
                                                           (100, 50), # 2 hidden layers
                                                           (100, 50, 20)], # 3 hidden
                                    "activation": ["relu", "tanh"],
                                    "solver": ["adam", "sgd"],
                                    "learning_rate": ["constant", "invscaling"],
                                    "max_iter": [500, 1000, 2000],
                                    "learning_rate_init": [0.001, 0.005],
                                    "early_stopping": [False, True]}},
             # ---
             "support-vector-machine": {0: {}, # no optimization
                                        1: {"C": [0.1, 1, 10],
                                            "kernel": ["linear", "poly", "rbf"],
                                            "degree": [2, 3, 4],
                                            "decision_function_shape": ["ovo", "ovr"]}},
             # ---
             "decision-tree": {0: {}, # no optimization
                               1: {"criterion": ["gini", "entropy"],
                                   "min_samples_split": [2, 3, 4],
                                   "min_samples_leaf": [1, 2]}}}
    
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
            The proportion of the data to be used for testing
        random_state: int (default=0)
            Controls the split of dataset into training, testing and validation
        init: bool (default=True)
            Whether the model being built belongs to the first set of models (note
            that the second set of models are fed on datasets which encompass
            predictions made by the first set of models)
        final_model: bool (deafult=False)
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
    
    @staticmethod
    def _read_csv(data_path: str) -> pd.DataFrame:
        """
        Reads and returns the dataset pointed to by <data_path>. All columns except
        for the last are casted to np.float32, reducing RAM usage by approximately 50%.
        Returns the newly created dataset.
        
        Parameters
        ----------
        data_path: str
            The path to the .csv file containing the data feeding the model
        """
        # display state of the process on screen (read data)
        print(f"\nReading data from {data_path!r}...")
        # construct mapping for data types
        map_ = collections.defaultdict(lambda: np.float32)
        map_["Function"] = object
        # read .csv (with types specified by <map_>) and return the resulting dataset
        return pd.read_csv(data_path, dtype=map_)
    
    def __repr__(self) -> str:
        """
        Returns the string representation of the object.
        """
        c = self.__class__.__name__
        d, m, a, t = self.data_path, self.models_dir, self.algorithm, self.test_size
        r, i, f = self.random_state, self.init, self.final_model
        res = f"{c}({d!r}, {m!r}, {a!r}, {t!r}, {r!r}, {i!r}, {f!r})"
        return res
        
    def _save_xtest_indices(self, x_test: pd.DataFrame, txt_name: str) -> None:
        """
        Saves the indices of the rows assigned for testing purposes. It is only
        called when <self.final_model> and <self.init> are set to True. Its purpose
        is to test the second set of models (the ones fed on an expanded dataset
        encompassing predictions made by the first set of models) on the same data
        used to test the first set of models.
        
        Parameters
        ----------
        x_test: pd.DataFrame
            The feature vectors of the testing portion of the data
        txt_name: str
            The name to be given to the .txt file storing the indices
        """
        # save indices to a .txt file for future use (case study)
        with open(txt_name + ".txt", "w") as fout:
            for ind in x_test.index:
                fout.write(str(ind) + "\n")
        
    def _train_test_split_init(self) -> tuple:
        """
        Executes the division of the dataset into training and testing. Returns a
        tuple containing 4 numpy arrays: 2 arrays of features (one for training,
        other for testing) and 2 arrays of labels (one for training, other for
        testing). It is used when <self.init> is set to True.
        """
        # display state of the process on screen (train-test split)
        print("Performing 'init' train-test split on the dataset...")
        # split the dataset into features and labels
        df_feats = self._dataset.iloc[:, 1:-1]
        df_labels = self._dataset.iloc[:, -1]
        # split the data into training and testing sets
        x_trn, x_tst, y_trn, y_tst = train_test_split(df_feats,
                                                      df_labels,
                                                      stratify=df_labels,
                                                      test_size=self.test_size,
                                                      random_state=self.random_state)
        # save <x_tst> indices to a .txt file if <self.final_model> is set to True
        if self.final_model:
            txt_name = self.models_dir + "/test-indices-" + self._name_set
            self._save_xtest_indices(x_tst, txt_name)
        # return tuple of numpy arrays
        return x_trn, x_tst, y_trn, y_tst
    
    def _train_test_split_cs(self) -> tuple:
        """
        Executes the division of the dataset into training and testing. Returns a
        tuple containing 4 numpy arrays: 2 arrays of features (one for training,
        other for testing) and 2 arrays of labels (one for training, other for
        testing). It is used when <self.init> is set to False.
        """
        # display state of the process on screen (train-test split)
        print("Performing 'cs' train-test split on the dataset...")
        # read test indices
        txt_name = "../models" + "/test-indices-" + self._name_set
        indices = [int(line.strip()) for line in open(txt_name+".txt").readlines()]
        # split the data into training and testing sets
        train, test = self._dataset.drop(indices), self._dataset.take(indices)
        x_trn, y_trn = train.iloc[:, 1:-1], train.iloc[:, -1]
        x_tst, y_tst = test.iloc[:, 1:-1], test.iloc[:, -1]
        # return tuple of numpy arrays
        return x_trn, x_tst, y_trn, y_tst #np.ravel???
    
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
        # display state of the process on screen (scale data)
        print("Scaling training and testing data...")
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
        # display state of the process on screen (select features)
        print("Selecting features...")
        # fit selector on a RandomForestClassifier and retrieve support vector
        selector = SelectFromModel(estimator=RFC()).fit(x_train, y_train)
        support_vector = selector.get_support()
        # display number of selected features
        print(f"- num_selected = {len(support_vector[support_vector==True])}")
        # save support vector to a .joblib file if <final_model> is set to True
        if self.final_model:
            support_name = self.models_dir + "/support-" + self._name_set
            joblib.dump(support_vector, support_name+".joblib")
        # return shrunken versions of <x_train> and <x_test>
        return x_train[:, support_vector], x_test[:, support_vector]
    
    def _optimize_hyperparameters(self,
                                  x_train: np.ndarray,
                                  y_train: np.ndarray) -> dict:
        """
        Computes the best combination of hyperparameters, for a given estimator and
        training set, by fitting a GridSearchCV object. Returns a dictionary object
        whose keys are the names of the hyperparameters and whose values are the
        best value found for the respective hyperparameter.
        
        Parameters
        ----------
        x_train: np.ndarray
            The feature vectors of training portion of the data
        y_train: np.ndarray
            The label vector of the training portion of the data
        """
        # display state of the process on screen (optimize hyperparameters)
        print("Optimizing hyperparameters...")
        # select hyperparameter grid based on <algorithm> and <final_model>
        param_grid = MLModel.HYPER[self.algorithm][self.final_model]
        # skip optimization if <param_grid> is empty (avoids unnecessary fit)
        if not param_grid:
            best_params = {}
        # otherwise, determine optimal combination of hyperparameters
        else:
            # (cv: "For integer/None inputs, if the estimator is a classifier and y is
            # either binary or multiclass, StratifiedKFold is used. In all other cases,
            # KFold is used. These splitters are instantiated with shuffle=False so the
            # splits will be the same across calls.")
            grid = GridSearchCV(estimator=self._estimator(),
                                param_grid=param_grid,
                                scoring="f1_macro",
                                refit=False).fit(x_train, y_train)
            # assign <grid.best_params_> to <best_params> to match the first case
            best_params = grid.best_params_
        # display and return best combination of hyperparameters
        print(f"- {best_params = }")
        return best_params
    
    def _test_model(self,
                    model: BaseEstimator,
                    x_test: np.ndarray,
                    y_test: np.ndarray) -> None:
        """
        Tests the model <model> on the testing data <x_test> and <y_test>. Computes
        the accuracy score, precision score, recall score and f1-score of the model.
        
        Parameters
        ----------
        model: BaseEstimator
            An estimator fitted on training data
        x_test: np.ndarray
            The feature vectors of testing portion of the data
        y_test: np.ndarray
            The label vector of the testing portion of the data
        """
        # display state of the process on screen (test model)
        print(f"\nTesting {self.algorithm.upper()} model ({self._name_set})\n---")
        # compute predictions
        y_pred = model.predict(x_test)
        # compute metrics and display them on screen
        accuracy = accuracy_score(y_test, y_pred)
        precision = precision_score(y_test, y_pred, average="macro")
        recall = recall_score(y_test, y_pred, average="macro")
        f1_scr = f1_score(y_test, y_pred, average="macro")
        print(f"Metrics on testing data:\n{accuracy = :.2%}\n{precision = :.2%}\n"
              f"{recall = :.2%}\nf1-score = {f1_scr:.2%}\n")
    
    def build_model(self) -> None:
        """
        Builds an ML model given the algorithm and the data provided by the user.
        To do so, splits the dataset in training and testing sets, scales them using
        a MinMaxScaler, selects the most important features by fitting a sklearn's
        SelectFromModel instance on a RandomForestClassifier (fitted on the training
        data) and optimizes the hyperparameters of the estimator chosen by the user.
        The model is then built using these hyperparameters and the training data.
        Finally, the model is tested on the testing data.
        """
        # display state of the process on screen (build model)
        print(f"\nBuilding {self.algorithm.upper()} model ({self._name_set})\n---")
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
        print("Fitting estimator on best combination of hyperparameters...")
        model = self._estimator(**params).fit(x_train, y_train)
        # save model to a .joblib file if <final_model> is set to True
        if self.final_model:
            model_name = self.models_dir + "/model-" + self._name_set
            joblib.dump(model, model_name+".joblib")
        # test the fitted estimator (model) on the testing data
        self._test_model(model, x_test, y_test)
