# -*- coding: utf-8 -*-

# basic imports
import joblib
import numpy as np
import os
import pandas as pd
# sklearn estimators
from sklearn.ensemble import HistGradientBoostingClassifier as HGBC
from sklearn.ensemble import RandomForestClassifier as RFC
from sklearn.naive_bayes import CategoricalNB as NaiveBayes
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
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.preprocessing import MinMaxScaler

class MLModel:
    
    """
    Constructs a machine learning (ML) model based on an ML algorithm and data
    provided by the user. The validation and test sizes are also adjustable.
    """
    
    ALGOS = {"gradient-boosting": HGBC,
             "random-forest": RFC,
             "naive-bayes": NaiveBayes,
             "k-nearest-neighbors": KNC,
             "neural-network": MLPC,
             "support-vector-machine": SVC,
             "decision-tree": DTC}
    
    HYPER = {"gradient-boosting": {"learning_rate": [0.05, 0.1],
                                   "max_iter": [100, 150, 200]},
             # ---
             "random-forest": {"hyper1": ["val1", "val2", "val3"],
                               "hyper2": ["val1", "val2", "val3"]},
             # ---
             "naive-bayes": {"hyper1": ["val1", "val2", "val3"],
                             "hyper2": ["val1", "val2", "val3"]},
             # ---
             "k-nearest-neighbors": {"hyper1": ["val1", "val2", "val3"],
                                     "hyper2": ["val1", "val2", "val3"]},
             # ---
             "neural-network": {"hyper1": ["val1", "val2", "val3"],
                                "hyper2": ["val1", "val2", "val3"]},
             # ---
             "support-vector-machine": {"hyper1": ["val1", "val2", "val3"],
                                        "hyper2": ["val1", "val2", "val3"]},
             # ---
             "decision-tree": {"hyper1": ["val1", "val2", "val3"],
                               "hyper2": ["val1", "val2", "val3"]},}
    
    def __init__(self,
                 data_path: str,
                 models_dir: str,
                 algorithm: str,
                 test_size: float,
                 random_state: int = 0,
                 final_model: bool = False) -> None:
        """
        Initializes an instance of MLModel.
        
        Parameters
        ----------
        data_path: str
            The path to a .csv file containing the data that will feed the model
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
        final_model: bool (deafult=False)
            Whether the model to be trained is definitive. If set to True, the
            scaler, the support vector (used for feature selection) and the model
            itself are saved in .joblib files for further use (make predictions on
            new data)
            
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
        self.final_model = final_model
        # attributes
        self._dataset = pd.read_csv(data_path)
        self._name_set = data_path.split("/")[-1].split(".")[0].replace("_","-")
        self._estimator = MLModel.ALGOS[algorithm]
        
    def _train_test_split(self) -> tuple:
        """
        Executes the division of the dataset into training and testing. Returns a
        tuple containing 4 numpy arrays: 2 arrays of features (one for training,
        other for testing) and 2 arrays of labels (one for training, other for
        testing).
        """
        # display state of the process on screen
        print("Performing train-test split on the dataset...")
        # split the dataset into features and labels
        df_feats = self._dataset.iloc[:, 1:-1]
        df_labels = self._dataset.iloc[:, -1]
        # split the data into training and testing sets
        x_trn, x_tst, y_trn, y_tst = train_test_split(df_feats,
                                                      df_labels,
                                                      stratify=df_labels,
                                                      test_size=self.test_size,
                                                      random_state=self.random_state)
        # returns tuple of numpy arrays
        return x_trn, x_tst, np.ravel(y_trn), np.ravel(y_tst)
    
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
        # display state of the process on screen
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
        # display state of the process on screen
        print("Selecting features...")
        # fit selector on a RandomForestClassifier and retrieve support vector
        selector = SelectFromModel(estimator=RFC()).fit(x_train, y_train)
        support_vector = selector.get_support()
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
        # display state of the process on screen
        print("Optimizing hyperparameters...")
        # select hyperparameter grid based on the selected algorithm
        param_grid = MLModel.HYPER[self.algorithm]
        # initialize estimator
        estimator = self._estimator()
        # compute optimal combination of hyperparameters
        grid = GridSearchCV(estimator=estimator,
                            param_grid=param_grid,
                            # by default, n_splits=5
                            cv=StratifiedKFold(shuffle=True,
                                               random_state=self.random_state))
        grid.fit(x_train, y_train)
        # return best combination of hyperparameters
        return grid.best_params_
    
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
        # display state of the process on screen
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
        # display state of the process on screen
        print(f"\nBuilding {self.algorithm.upper()} model ({self._name_set})\n---")
        # train-test split
        x_train, x_test, y_train, y_test = self._train_test_split()
        # scale <x_train> and <x_test>
        x_train, x_test = self._scale_data(x_train, x_test)
        # select most important features
        x_train, x_test = self._select_features(x_train, y_train, x_test)
        # optimize hyperparameters and fit the model based on them
        params = self._optimize_hyperparameters(x_train, y_train)
        model = self._estimator(**params).fit(x_train, y_train)
        # save model to a .joblib file if <final_model> is set to True
        if self.final_model:
            model_name = self.models_dir + "/model-" + self._name_set
            joblib.dump(model, model_name+".joblib")
        # test the fitted estimator on the testing data
        self._test_model(model, x_test, y_test)
