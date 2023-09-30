# -*- coding: utf-8 -*-

import collections
import numpy as np
import os
import pandas as pd
import sys
sys.path.append("..")
from ml_model import MLModel

def get_dataset() -> None:
    """
    Constructs dataset by concatenating all datasets present in
    the directory "../../database".
    """
    # ---
    def _read_csv(path: str) -> pd.DataFrame:
        map_ = collections.defaultdict(lambda: np.float32)
        map_["Description"] = object
        map_["Function"] = object
        return pd.read_csv(path, index_col=[0], dtype=map_)
    # ---
    dfs = []
    common_path = "../../database"
    for file in os.scandir(common_path):
        if not file.name.endswith(".csv"):
            continue
        if "all" in file.name:
            continue
        dfs.append(_read_csv(f"{common_path}/{file.name}"))
    data = pd.concat(dfs)
    data.to_csv("linear/database/all.csv")

def build_model():
    """
    Fits "gradient-boosting" model. The model is fed on the data
    present in "linear/database/all.csv".
    """
    MLModel(data_path="linear/database/all.csv",
            models_dir="linear/models",
            algorithm="gradient-boosting",
            test_size=0.2,
            init=True,
            final_model=True).build_model()
    
if __name__ == "__main__":
    
    from pathlib import Path
    
    # initialize directories
    for dir_ in ("linear", "linear/database", "linear/models"): 
        Path(dir_).mkdir(exist_ok=True)
        
    # construct dataset
    print("\nConstructing dataset...")
    get_dataset()
    print("---\nDONE")
    
    # build "gradient-boosting" model
    build_model()
