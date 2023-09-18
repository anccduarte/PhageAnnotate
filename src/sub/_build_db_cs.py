# -*- coding: utf-8 -*-

import os
import pandas as pd
import sys
sys.path.append("..")
from ml_dataset import MLDataset

def build_db(folder: str) -> None:
    """
    Builds the database of training sequences that will eventually feed the
    case-study (cs) models.
    
    Parameters
    ----------
    folder: str
        The (relative) path to the .fasta files containing the sequences to
        be featurized
    """
    # ---
    def get_prot_name(file_name: str) -> str:
        _, *temp = nfile.split("_")
        temp[-1] = temp[-1].split(".")[0]
        return "-".join(temp)
    # ---
    def get_class_name() -> str:
        *_, class_name = folder.split("/")
        return class_name
    # ---
    # initialize "ttable", "icodons" and "df_class"
    ttable = "11"
    icodons = ("TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG")
    df_class = pd.DataFrame()
    # iterate through the sequences in "folder"
    for file in os.scandir(folder):
        if not file.name.endswith(".fasta"):
            continue
        # featurize DNA sequences present in folder/<file>
        df = MLDataset(file=f"{folder}/{file.name}",
                       prot_name=get_prot_name(file.name),
                       ttable=ttable,
                       icodons=icodons).build_dataset()
        # concat "df" to "df_class"
        df_class = pd.concat([df_class, df])
    # save "df_class" to a .csv file
    df_class.to_csv(f"{DATABASE}/{get_class_name()}.csv")
        
        
if __name__ == "__main__":
    
    """
    Note: This constitutes the final step before constructing the case-study
    (cs) models. Remember that the sequences used to construct the datasets
    are the ones present in "../../sequences_train_cs" attained after running
    the software "cd-hit-est-2d" on the provisory training sequence database
    to remove sequences similar to the ones present in the homologous files in
    the database "../../sequences_test_cs".
    """
    
    DATABASE = "../../database_cs"
    os.mkdir(DATABASE)
    
    common = "../../sequences_train_cs"
    # ---
    func_classes = ("dna_modification", "dna_replication", "lysis",
                    "lysogeny_repressor", "packaging", "structural")
    # ---
    folders = [f"{common}/func_classes/{fc}" for fc in func_classes]
    folders.append(f"{common}/all")
    
    for folder in folders:
        build_db(folder=folder)
