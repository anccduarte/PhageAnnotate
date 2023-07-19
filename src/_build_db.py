# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
from ml_dataset import MLDataset

"""
IMPORTANT
---
The module must be executed in the command line. Importing individual functions
may throw an error, as some of them include global variables which are only
accessible to this particular module. These variables are: "DF_ALL", "TO_ADD",
"SEQUENCES" and "DATABASE".
"""

def get_protein_name(nfile: str) -> str:
    """
    Returns the protein name contained in <nfile>.
    
    Parameters
    ----------
    nfile: str
        The name of the file containing the protein's name
    """
    _, *temp = nfile.split("_")
    temp[-1] = temp[-1].split(".")[0]
    return "-".join(temp)

def build_df_from_dir(dir_: str, ttable: str, icodons: tuple) -> pd.DataFrame:
    """
    Constructs and returns a pd.DataFrame containing a featurized version
    of the sequences present in the directory "<SEQUENCES>/<dir_>". It also
    updates the global variable "DF_ALL".
    
    Parameters
    ----------
    dir_: str
        The subdirectory containing the sequences
    ttable: str
        The identifier of the translation table to be used
    icodons: tuple
        A tuple of possible initiation codons given the table <ttable>
    """
    # allow to mutate "DF_ALL"
    global DF_ALL
    # initialize df for functional class
    df_class = pd.DataFrame()
    # main loop -> iterate through files in <dir_>
    for file in os.scandir(f"{SEQUENCES}/{dir_}"):
        # guard clause -> ignore files not ending in ".fasta" 
        if not file.name.endswith("fasta"):
            continue
        # get name of the protein coded by the DNA sequences in "file"
        prot_name = get_protein_name(nfile=file.name)
        # featurize DNA sequences in <file>
        df = MLDataset(file=f"{SEQUENCES}/{dir_}/{file.name}",
                       prot_name=prot_name,
                       ttable=ttable,
                       icodons=icodons).build_dataset()
        # concat "df" to "df_class" and "df_all"
        df_class = pd.concat([df_class, df])
        DF_ALL = pd.concat([DF_ALL, df])
        print(f"- Added {prot_name!r} sequences")
    # return "df_class"
    return df_class
        
def build_dbs(dirs: dict, ttable: str, icodons: tuple) -> None:
    """
    Constructs pd.DataFrames containing featurized versions of the sequences
    present in the directories <dirs>, and saves them to .csv files. It also
    updates the global variable "TO_ADD".
    
    Parameters
    ----------
    dirs: dict
        Dictionary containing the names of the directories where the DNA
        sequences to be featurized are stored
    ttable: str
        The identifier of the translation table to be used
    icodons: tuple
        A tuple of possible initiation codons given the table <ttable>
    """
    # allow to mutate "TO_ADD"
    global TO_ADD
    # main loop -> iterate through directories
    for dir_ in dirs:
        func_class = dirs[dir_]
        print(f"Building dataset of {func_class!r} proteins...")
        # build dataframe (featurized sequences in <dir_>)
        df_class = build_df_from_dir(dir_, ttable, icodons)
        # save "df_class" to a .csv file
        df_class.to_csv(f"{DATABASE}/{dir_}.csv")
        # update "to_add" with the name of the functional class
        TO_ADD += [func_class] * df_class.shape[0]
        # separate logs of different functional classes
        print("---")
        

if __name__ == "__main__":
    
    import utils
    from pathlib import Path
    
    # get command line argument
    args = utils.get_args(("-init",))
    init = args.init
    
    # verify validity of <init> (1)
    if init is None:
        raise Exception("<init> has no default value. Please, do:\n"
                        ">>> python _build_db.py -init <init>")
    
    # verify validity of <init> (2)
    if init not in {"yes", "no"}:
        raise ValueError(f"{init!r} is not valid for 'init'."
                         "Choose one of {{'yes', 'no'}}.")
    
    # initialize global variables
    SEQUENCES = "../sequences" if init == "yes" else "../sequences_cs"
    DATABASE = "../database" if init == "yes" else "../database_cs"
    
    # initialize database (if not already initialized)
    Path(DATABASE).mkdir(exist_ok=True)
    
    # directories where the DNA sequences are stored
    directories = {"dna_modification": "dna-modification",
                   "dna_replication": "dna-replication",
                   "lysis": "lysis",
                   "lysogeny_repressor": "lysogeny-repressor",
                   "packaging": "packaging",
                   "structural": "structural"}
    
    # add the dataset "other" to <directories> if not <init>
    if init == "no":
        directories["other"] = "other"
    
    # tuple of initiation codons for translation table 11
    icodons = ("TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG")
    
    # initialize second set of global variables
    DF_ALL = pd.DataFrame()
    TO_ADD = []
    
    # build database
    build_dbs(dirs=directories, ttable="11", icodons=icodons)
    
    # delete column with protein names from "DF_ALL"
    del DF_ALL["Function"]
    
    # add new column to "DF_ALL" (functional classes)
    DF_ALL["Func-Class"] = TO_ADD

    # save "DF_ALL" to a .csv file
    print("Saving 'DF_ALL' to 'all.csv'...")
    DF_ALL.to_csv(f"{DATABASE}/all.csv")
    print("DONE")
    