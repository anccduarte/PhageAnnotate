# -*- coding: utf-8 -*-

import os
import pandas as pd
from ml_dataset import MLDataset
from pathlib import Path

# initialize database (if not already initialized)
DATABASE = "../database"
Path(DATABASE).mkdir(exist_ok=True)

# initialize global variables
DF_ALL = pd.DataFrame()
TO_ADD = []

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

def build_df_from_dir(dir_: str) -> pd.DataFrame:
    """
    Constructs and returns a pd.DataFrame containing a featurized version
    of the sequences present in the directory "../sequences/<dir_>".
    
    Parameters
    ----------
    dir_: str
        The subdirectory containing the sequences
    """
    # allow to mutate "DF_ALL"
    global DF_ALL
    # initialize df for functional class
    df_class = pd.DataFrame()
    # main loop -> iterate through files in <dir_>
    for file in os.scandir(f"../sequences/{dir_}"):
        # guard clause -> ignore files not ending in ".fasta" 
        if not file.name.endswith("fasta"):
            continue
        # get name of the protein coded by the DNA sequences in "file"
        prot_name = get_protein_name(nfile=file.name)
        # featurize DNA sequences in <file>
        df = MLDataset(file=f"../sequences/{dir_}/{file.name}",
                       protein_name=prot_name).build_dataset()
        # concat "df" to "df_class" and "df_all"
        df_class = pd.concat([df_class, df])
        DF_ALL = pd.concat([DF_ALL, df])
        print(f"- Added {prot_name!r} sequences")
    # return "df_class"
    return df_class
        
def build_dbs(dirs: dict) -> None:
    """
    Constructs pd.DataFrames containing featurized versions of the
    sequences present in the directories <dirs>, and saves them to .csv
    files.
    
    Parameters
    ----------
    dirs: dict
        Dictionary containing the names of the directories where the DNA
        sequences to be featurized are stored
    """
    # allow to mutate "TO_ADD"
    global TO_ADD
    # main loop -> iterate through directories
    for dir_ in dirs:
        func_class = dirs[dir_]
        print(f"Building dataset of {func_class!r} proteins...")
        # build dataframe (featurized sequences in <dir_>)
        df_class = build_df_from_dir(dir_=dir_)
        # save "df_class" to a .csv file
        df_class.to_csv(f"{DATABASE}/{dir_}.csv")
        # update "to_add" with the name of the functional class
        TO_ADD += [func_class] * df_class.shape[0]
        # separate logs of different functional classes
        print("---")
        

if __name__ == "__main__":
    
    # directories where the DNA sequences are stored
    directories = {"dna_modification": "dna-modification",
                   "dna_replication": "dna-replication",
                   "lysis": "lysis",
                   "lysogeny_repressor": "lysogeny-repressor",
                   "packaging": "packaging",
                   "structural": "structural",
                   "other": "other"}
    
    # build databases
    build_dbs(dirs=directories)
    
    # delete column with protein names from "DF_ALL"
    del DF_ALL["Function"]
    
    # add new column to "DF_ALL" (functional classes)
    DF_ALL["Func-Class"] = TO_ADD

    # save "DF_ALL" to a .csv file
    print("Saving 'DF_ALL' to 'all.csv'...")
    DF_ALL.to_csv(f"{DATABASE}/all.csv")
    print("DONE")
    