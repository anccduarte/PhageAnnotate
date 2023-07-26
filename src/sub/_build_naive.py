# -*- coding: utf-8 -*-

import collections
import functools as ft
import numpy as np
import pandas as pd
import sys
sys.path.append("..")
from Bio import SeqIO
from collect_sequences import CollectSequences
from ml_dataset import MLDataset
from ml_model import MLModel
from pathlib import Path

# COLLECT & FILTER SEQUENCES
# ---

def collect() -> None:
    """
    Collects all sequences present in the genome records of bacteriophages.
    The resulting .fasta file is then filtered to remove duplicates, and a
    new .fasta file containing non-redundant sequences is generated.
    """
    # inititialize directory for storing the .fasta file
    Path("naive/sequences").mkdir(exist_ok=True)
    # collect sequences
    CollectSequences(db="naive/sequences/all",
                     cname="all",
                     terms=["not a protein common name"],
                     negatives=True).get_sequences(taxid="2731619")

# BUILD 1ST DATASET OF FEATURIZED SEQUENCES
# (all phage sequences - excluding duplicates - are represented in this 1st
# dataset, amounting to 3179455 sequences)
# ---

def _retrieve_protein_names(nfile: str) -> list:
    """
    Given the path to a .fasta file containing DNA sequences, it retrieves
    the names of the proteins associated to all sequences and stores them
    in a list object.

    Parameters
    ----------
    nfile: str
        The name of the .fasta file containing the DNA sequences
    """
    # initialize list object for storing protein names
    protein_names = []
    # loop through the records in <nfile> and add descriptions
    records = SeqIO.parse(nfile, "fasta")
    for rec in records:
        name = " | ".join(rec.description.split(" | ")[2:])
        protein_names.append(name)
    # return list of protein names
    return protein_names
    
def build_db() -> None:
    """
    Builds a database consisting of a dataset of featurized sequences (the
    sequences whose featurization is performed are the ones present in the
    file generated by "collect").
    """
    # specify name of the file
    nfile = "naive/sequences/all/txid2731619_all.fasta"
    # retrieve the protein names associated to the sequences in <nfile>
    protein_names = _retrieve_protein_names(nfile)
    # featurize sequences present in <nfile>
    icodons = ("TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG")
    mld = MLDataset(nfile, "unknown", "11", icodons) # "11" -> ttable
    dframe = mld.build_dataset()
    # remove bad indices from <protein_names>
    bad_inds = mld._bad_inds
    protein_names = np.delete(np.array(protein_names), bad_inds)
    # substitute last column of <dframe> (all "unknown") by <protein_names>
    dframe["Function"] = protein_names
    # initialize directory for storing the dataset
    Path("naive/database").mkdir(exist_ok=True)
    # save dataframe to a .csv file
    dframe.to_csv("naive/database/all.csv")
    
# BUILD 2ND DATASET OF FEATURIZED SEQUENCES
# (only containing entries where "Function" is a protein name whose number of
# occurrences across all bacteriophage genome records >= 100)
# ---

def _read_csv(csv_path: str) -> pd.DataFrame:
    """
    Reads and returns the dataset pointed to by <csv_path>. All columns except
    for the last are casted to np.float32, reducing RAM usage by approximately 50%.
    Returns the newly created dataset.

    Parameters
    ----------
    data_path: str
        The path to the .csv file containing the data feeding the model
    """
    # construct mapping for data types
    map_ = collections.defaultdict(lambda: np.float32)
    map_["Function"] = object
    # read .csv (with types specified by <map_>) and return dataset
    # (index_col=[0] gets rid of "Unnamed" column -> keeping the "Unnamed" column
    # would originate a problem, since the dataset is saved again in "build_db2",
    # that is, the final dataset would contain 2 "Unnamed" columns)
    return pd.read_csv(csv_path, index_col=[0], dtype=map_)

def _get_prot_names_min100(function_col: pd.Series) -> set:
    """
    Returns a set of all protein names whose number of occurrences across all
    entries of <function_col> is greater than or equal to 100.
    
    Parameters
    ----------
    function_col: pd.Series
        A pd.Series of protein common names
    """
    # build dictionary of names and occurrences
    names_all = {}
    for func in function_col:
        if func not in names_all:
            names_all[func] = 1
        else:
            names_all[func] += 1
    # return set of names whose number of occurrences >= 100
    return {name for name, count in names_all.items() if count >= 100}
        
def build_db2() -> None:
    """
    Constructs a second version of the database by filtering "all.csv" to
    remove entries where "Function" is not included in a list of protein
    names whose number of occurrences is greater than or equal to 100.
    """
    # read original .csv and convert all functions to lower case
    print("Reading 'all.csv'...")
    data = _read_csv("naive/database/all.csv")
    data["Function"] = data["Function"].str.lower()
    # get the common names of proteins whose number of occurrences >= 100
    print("Collecting protein names...")
    prot_names = _get_prot_names_min100(function_col=data["Function"])
    # build and save subset of <data> where "Function" in <prot_names>
    print("Subsetting 'all.csv' and saving result...")
    subset = data.loc[data["Function"].isin(prot_names), :]
    subset.to_csv("naive/database/all_tol100.csv")

# BUILD ML MODEL
# (adequate for the construction of an ML model fed on the data stored
# either in "all.csv" (build_db) or "all_tol100.csv" (build_db2))
# ---
    
def build_model(nfile: str) -> None:
    """
    Builds a gradient boosting model (HGBC) fed on the dataset generated by
    "build_db" ("all.csv") or "build_db2" (all_tol100.csv). The model is
    saved for further use (make predictions on new DNA sequences). It then
    tests the model on 20% of the data and displays the attained results.
    
    Parameters
    ----------
    nfile: str
        The name of the .csv file
    """
    # initialize directory for storing .joblib files
    Path("naive/models").mkdir(exist_ok=True)
    # build model and display results (metrics on testing data)
    MLModel(data_path=f"naive/database/{nfile}",
            models_dir="naive/models",
            algorithm="gradient-boosting",
            test_size=0.2,
            init=True,
            final_model=True).build_model()


if __name__ == "__main__":
    
    import utils
    
    args = utils.get_args(("-action",))
    action = args.action
    
    if action is None:
        raise Exception("<action> has no default value. Please do:\n"
                        ">>> python _build_naive.py -action <action>")
        
    options = {"collect": collect,
               "build-db": build_db,
               "build-db-v2": build_db2,
               "build-model": ft.partial(build_model, "all.csv"),
               "build-model-v2": ft.partial(build_model, "all_tol100.csv")}
    
    if action not in options:
        raise ValueError(f"{action!r} is not valid for 'action'. "
                         f"Choose one of {{{', '.join(options)}}}.")
    
    Path("naive").mkdir(exist_ok=True)
        
    func = options[action]
    func()
