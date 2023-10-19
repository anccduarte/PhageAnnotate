# -*- coding: utf-8 -*-

import os
import pandas as pd
import sys
sys.path.append("..")
from Bio import SeqIO
from ml_dataset import MLDataset
from tqdm.auto import tqdm

def build_test_db_func_classes() -> None:
    """
    Assembles a correct database of testing sequences (pertaining functional classes).
    The script "_build_test_db_cs.py" does not accomplish this task correctly, since,
    unexpectedly, sequence descriptions do not univocally characterize sequence records.
    """
    # ---
    def get_protein_name(nfile: str) -> str:
        # get and return the name of the protein
        rest, ext = nfile.split(".")
        _, *name = rest.split("_")
        return "-".join(name)
    # ---
    def get_name_new_file(new_path: str, dataset: str, prot_name: str) -> str:
        # the name of the file to be opened for writing
        prot_name = prot_name.replace("-","_")
        return f"{new_path}/{dataset}/txid2731619_{prot_name}.fasta"
    # ---
    def get_sequences(nfile: str, prot_name: str) -> str:
        # parse records in <nfile> and add information to <sequences> 
        sequences = []
        records = SeqIO.parse(nfile, "fasta")
        for rec in records:
            sequences.append((rec.description, rec.seq, prot_name))
        # return list of tuples (descrip, seq, prot_name)
        return sequences
    # ---
    # initialize paths to the incorrectly constructed database of testing sequences, to the
    # testing data (attained upon the train-test split), and a tuple of the names of the
    # functional classes
    path_to_seqs = "../../_sequences_test_cs/func_classes"
    path_to_test_data = "../../models/test-data"
    datasets = ("dna-modification", "dna-replication", "lysis",
                "lysogeny-repressor", "packaging", "structural")
    # create new path (and directory) for storing the correct testing sequences
    new_path = "../../sequences_test_cs/func_classes"
    os.makedirs(new_path)
    # main loop -> iterate through the datasets (functional classes)
    for dataset in datasets:
        # log info to screen
        print(f"Processing {dataset.upper()!r} sequences...")
        # read testing data attained upon the train-test split
        data_test_original = pd.read_csv(f"{path_to_test_data}-{dataset}.csv",
                                         index_col=False)
        data_test_original = data_test_original.iloc[:, 1:-1].astype("float32")
        # create new path for the specific dataset (functional class)
        ds_under = dataset.replace("-","_")
        os.makedirs(f"{new_path}/{ds_under}")
        # intitialize dictionary ought to contain file handles, list representing sequence
        # sequence records of a particular functional class, and list containing dataframes
        # pertaining that functional class (built with MLDataset)
        fwrite_dic = {}
        putative_testing_sequences = []
        putative_testing_dataframe = []
        # inner loop -> iterate through the files of a given functional class
        for file in os.scandir(f"{path_to_seqs}/{ds_under}"):
            if not file.name.endswith(".fasta"):
                continue
            # get the protein name associated to the file
            prot_name = get_protein_name(nfile=file.name)
            # add sequence records to the list of records referred above
            sequences = get_sequences(nfile=f"{path_to_seqs}/{ds_under}/{file.name}",
                                      prot_name=prot_name)
            putative_testing_sequences += sequences
            # add file handle to the dictionary referred above
            fwrite_dic[prot_name] = open(get_name_new_file(new_path,
                                                           ds_under,
                                                           prot_name), "w")
            # construct dataset pertaining to a functional role belonging to the current
            # functional class and append it to the list referred above
            nfile = f"{path_to_seqs}/{ds_under}/{file.name}"
            data = MLDataset(nfile, "unknown", TTABLE, ICODONS).build_dataset()
            putative_testing_dataframe.append(data)
        # concatenate the datasets to create a single dataframe pertaining the current
        # functional class
        df_putative = pd.concat(putative_testing_dataframe).reset_index(drop=True)
        df_putative = df_putative.iloc[:, :-2]
        # juxtapose rows of the original testing data and the newly created dataframe
        # if a row matches, the respective sequence record is added to the appropriate fasta
        # file in the dictionary of file handles
        for rec, (_, vec1) in tqdm(zip(putative_testing_sequences, df_putative.iterrows())):
            for _, vec2 in data_test_original.iterrows():
                if (vec1 == vec2).all():
                    description, sequence, prot_name = rec
                    f = fwrite_dic[prot_name]
                    f.write(f">{description}\n{sequence}\n\n")
                    break
        # close all file handles in <fwrite_dic>
        for prot_name in fwrite_dic:
            f = fwrite_dic[prot_name]
            f.close()
            
# ----------------------------------------------------------------------------------------------------
            
def build_test_db_all() -> None:
    """
    Assembles a correct database of testing sequences (pertaining ALL). The script
    "_build_test_db_cs.py" does not accomplish this task correctly, since, unexpectedly,
    sequence descriptions do not univocally characterize sequence records.
    """
    # ---
    def get_name_new_file(new_path: str, func_class: str) -> str:
        # the name of the file to be opened for writing
        func_class = func_class.replace("-","_")
        return f"{new_path}/txid2731619_{func_class}.fasta"
    # ---
    def get_sequences(path_to_seqs: str, func_class: str) -> str:
        # construct name of the file
        nfile = f"{path_to_seqs}/txid2731619_{func_class.replace('-','_')}.fasta"
        # parse records in <nfile> and add information to <sequences> 
        sequences = []
        records = SeqIO.parse(nfile, "fasta")
        for rec in records:
            sequences.append((rec.description, rec.seq))
        # return list of tuples (descrip, seq, func_class)
        return sequences
    # ---
    path_to_seqs = "../../_sequences_test_cs/all"
    path_to_test_data = "../../models/test-data-all.csv"
    # ---
    data_test_original = pd.read_csv(path_to_test_data, index_col=False)
    # ---
    new_path = "../../sequences_test_cs/all"
    os.makedirs(new_path)
    # ---
    datasets = ("dna-modification", "dna-replication", "lysis",
                "lysogeny-repressor", "packaging", "structural")
    # ---
    dataframes = [data_test_original[data_test_original["Function"] == dataset]
                  for dataset in datasets]
    # ---
    for i, dataset in enumerate(datasets):
        # ---
        print(f"Processing {dataset.upper()!r} sequences...")
        # ---
        original_df = dataframes[i]
        original_df = original_df.iloc[:, 1:-1].astype("float32")
        # ---
        putative_seqs = get_sequences(path_to_seqs=path_to_seqs, func_class=dataset)
        # ---
        nfile = f"{path_to_seqs}/txid2731619_{dataset.replace('-','_')}.fasta"
        putative_df = MLDataset(nfile, "unknown", TTABLE, ICODONS).build_dataset()
        putative_df = putative_df.iloc[:, :-2]
        # ---
        with open(get_name_new_file(new_path, dataset), "w") as fout:
            # ---
            for rec, (_, vec1) in tqdm(zip(putative_seqs, putative_df.iterrows())):
                for _, vec2 in original_df.iterrows():
                    if (vec1 == vec2).all():
                        description, sequence = rec
                        fout.write(f">{description}\n{sequence}\n\n")
                        break
        

if __name__ == "__main__":
    
    """
    Note: The present script (very inefficiently) fixes "sequences_test_cs". This
    database (remember that it constitutes the sequence database produced via
    executing the script "_build_test_db_cs.py", that is, it is the reconstruction
    of the testing sequence database from the .txt file produced upon performing an
    "init" train-test split on the data), quite mysteriously, contains sequences
    whose featurized representation was not used for testing purposes. This rather
    intriguing mystery emerges from a yet to be confirmed (but very plausible) cause:
    contrary to the expected, the descriptions utilized to identify sequences do not
    univocally characterize them. (Remember that sequence descriptions are built as
    f"> {record_id} | {record_source} | {gene_product}".) This obviously brings
    "intruders" to the reconstructed testing sequence database. This script gets rid
    of such intruders.
    ---
    Note 2: The database originated by "_build_test_db_cs.py" becomes obsolete, but it
    is still needed to execute the present script. Before executing it, rename
    "sequences_test_cs" to "_sequences_test_cs". The script uses the latter to get rid
    of intruders and originates a new database precisely called "sequences_test_cs".
    """
    
    import utils
    
    args = utils.get_args(("-data",))
    data = args.data
    
    if data is None:
        raise ValueError("'data' has no default value")
        
    options = {"func_classes": build_test_db_func_classes,
               "all": build_test_db_all}
    
    if data not in options:
        opts = ", ".join(list(options.keys()))
        raise ValueError(f"'data' must be in {{{opts}}}")
    
    TTABLE = "11"
    ICODONS = ("TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG")
    
    func = options[data]
    func()
