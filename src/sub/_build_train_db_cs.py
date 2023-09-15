# -*- coding: utf-8 -*-

import os
import shutil
import system
sys.path.append("..")
import utils
from Bio import SeqIO

# --------------------------------------------------------------------------------

def build_temp_train_cs() -> None:
    """
    Builds a temporary database of sequences "sequences_train_cs". This database
    consists of two separate databases:
    ---
    1. a database identical to "sequences_cs" called "func_classes";
    2. a database, called "all", attained by concatenating the .fasta files present
    in each folder of "sequences_cs"
    ---
    The construction of this provisory database is vital to replicate the structure
    of the previously built "sequences_test_cs" (makes the filtering process of
    "cd-hit-est-2d" much more straightforward).
    """
    # ---
    def create_seqs_cs_copy() -> None:
        destination = f"{BASE}/sequences_train_cs/func_classes"
        shutil.copytree(f"{BASE}/sequences_cs", destination)
    # ---
    def concatenate_seqs_cs_by_func_class() -> None:
        # ---
        def concatenate_sequences(path_to_all: str, entry: str) -> None:
            func_class_folder = f"{BASE}/sequences_cs/{entry}"
            new_fasta_file = f"{path_to_all}/txid2731619_{entry}.fasta"
            # ---
            with open(new_fasta_file, "w") as f:
                for file in os.scandir(func_class_folder):
                    if not file.name.endswith(".fasta"):
                        continue
                    file_name = f"{func_class_folder}/{file.name}"
                    records = SeqIO.parse(file_name, "fasta")
                    for rec in records:
                        f.write(f">{rec.description}\n{rec.seq}\n\n")
        # ---
        path_to_all = f"{BASE}/sequences_train_cs/all"
        os.mkdir(path_to_all)
        for entry in os.scandir(f"{BASE}/sequences_cs"):
            if not os.path.isdir(entry):
                continue
            concatenate_sequences(path_to_all, entry.name)
    # ---
    def main() -> None:
        create_seqs_cs_copy()
        concatenate_seqs_cs_by_func_class()
    # ---
    main()
    
# --------------------------------------------------------------------------------

def filter_seqs_train_cs() -> None:
    """
    Filters the sequences present in the provisory database "sequences_train_cs" by
    running the software "cd-hit-est-2d" against homologous sequences in the database
    "sequences_test_cs". The filtering process is performed twice:
    ---
    1. the local function "filter_func_classes" acts on the files present in the
    subdirectories "func_classes";
    2. "filter_all" deals with the files present in the subdirectories "all"
    ---
    The attained sequence database is solely consisted of sequences that, upon the
    featurization process, will be used for training (see loading of the training data
    in "src/ml_model.py" -> "train_test_split_cs").
    """
    # ---
    def remove_txts(folder: str) -> None:
        for file in os.scandir(folder):
            if not file.name.endswith(".txt"):
                continue
            os.remove(f"{folder}/{file.name}")
    # ---
    def filter_func_classes() -> None:
        # modify cwd to "../../../cd-hit"
        with utils.change_dir(CD_HIT):
            # initialize variables used throughout the function
            cd_hit_save = "results_train_cs/func_classes"
            os.mkdir(cd_hit_save)
            sequences_train = f"../PhageAnnotate/sequences_train_cs/func_classes"
            sequences_test = f"../PhageAnnotate/sequences_test_cs/func_classes"
            # iterate through the directories in "sequences_train"
            for func_class in os.scandir(sequences_train):
                if not os.path.isdir(func_class):
                    continue
                # if <func_class> is a directory, create a new subdirectory
                os.mkdir(f"{cd_hit_save}/{func_class}")
                # iterate through the .fasta files in sequences_train/<func_class>
                for file in os.scandir(f"{sequences_train}/{func_class}"):
                    if not file.name.endswith(".fasta"):
                        continue
                    # run cd-hit-est-2d
                    train_fasta = f"{sequences_train}/{func_class}/{file.name}"
                    test_fasta = f"{sequences_test}/{func_class}/{file.name}"
                    out_fasta = f"{cd_hit_save}/{func_class}/{file.name}"
                    cmd = f"cd-hit-est-2d -i {test_fasta} -i2 {train_fasta} -o {out_fasta}"
                    os.system(cmd)
                # remove .txt files that are also output of cd-hit-est-2d
                remove_txts(f"{CD_HIT}/{cd_hit_save}/{func_class}")
    # ---
    def filter_all() -> None:
        # modify cwd to "../../../cd-hit"
        with utils.change_dir(CD_HIT):
            # initialize variables used throughout the function
            cd_hit_save = "results_train_cs/all"
            os.mkdir(cd_hit_save)
            sequences_train = f"../PhageAnnotate/sequences_train_cs/all"
            sequences_test = f"../PhageAnnotate/sequences_test_cs/all"
            # iterate through the files in "sequences_train"
            for file in os.scandir(sequences_train):
                if not file.name.endswith(".fasta"):
                    continue
                # run cd-hit-est-2d
                train_fasta = f"{sequences_train}/{file.name}"
                test_fasta = f"{sequences_test}/{file.name}"
                out_fasta = f"{cd_hit_save}/{file.name}"
                cmd = f"cd-hit-est-2d -i {test_fasta} -i2 {train_fasta} -o {out_fasta}"
                os.system(cmd)
            # remove .txt files that are also output of cd-hit-est-2d
            remove_txts(f"{CD_HIT}/{cd_hit_save}")
    # ---
    def main() -> None:
        os.mkdir(f"{CD_HIT}/results_train_cs")
        filter_func_classes()
        filter_all()
    # ---
    main()
        
        
if __name__ == "__main__":
    
    """
    Note: Upon execution of this module, a folder "results_train_cs" will be created at
    "../../../cd-hit". The folder "../../sequences_train_cs" must be renamed to
    "../../_sequences_train_cs" and "results_train_cs" must be copied to "../.." with
    the name "sequences_train_cs".
    """
    
    # used in "build_temp_train_cs"
    BASE = "../.."
    # used in "filter_seqs_train_cs"
    CD_HIT = "../../../cd-hit"
    
    # build temporary database of training sequences (this database contains sequences
    # similar to the ones present in the testing data)
    build_temp_train_cs()
    
    # after building the temporary database, filter the sequences to remove the ones
    # which show similarity to the testing sequences
    filter_seqs_train_cs()
