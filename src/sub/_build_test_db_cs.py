# -*- coding: utf-8 -*-

import os
from Bio import SeqIO

# --------------------------------------------------------------------------------

def readlines_file(file: str) -> list:
    """
    Reads and returns the lines present in the file pointed to by <file>.
    
    Parameters
    ----------
    file: str
        The path to the file to be read
    """
    with open(file, "r") as fin:
        lines = fin.readlines()
    return [line[:-1] for line in lines]

def build_labels_dict(lines: list) -> dict:
    """
    Receives a list of lines each containing a label and a corresponding
    description. Builds and returns a dictionary instance whose keys are unique
    labels and whose values are lists of descriptions associated to the labels.
    
    Parameters
    ----------
    lines: list
        The list of lines containing labels and descriptions
    """
    # ---
    def get_label(line: str) -> str:
        label, *_ = line.split(":")
        return label[2:].replace("-", "_")
    # ---
    def get_description(line: str) -> str:
        _, *description = line.split(":")
        return ":".join(description)[1:]
    # ---
    def build_dict() -> dict:
        dict_ = {}
        for line in lines:
            label = get_label(line)
            description = get_description(line)
            if label not in dict_:
                dict_[label] = {description}
            else:
                dict_[label].add(description)
        return dict_
    # ---
    return build_dict() 

# --------------------------------------------------------------------------------

def build_seq_db_func_class(file: str, name_dir: str) -> None:
    """
    Builds a set of .fasta files by inspecting the contents of the .txt pointed
    to by <file>. The resulting .fasta files contain the testing sequences
    attained during the train-test split performed upon the construction of the
    initial (init) set of ML models. The .fasta files are stored in <name_dir>.
    Note: to be used for the test data info generated during the train-test
    split of the datasets distinguishing specific functions.
    
    Parameters
    ----------
    file: str
        The name of the .txt file generated during the train-test split
    name_dir: str
        The name of the directory where the .fasta files are to be saved
    """
    # ---
    def get_name_func_class() -> str:
        file_name = file.split("/")[-1]
        func_class = file_name[10:-4] # test-data-[10] AND .txt[-4]
        return func_class.replace("-", "_")
    # ---
    def build_fasta_files() -> None:
        # ---
        func_class = get_name_func_class()
        path_test_seqs = f"{name_dir}/{func_class}"
        os.mkdir(path_test_seqs)
        # ---
        lines_file = readlines_file(file)
        labels_dict = build_labels_dict(lines_file)
        # ---
        for label in labels_dict:
            descriptions = labels_dict[label]
            seqs_file = f"../../sequences/{func_class}/txid2731619_{label}.fasta"
            with open(f"{path_test_seqs}/txid2731619_{label}.fasta", "w") as f:
                records = SeqIO.parse(seqs_file, "fasta")
                for rec in records:
                    if rec.description[1:] not in descriptions:
                        continue
                    f.write(f">{rec.description}\n{rec.seq}\n\n")
    # ---
    build_fasta_files()
    
# --------------------------------------------------------------------------------
    
def build_seq_db_all(file: str, name_dir: str) -> None:
    """
    Builds a set of .fasta files by inspecting the contents of the .txt pointed
    to by <file>. The resulting .fasta files contain the testing sequences
    attained during the train-test split performed upon the construction of the
    initial (init) set of ML models. The .fasta files are stored in <name_dir>.
    Note: to be used for the test data info generated during the train-test
    split of the datasets distinguishing functional classes (ALL).
    
    Parameters
    ----------
    file: str
        The name of the .txt file generated during the train-test split
    name_dir: str
        The name of the directory where the .fasta files are to be saved
    """
    # ---
    lines_file = readlines_file(file)
    labels_dict = build_labels_dict(lines_file)
    # ---
    for label in labels_dict:
        descriptions = labels_dict[label]
        label_folder = f"../../sequences/{label}"
        with open(f"{name_dir}/txid2731619_{label}.fasta", "w") as f:
            for file in os.scandir(label_folder):
                if not file.name.endswith(".fasta"):
                    continue
                records = SeqIO.parse(f"{label_folder}/{file.name}", "fasta")
                for rec in records:
                    if rec.description[1:] not in descriptions:
                        continue
                    f.write(f">{rec.description}\n{rec.seq}\n\n")
    
            
if __name__ == "__main__":
    
    """
    Note 1: This module may be executed after the construction of the "final initial"
    models ("initial" as they refer to the "init" models; "final" since they refer
    to the models built using the optimal algorithm for each dataset). It builds a
    directory "../../sequences_test_cs", consisting of two subdirectories "all" and
    "func_classes", which contain .fasta files corresponding to the test sequences
    determined during the "final initial" train-test split for each dataset.
    ---
    Note 2: Two separate databases ("func_classes" and "all") must be created since
    the testing sequences selected for each functional class are not the same as the
    testing sequences selected for the dataset "ALL", upon the train-test split.
    """
    
    os.mkdir("../../sequences_test_cs")
    
    func_class_dir = "../../sequences_test_cs/func_classes"
    os.mkdir(func_class_dir)
    
    all_dir = "../../sequences_test_cs/all"
    os.mkdir(all_dir)
    
    files = ("test-data-dna-modification.txt",
             "test-data-dna-replication.txt",
             "test-data-lysis.txt",
             "test-data-lysosgeny-repressor.txt",
             "test-data-packaging.txt",
             "test-data-structural.txt",
             "test-data-all.txt")
    # ---
    files_with_dir = (f"../../models/{file}" for file in files)
    
    for file in files_with_dir:
        if "all" in file:
            build_seq_db_all(file, all_dir)
        else:
            build_seq_db_func_class(file, func_class_dir)
