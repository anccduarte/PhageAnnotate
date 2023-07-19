# -*- coding: utf-8 -*-

import pandas as pd
from Bio import SeqIO, SeqRecord

def write_other(record: SeqRecord) -> None:
    """
    It appends the information contained within <record> to the file
    "../../sequences_cs/txid2731619_other.fasta".
    
    Parameters
    ----------
    record: SeqRecord
        The .fasta record of the putative hypothetical protein whose
        function was successfully predicted
    """
    # construct new description for the record
    keep = record.description.split(" | ")[:2]
    description = " | ".join(keep) + " | " + "other"
    # open appropriate .fasta file and append new sequence to it
    sub_path = f"other/txid2731619_other.fasta"
    with open(f"../../sequences_cs/{sub_path}", "a") as fout:
        fout.write(f">{description}\n{record.seq}\n\n")

def write_function(func_class: str,
                   function: str,
                   record: SeqRecord) -> None:
    """
    Given <func_class> and <function>, it appends the information contained
    within <record> to the appropriate file. The file is pointed to by the
    path "../../sequences_cs/<func_class>/txid2731619_<function>.fasta"
    
    Parameters
    ----------
    func_class: str
        The predicted functional class
    function: str
        The predicted specific function
    record: SeqRecord
        The .fasta record of the putative hypothetical protein whose
        function was not successfully predicted
    """
    # construct new description for the record
    keep = record.description.split(" | ")[:2]
    description = " | ".join(keep) + " | " + function
    # rename <func_class> and <function>
    func_class = func_class.replace("-","_")
    function = function.replace("-","_")
    # open appropriate .fasta file and append new sequence to it
    sub_path = f"{func_class}/txid2731619_{function}.fasta"
    with open(f"../../sequences_cs/{sub_path}", "a") as fout:
        fout.write(f">{record.description}\n{record.seq}\n\n")
        
def distribute_all(path_to_fasta: str, path_to_csv: str) -> tuple:
    """
    Ditributes all .fasta records in <path_to_fasta> to the appropriate
    .fasta files according to the prediction contained in <path_to_csv>.
    
    Parameters
    ----------
    path_to_fasta: str
        The path to file of putative hypothetical protein records
    path_to_csv: str
        The path to the file containing the function predictions of the
        sequences in <path_to_fasta>
    """
    # set counters
    other, other_function = 0, 0
    # read files
    fasta_records = SeqIO.parse(path_to_fasta, "fasta")
    predictions = pd.read_csv(path_to_csv)
    # main loop -> parse files and distribute sequences
    for rec, (_, pred) in zip(fasta_records, predictions.iterrows()):
        # get functional class prediction
        func_class = pred["Functional-Class"]
        # write to the file "../../sequences_cs/txid2731619_other.fasta"
        if func_class == "other":
            write_other(rec)
            other += 1
        else:
            # get function and write new sequence to the appropriate file
            function = pred["Function"]
            if function == "other-function":
                other_function += 1
            write_function(func_class, function, rec)
    # return counts
    return other, other_function
        
             
if __name__ == "__main__":
        
    import sys
    sys.path.append("..")
    import utils
    from pathlib import Path
    
    # initialize directory for storing sequences whose functional class is "other"
    Path("../../sequences_cs/other").mkdir(exist_ok=True)
    
    # get command line arguments
    args = utils.get_args(("-path_to_fasta",),
                          ("-path_to_csv",))
    fasta = args.path_to_fasta
    csv = args.path_to_csv
    
    # check validity of command line arguments
    if fasta is None or csv is None:
        fa, csv = "path_to_fasta", "path_to_csv"
        raise Exception(f"<{fa}> and <{csv}> have no default values. Please do:"
                        f">>> python _distribute_seqs.py -{fa} <{fa}> -{csv} <{csv}>")
    
    # distribute sequences
    o, of = distribute_all(path_to_fasta=fasta, path_to_csv=csv)
    print(f"other - {o}\nother-function - {of}")
