# -*- coding: utf-8 -*-

import pandas as pd
from Bio import SeqIO, SeqRecord

def update_dict(dict_: dict, *args) -> dict:
    """
    Updates the dictionary it takes as argument. Its behavior is dependent on
    the number of optional arguments it takes as input. Raises a ValueError if
    the number of oprional arguments is not in {1, 2}.
    
    Parameters
    ----------
    dict_: dict
        The dictionary instance to be updated
    """
    # ---
    if len(args) == 1:
        func_class = args[0]
        if func_class in dict_:
            dict_[func_class] += 1
        else:
            dict_[func_class] = 1
    # ---
    elif len(args) == 2:
        func_class, function = args
        if func_class not in dict_:
            dict_[func_class] = {function: 1}
        else:
            if function in dict_[func_class]:
                dict_[func_class][function] += 1
            else:
                dict_[func_class][function] = 1
    # ---
    else:
        msg = "The number of optional arguments must be either 1 or 2."
        raise ValueError(msg)
    # ---
    return dict_  

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
        fout.write(f">{description}\n{record.seq}\n\n")
        
def distribute_all(path_to_fasta: str, path_to_csv: str) -> tuple:
    """
    Ditributes all .fasta records in <path_to_fasta> to the appropriate .fasta
    files according to the prediction contained in <path_to_csv>.
    ---
    Returns a tuple of 3 "counters":
    - ---
    - other: number of sequences whose functional class prediction was not
      successful (predicted as "other")
    - ---
    - other_function: number of sequences whose functional class was predicted
      successfully, but whose function was predicted as "other-function"
    - ---
    - good_preds: number of sequences whose functional class and specific
      function were successfully predicted
    
    Parameters
    ----------
    path_to_fasta: str
        The path to file of putative hypothetical protein records
    path_to_csv: str
        The path to the file containing the function predictions of the
        sequences in <path_to_fasta>
    """
    # set counters
    # ("other_function" and "good_preds" are dictionaries of counters)
    other, other_function, good_preds = 0, {}
    # read files
    fasta_records = SeqIO.parse(path_to_fasta, "fasta")
    predictions = pd.read_csv(path_to_csv)
    # main loop -> parse files and distribute sequences
    for rec, (_, pred) in zip(fasta_records, predictions.iterrows()):
        # get functional class prediction
        func_class = pred["Functional-Class"]
        # max_prob(func_class) < thresh_func_class
        if func_class == "other":
            other += 1
        else:
            # get function and write new sequence to the appropriate file
            function = pred["Function"]
            # max_prob(function) < thresh_function
            if function == "other-function":
                other_function = update_dict(other_function, func_class)
            else:
                good_preds = update_dict(good_preds, func_class, function)
                write_function(func_class, function, rec)
    # return counters
    return other, other_function, good_preds
        
             
if __name__ == "__main__":
    
    """
    Note: Before running this module, manually copy the contents of the folder
    "../../sequences" to a new folder "../../sequences_cs". Also remember that
    the sequences present in the latter will be filtered by CD-HIT-EST. Only
    then, run the module "_build_train_db_cs.py".
    """
        
    import pprint
    import sys
    sys.path.append("..")
    import utils
    from pathlib import Path
    
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
    res = distribute_all(path_to_fasta=fasta, path_to_csv=csv)
    other, other_function, good_preds = res
    
    # display results
    print(f"predicted as 'other': {other}")
    print("---")
    print("predicted as 'other-function':")
    pprint.pprint(other_function)
    print("---")
    print("good/useful predictions:")
    pprint.pprint(good_preds)
