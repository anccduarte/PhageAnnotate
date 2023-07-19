# -*- coding: utf-8 -*-

import random
from Bio import SeqIO

def get_num_records(path_to_file: str) -> int:
    """
    Returns the number of records contained in <path_to_file>.
    
    Parameters
    ----------
    path_to_file: str
        The path to the file containing the records.
    """
    records = SeqIO.parse(path_to_file, "fasta")
    counter = 0
    for _ in records:
        counter += 1
    return counter

def get_random_ints(num: int, lower_lim: int, upper_lim: int) -> list:
    """
    Returns a list of size min(num, [upper_lim - lower_lim]) containing
    unique randomly generated integers belonging to the interval
    [lower_lim, upper_lim).
    
    Parameters
    ----------
    num: int
        The number of elements to be sampled
    lower_lim: int
        The "lower_lim" in [lower_lim, upper_lim)
    upper_lim: int
        The "upper_lim" in [lower_lim, upper_lim)
    """
    interval = range(lower_lim, upper_lim)
    if num < len(interval):
        chosen = random.sample(interval, k=num)
    else:
        chosen = interval
    return chosen
    
def sample_sequences(func_class: str, function: str, num_seqs: int) -> None:
    """
    Samples a maximum of <num_seqs> sequences from the pool of sequences
    present in the file "<func_class>/txid2731619_<function>.fasta". The
    sampled sequences are stored in a new .fasta file.
    
    Parameters
    ----------
    num_seqs: int
        The number of sequences to be sampled
    """
    # rename <func_class> and <function>
    func_class = func_class.replace(" ","_").replace("-","_")
    function = function.replace(" ","_").replace("-","_")
    # construct path to file
    base = "../../sequences_cs"
    path_to_file = f"{base}/{func_class}/txid2731619_{function}.fasta"
    # get indices of the sequences to be chosen
    num_records = get_num_records(path_to_file)
    indices = set(get_random_ints(num_seqs, 0, num_records))
    # get records and sample
    records = SeqIO.parse(path_to_file, "fasta")
    with open(f"{base}/{func_class}/{function}_sample.fasta", "w") as fout:
        for i, record in enumerate(records):
            if i not in indices:
                continue
            fout.write(f">{record.description}\n{record.seq}\n\n")
    
    
if __name__ == "__main__":
    
    import sys
    sys.path.append("..")
    import utils
    
    random.seed(0)
    
    args = utils.get_args(("-func_class",),
                          ("-function",),
                          ("-num_seqs",))
    func_class = args.func_class
    function = args.function
    num_seqs = args.num_seqs
    
    if func_class is None or function is None or num_seqs is None:
        fc, f, ns = "func_class", "function", "num_seqs"
        raise Exception(f"<{fc}>, <{f}> and <{ns}> have no default values. Please, do:\n"
                        f">>> python _sample_seqs.py -{fc} <{fc}> -{f} <{f}> -{ns} <{ns}>")
        
    sample_sequences(func_class, function, int(num_seqs))
    