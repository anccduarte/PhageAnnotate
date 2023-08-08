# -*- coding: utf-8 -*-

import random
random.seed(0)
from Bio import SeqIO

def get_num_records(path_to_file: str) -> int:
    """
    Returns the number of records contained in <path_to_file>.
    
    Parameters
    ----------
    path_to_file: str
        The path to the file containing the records
    """
    records = SeqIO.parse(path_to_file, "fasta")
    count = 0
    for _ in records:
        count += 1
    return count

def get_indices(num: int, lower_lim: int, upper_lim: int) -> set:
    """
    Returns a set of size min(num, [upper_lim - lower_lim]) containing
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
    return set(chosen)

def write_to_file(original_path: str, new_path: str, indices: set) -> int:
    """
    Writes sampled sequences from <original_path> to a new .fasta file
    pointed to by <new_path>. Returns the number of DNA sequences present
    in the original file.
    
    Parameters
    ----------
    original_path: str
        The path to the original file
    new_file: str
        The path to the new file
    indices: set
        The indices of the sequences to be sampled from <original_path>
    """
    # get sequence records
    records = SeqIO.parse(original_path, "fasta")
    # write sampled sequences to the new file
    with open(new_path, "w") as fout:
        for i, rec in enumerate(records):
            if i not in indices:
                continue
            fout.write(f">{rec.description}\n{rec.seq}\n\n")
    # return number of sequences in the original file
    return i+1
    
def sample_sequences(base_dir: str,
                     func_class: str,
                     function: str,
                     num_seqs: int) -> int:
    """
    Samples a maximum of <num_seqs> sequences from the pool of sequences
    present in the file "<func_class>/txid2731619_<function>.fasta". The
    sampled sequences are stored in a new .fasta file. Returns the number
    of sequences present in the original file.
    
    Parameters
    ----------
    base_dir: str
        The base directory where the sequences cam be found
    func_class: str
        The functional class associated to the target sequences
    function: str
        The specific function associated to the target sequences
    num_seqs: int
        The number of sequences to be sampled
    """
    # rename <func_class> and <function>
    func_class = func_class.replace(" ","_").replace("-","_")
    function = function.replace(" ","_").replace("-","_")
    # construct path to the original and new files
    original_path = f"{base_dir}/{func_class}/txid2731619_{function}.fasta"
    new_path = f"{base_dir}/{func_class}/{function}_sample.fasta"
    # get indices of the sequences to be sampled from <original_path>
    num_records = get_num_records(original_path)
    indices = get_indices(num_seqs, 0, num_records)
    # write sampled sequences to a new file
    num_original = write_to_file(original_path, new_path, indices)
    # return number of sequences in the original file
    return num_original
    
    
if __name__ == "__main__":
    
    import sys
    sys.path.append("..")
    import utils
    
    args = utils.get_args(("-base_dir",),
                          ("-func_class",),
                          ("-function",),
                          ("-num_seqs",))
    base_dir = args.base_dir
    func_class = args.func_class
    function = args.function
    num_seqs = int(args.num_seqs)
    
    if any(arg is None for arg in (base_dir, func_class, function, num_seqs)):
        b, fc, f, ns = "base_dir", "func_class", "function", "num_seqs"
        raise Exception(f"<{b}>, <{fc}>, <{f}> and <{ns}> have no default values. Please, do:\n"
                        f">>> python _sample_seqs.py -{b} <{b}> -{fc} <{fc}> -{f} <{f}> -{ns} <{ns}>")
        
    num_original = sample_sequences(base_dir, func_class, function, num_seqs)
    print(f"{num_original = }\n{num_seqs = }")
    