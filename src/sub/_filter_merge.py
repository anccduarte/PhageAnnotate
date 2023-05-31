# -*- coding: utf-8 -*-

import re
from Bio import SeqIO

def get_num_tms_per_seq(protein: str, num_files: int) -> list:
    """
    Returns a list object storing the number of transmembrane domains of each
    protein sequence present in the files whose name belongs to the set:
    "results_tmhmm/TMHMM_<protein>s_<n>.html", n in [1, num_files].
    
    Parameters
    ----------
    protein: str
        The common name of the target protein
    num_files: int
        The number of .html files containing the TMHMM results
    """
    # read contents of files
    common = f"results_tmhmm/TMHMM_{protein}s"
    temp = [open(f"{common}_{n}.html").read() for n in range(1, num_files+1)]
    # conacatenate contents of files
    text = "\n".join(temp)
    # get matches of the expression "Number of predicted TMHs:" in "text"
    pattern = re.compile(r"Number of predicted TMHs:")
    matches = re.finditer(pattern, text)
    # compute and return the number of TMs for each protein sequence
    tms = [int(text[match.end()+2]) for match in matches]
    return tms

def filter_merge(protein: str, num_files: int, original_path: str) -> tuple:
    """
    Filters the sequences present in <original_path> based on the results
    attained by running TMHMM. Only the DNA sequences whose protein product
    have at least one transmembrane domain are kept in a new .fasta file.
    Returns a tuple containing the number of sequences present in the
    original file and in the new file.
    
    Parameters
    ----------
    protein: str
        The common name of the target protein
    num_files: int
        The number of .html files containing the TMHMM results
    original_path: str
        The path to the original .fasta file
    """
    # get records in the original file and matches in "tmhmm_results"
    records = SeqIO.parse(original_path, "fasta")
    num_tms_per_seq = get_num_tms_per_seq(protein, num_files)
    # open new file for writing
    file_name = original_path.split("/")[-1]
    with open(file_name, "w") as f:
        # main loop -> create new file (filter sequences)
        num_seqs, total_tms = 0, 0
        for num_tms, record in zip(num_tms_per_seq, records):
            num_seqs += 1
            if num_tms == 0:
                continue
            total_tms += 1
            f.write(f">{record.description}\n{record.seq}\n\n")
    # return number of sequences in both files
    return num_seqs, total_tms
            
    
if __name__ == "__main__":
    
    import sys
    sys.path.append("..")
    import utils
    
    args = utils.get_args(("-protein",),
                          ("-num_files",),
                          ("-original_path",))
    
    protein = args.protein
    num_files = args.num_files
    original_path = args.original_path
    
    if any(arg is None for arg in [protein, num_files, original_path]):
        p, n, o = "protein", "num_files", "original_path"
        e = (f"<{p}>, <{n}> and <{o}> have no default values. Please do:\n"
             f">>> python _split_holins.py -{p} <{p}> -{n} <{n}> -{o} <{o}>")
        raise ValueError(e)
        
    num_original, num_filtered = filter_merge(protein=protein,
                                              num_files=int(num_files),
                                              original_path=original_path)
    
    rem = 1 - (num_filtered / num_original)
    print(f"Number of sequences in the original file: {num_original}")
    print(f"Number of sequences in the new file: {num_filtered}")
    print(f"Percentage of sequences removed: {rem:.2%}")
    