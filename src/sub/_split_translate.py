# -*- coding: utf-8 -*-

import io
from Bio import Seq, SeqIO

def get_handle(protein: str, suffix: int) -> io.TextIOWrapper:
    """
    Returns a file handle to '<protein>s/<protein>s_<suffix>.fasta'.
    
    Parameters
    ----------
    protein: str
        The common name of the protein
    suffix: int
        The suffix in 'holins/holins_<suffix>.fasta'
    """
    return open(f"{protein}s/{protein}s_{suffix}.fasta", "w")

def translate(seq: Seq.Seq, ttable: str, icodons: tuple) -> Seq.Seq:
    """
    Translates the DNA sequence <seq>, and returns the resulting protein
    sequence.
    
    Parameters
    ----------
    seq: Seq.Seq
        The DNA sequence
    ttable: str
        The identifier of the translation table to be used
    icodons: tuple
        A tuple of possible initiation codons given the table <ttable>
    """
    valid = lambda s: len(s) % 3 == 0
    if len(seq) % 3 != 0:
        icodon = seq[:3]
        if icodon in icodons:
            while not valid(seq): seq = f"{seq}N"
        else:
            while not valid(seq): seq = f"N{seq}"
    tseq = Seq.Seq(seq).translate(table=ttable)[:-1]
    return tseq

def split_translate(cutoff: int,
                    func_class: str,
                    protein: str,
                    ttable: str,
                    icodons: tuple) -> int:
    """
    Breaks '../../sequences/<func_class>/txid2731619_<protein>.fasta'
    into a number of files equal to ceil(num_sequences/<cutoff>). Also,
    the sequences in the new files correspond to the translated versions
    of the sequences present in the original .fasta file. Returns the
    number of files created.
    
    Parameters
    ----------
    cutoff: int
        The maximum number of sequences to be contained in each file
    func_class: str
        The functional class of the target protein
    protein: str
        The common name of the target protein
    ttable: str
        The identifier of the translation table to be used
    icodons: tuple
        A tuple of possible initiation codons given the table <ttable>
    """
    suffix = 1
    curr = 0
    f = get_handle(protein, suffix)
    path = f"../../sequences/{func_class}/txid2731619_{protein}.fasta"
    for record in SeqIO.parse(f"{path}", "fasta"):
        if curr == cutoff:
            curr = 0
            suffix += 1
            f.close()
            f = get_handle(protein, suffix)
        else:
            amino_seq = translate(record.seq, ttable, icodons)
            f.write(f">{record.description}\n{amino_seq}\n\n")
            curr += 1
    f.close()
    return suffix
        
    
if __name__ == "__main__":
    
    import sys
    sys.path.append("..")
    import utils
    from pathlib import Path
    
    # get arguments (sys.argv)
    args = utils.get_args(("-cutoff",),
                          ("-func_class",),
                          ("-protein",))
    
    # "unpack" arguments
    cutoff = args.cutoff
    func_class = args.func_class
    protein = args.protein
    
    # check validity of arguments
    if any(arg is None for arg in [cutoff, func_class, protein]):
        c, f, p = "cutoff", "func_class", "protein"
        e = (f"<{c}>, <{f}> and <{p}> have no default values. Please do:\n"
             f">>> python _split_holins.py -{c} <{c}> -{f} <{f}> -{p} <{p}>")
        raise ValueError(e)
    
    # remaining arguments
    ttable = "11"
    icodons = ("TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG")
    
    # initialize directory (if not already initialized)
    Path(f"{protein}s").mkdir(exist_ok=True)

    # split and translate sequences
    num_files = split_translate(cutoff=int(cutoff),
                                func_class=func_class,
                                protein=protein,
                                ttable=ttable,
                                icodons=icodons)
    
    print(f"Number of new .fasta files: {num_files}")
