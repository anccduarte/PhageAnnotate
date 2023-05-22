# -*- coding: utf-8 -*-

from Bio import SeqIO

def verify_all_cds_from_file(func_class: str, protein: str) -> int:
    """
    Verifies the validity of all the coding regions (CDS) present in the file
    pointed to by "../sequences/<func_class>/txid2731619_<protein>.fasta" by
    printing to screen some useful information about a CDS whenever the length
    of the corresponding sequence is not divisible by 3:
    - first codon
    - last codon
    - length of the sequence
    - length of the sequence modulo 3
    
    Parameters
    ----------
    func_class: str
        The subdirectory where <file> is stored
    protein: str
        The name of the protein
    """
    path = f"../sequences/{func_class}/txid2731619_{protein}.fasta"
    records = SeqIO.parse(path, "fasta")
    counter = 0
    print("---")
    for record in records:
        seq = record.seq
        l = len(seq)
        if l%3 != 0:
            counter += 1
            print(record.description[1:])
            print(f"first codon: {seq[:3]}")
            print(f"last codon: {seq[-3:]}")
            print(f"len sequence: {l}")
            print(f"len sequence mod 3: {l%3}")
            print("---")
    return counter


if __name__ == "__main__":
        
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-func_class")
    parser.add_argument("-protein")
    
    args = parser.parse_args()
    func_class = args.func_class
    protein = args.protein
    
    if func_class is None or protein is None:
        e = ("<func_class> and <file> have no default values. Please do:\n"
             ">>> python _verify_cds.py -func_class <func_class> -protein <protein>")
        raise ValueError(e)
        
    counter = verify_all_cds_from_file(func_class, protein)
    print(f"Number of non valid CDS for {protein!r}: {counter}")
    print("---")
    