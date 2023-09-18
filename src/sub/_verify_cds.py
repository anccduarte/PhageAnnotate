# -*- coding: utf-8 -*-

from Bio import SeqIO

def verify_all_cds_from_file(func_class: str, protein: str) -> tuple:
    """
    Verifies the validity of all the coding regions (CDS) present in the file
    pointed to by "../../sequences/<func_class>/txid2731619_<protein>.fasta" by
    printing to screen some useful information about a CDS whenever the length
    of the corresponding sequence is not divisible by 3:
    - ---
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
    protein = "_".join(protein.split())
    path = f"../../sequences/{func_class}/txid2731619_{protein}.fasta"
    records = SeqIO.parse(path, "fasta")
    counter = 0
    print("---")
    for i, record in enumerate(records):
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
    return counter, i+1


if __name__ == "__main__":
    
    """
    Note: The module is used for validating the quality of the sequences initially
    retrieved from NCBI's Nucleotide (via executing "_collect.py" with "-base_dir"
    set to "init"). It checks, for all sequences in a given .fasta file, the length
    of its sequences; if the length of a sequence is not divisible by 3 (meaning
    that there is something wrong with the sequence, since codons are constituted
    by 3 nucleotides), it prints its length, its length modulo 3 (to check how many
    nucleotides are in excess), its first codon and its last codon (if neither
    corresponds to an init codon or stop codon, respectively, then something went
    really bad upon the registering of the sequence in NCBI).
    ---
    Note 2: The module is only used to acess the quality of the sequence database
    in a very embrionary stage of the pipeline of PhageAnnotate and is not of vital
    importance for the tool's design. However, valuable insights may me drawn from
    it and hence the decision for its presence in the repository.
    """
        
    import sys
    sys.path.append("..")
    import utils
    
    args = utils.get_args(("-func_class",),
                          ("-protein",))
    
    func_class = args.func_class
    protein = args.protein
    
    if func_class is None or protein is None:
        e = ("<func_class> and <file> have no default values. Please do:\n"
             ">>> python _verify_cds.py -func_class <func_class> -protein <protein>")
        raise ValueError(e)
        
    not_valid, total = verify_all_cds_from_file(func_class, protein)
    perc = not_valid / total
    print(f"Number of non valid CDS regions for {protein!r}: {not_valid}")
    print(f"Percentage of non valid CDS regions for {protein!r}: {perc:.2%}")
    print("---")
    