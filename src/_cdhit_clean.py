
import os
import utils
from Bio import SeqIO

# GET NON-REDUNDANT
# ---

def get_non_redundant() -> None:
    """
    Runs "cd-hit-est" for all .fasta files in <COMMON>. For each input .fasta
    file, it generates a new .fasta file only containing non-redundant DNA
    sequences and a .clstr file containing information on the constructed
    clusters of sequences. [Note: all non-required parameters of "cd-hit-est"
    are left unchanged - including the sequence identity threshold (0.9)]
    """
    with utils.change_dir(CD_HIT):
        os.mkdir(RES)
        for fc in FUNC_CLASSES:
            print(f"{fc.upper()}")
            path = f"{COMMON}/{fc}"
            for file in os.scandir(path):
                if not file.name.endswith(".fasta"):
                    continue
                complete_path = f"{path}/{file.name}"
                cmd = f"cd-hit-est -i {complete_path} -o {RES}/{file.name}"
                os.system(cmd)

# DISTRIBUTE FILES
# ---

def distribute_files() -> None:
    """
    Distributes the .fasta files generated by "get_non_redundant". All files
    are placed in an appropriate directory according to the functional class
    of the sequences contained within them.
    """
    with utils.change_dir(CD_HIT):
        for fc in FUNC_CLASSES:
            os.mkdir(f"{RES}/{fc}")
            path = f"{COMMON}/{fc}"
            for file in os.scandir(path):
                if not file.name.endswith(".fasta"):
                    continue
                path_cdhit = f"{RES}/{file.name}"
                new_cdhit = f"{RES}/{fc}/{file.name}"
                os.rename(path_cdhit, new_cdhit)

# GET RESULTS
# ---

def _count_sequences(path: str):
    """
    Counts and returns the number of DNA sequences present in the .fasta file
    pointed to by <path>.

    Parameters
    ----------
    path: str
        The path to the .fasta file
    """
    count = 0
    records = SeqIO.parse(path, "fasta")
    for _ in records:
        count += 1
    return count

def get_results(name_out: str) -> tuple:
    """
    Generates a .txt file containing, for each specific function within each
    functional class, the number of sequences existing in the original file
    (generated by _collect.py) and in the filtered .fasta (generated by
    "get_non_redundant" --> "distribute_files"). Returns a tuple of length 2
    consisting of the total number of sequences present in the original and
    filtered files.

    Parameters
    ----------
    name_out: str
        The name to be given to the output file (containing the counts)
    """
    total_original = 0
    total_new = 0
    # ---
    with utils.change_dir(CD_HIT):
        # ---
        with open(f"{name_out}.txt", "w") as fout:
            # ---
            fout.write("\n")
            # ---
            for fc in FUNC_CLASSES:
                # ---
                fout.write(f"- {fc.upper()}\n")
                path = f"{COMMON}/{fc}"
                # ---
                for file in os.scandir(path):
                    if not file.name.endswith(".fasta"):
                        continue
                    # ---
                    temp = "-".join(file.name.split("_")[1:])
                    name_prot = temp.split(".")[0]
                    path_original = f"{path}/{file.name}"
                    path_new = f"{RES}/{fc}/{file.name}"
                    # ---
                    count_original = _count_sequences(path_original)
                    count_new = _count_sequences(path_new)
                    total_original += count_original
                    total_new += count_new
                    # ---
                    fout.write(f"  - {name_prot}\n")
                    fout.write(f"    - {count_original = }\n")
                    fout.write(f"    - {count_new = }\n")
                # ---
                fout.write("\n")
    # ---
    return (total_original, total_new)


if __name__ == "__main__":

    # initialize first set of global variables
    CD_HIT = "../../cd-hit"
    FUNC_CLASSES = ("dna_modification",
                    "dna_replication",
                    "lysis",
                    "lysogeny_repressor",
                    "packaging",
                    "structural")

    # get command line arguments
    args = utils.get_args(("-init",),
                          ("-name_out",))
    init = args.init
    name_out = args.name_out

    # check validity of args (1)
    if init is None or name_out is None:
        i, n = "init", "name_out"
        raise Exception(f"{i!r} and {n!r} have no default values. Please do:\n"
                        f">>> python _cdhit_clean.py -{i} <{i}> -{n} <{n}>")
        
    # check validity of args (2)
    if init not in {"yes", "no"}:
        raise ValueError(f"{init!r} is not valid for 'init'."
                         "Choose one of {{'yes', 'no'}}.")
        
    # convert <init> to bool and initialize second set of global variables
    init = (init == "yes")
    seqs_dir = "sequences" if init else "sequences_cs"
    COMMON = f"../PhageAnnotate/{seqs_dir}"
    RES = "results" if init else "results_cs"

    # run
    get_non_redundant()
    distribute_files()
    total_original, total_new = get_results(name_out=name_out)

    # display total counts
    print(f"\nCOUNTS\n{total_original = }\n{total_new = }")
