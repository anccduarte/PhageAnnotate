# -*- coding: utf-8 -*-

from pathlib import Path
from sequence_mining import SeqMining

BASE_DIR = "../sequences/"
Path(BASE_DIR).mkdir(exist_ok=True)

def get_sequences(db: str, terms: list[str], num_ids: int) -> None:
    """
    Wrapper function for SeqMining.
    
    Parameters
    ----------
    db: str
        The name of the directory where the retrieved sequences are stored
    terms: list[str]
        A list of search terms
    num_ids: int
        The number of sequence IDs to be collected for inspection
    negatives: bool
        Whether to collect positive or negative sequences (with respect to "terms")
    """
    # collect sequences
    print("---")
    for protein in terms:
        sm = SeqMining(db=db, taxid="2731619", terms=terms[protein], num_ids=num_ids)
        sm.get_sequences()
        print("---")
        
def get_structural(num_ids: int) -> None:
    """
    Attempt to replicate PhANNs' database. Almost all search terms match the ones used in the original
    work. However, the terms "putative" and "possible" are allowed when retrieving positive sequences
    from NCBI.
    
    Parameters
    ----------
    num_ids: int
        The number of sequence IDs to be collected for inspection
    """
    # terms used by the PhANNs team + "tail spike"
    terms = {"major capsid": ["major capsid"],
             "minor capsid": ["minor capsid"],
             "baseplate": ["baseplate"],
             "major tail": ["major tail"],
             "minor tail": ["minor tail"],
             "portal": ["portal"],
             "tail fiber": ["tail fiber", "tail fibre"],
             "tail spike": ["tail spike", "tailspike"],
             "collar": ["collar"],
             "shaft": ["tail shaft", "tail sheath"],
             "head-tail": ["head-tail"]}
    # get DNA sequences coding for structural proteins
    get_sequences(db=BASE_DIR+"structural", terms=terms, num_ids=num_ids)
    
def get_lysis(num_ids: int) -> None:
    """
    Retrieves DNA sequences coding for proteins belonging to the functional class of lysis proteins.
    Such proteins include endolysins, holins and spanins.
    
    Parameters
    ----------
    num_ids: int
        The number of sequence IDs to be collected for inspection
    """
    # terms for lysis proteins
    terms = {"endolysin": ["endolysin"],
             "holin": ["holin"],
             "spanin": ["spanin", "i-spanin", "o-spanin", "u-spanin", "Rz", "Rz1"]}
    # get DNA sequences coding for lysis proteins
    get_sequences(db=BASE_DIR+"lysis", terms=terms, num_ids=num_ids)

def get_modification_replication(num_ids: int) -> None:
    """
    Retrieves DNA sequences coding for proteins belonging to the functional class of DNA modification
    and replication proteins. Such proteins include, for example, polymerases and helicases.
    
    Parameters
    ----------
    num_ids: int
        The number of sequence IDs to be collected for inspection
    """
    # terms for DNA modification and replication proteins
    terms = {"DNA polymerase": ["DNA polymerase"],
             "RNA polymerase": ["RNA polymerase"],
             "ligase": ["ligase"],
             "helicase": ["helicase"],
             "primase": ["primase"],
             "topoisomerase": ["topoisomerase"],
             "DNA methylase": ["DNA methylase"],
             "nuclease": ["nuclease", "exonuclease", "endonuclease"],
             "recombinase": ["recombinase", "integrase"],
             "thymidylate synthase": ["thymidylate synthase"],
             "ribonucleotide reductase": ["ribonucleotide reductase"],
             "dihydrofolate reductase": ["dihydrofolate reductase"],
             "thymidine kinase": ["thymidine kinase"],
             "deoxynucleoside monophosphate kinase": ["deoxynucleoside monophosphate kinase"],
             "dCMP deaminase": ["dCMP deaminase", "deaminase"]}
    # get DNA sequences coding for DNA modification and replication proteins
    get_sequences(db=BASE_DIR+"modification_replication", terms=terms, num_ids=num_ids)
    
def get_packaging(num_ids: int) -> None:
    """
    Retrieves DNA sequences coding for proteins belonging to the functional class of DNA packaging.
    This functional class only includes small and large terminases.
    
    Parameters
    ----------
    num_ids: int
        The number of sequence IDs to be collected for inspection
    """
    # terms for DNA packaging
    terms = {"small terminase": ["small terminase", "terminase small subunit"],
             "large terminase": ["large terminase", "terminase large subunit"]}
    # get DNA sequences coding for DNA packaging proteins
    get_sequences(db=BASE_DIR+"packaging", terms=terms, num_ids=num_ids)
        
    
if __name__ == "__main__":
    
    import argparse
    
    # initialize parser
    parser = argparse.ArgumentParser()
    
    # add arguments
    parser.add_argument("-func_class")
    parser.add_argument("-num_ids")
    
    # read arguments
    args = parser.parse_args()
    class_ = args.func_class
    num_ids = args.num_ids
    
    # check whether <func_class> and <num_ids> were provided to sys.argv
    if class_ is None or num_ids is None:
        raise Exception("'func_class' and 'num_ids' have no default values. Please, do:\n"
                        ">>> python collect_seqs.py -func_class <func_class> -num_ids <num_ids>")
    
    # check validity of "num_ids"
    if not num_ids.isdigit():
        raise ValueError(f"'{num_ids}' is not a valid number.")
    
    # check validity of "class_"
    options = ["structural", "lysis", "modification-replication", "packaging"]
    if class_ not in options:
        raise ValueError(f"'{class_}' is not a valid class. Choose one of {{{', '.join(options)}}}.")
    
    # build dictionary of functions (each representing a functional class)
    funcs = {"structural": get_structural,
             "lysis": get_lysis,
             "modification-replication": get_modification_replication,
             "packaging": get_packaging}
    
    # call the appropriate function
    f = funcs[class_]
    f(num_ids=int(num_ids))
    