# -*- coding: utf-8 -*-

from collect_sequences import CollectSequences

def get_sequences(db: str, terms: dict, negatives: bool) -> None:
    """
    Wrapper function for CollectSequences.
    
    Parameters
    ----------
    db: str
        The name of the directory where the retrieved sequences are stored
    terms: list[str]
        A list of search terms
    negatives: bool
        Whether to include the terms in <terms> or terms not in <terms>
    """
    # collect sequences
    print("---")
    for protein in terms:
        print(f"Collecting {protein!r} sequences...")
        CollectSequences(db=db,
                         cname=protein,
                         terms=terms[protein],
                         negatives=negatives).get_sequences(taxid="2731619")
        print("---")
        
def get_dna_modification(only_terms: bool = False) -> None:
    """
    Retrieves DNA sequences coding for proteins belonging to the functional class of DNA modification
    proteins. If <only_terms> is set to True, the sequences are not retrieved and a function attribute is set
    to "terms".
    
    Parameters
    ----------
    only_terms: bool
        Whether to only get the terms associated to the functional class
    """
    # terms for DNA modification proteins
    terms = {"nuclease":             ["hnh endonuclease", "exonuclease", "hnh homing endonuclease",
                                      "cas4 family exonuclease", "putative exonuclease", "putative hnh endonuclease",
                                      "endonuclease", "homing endonuclease", "nuclease", "endonuclease vii",
                                      "putative homing endonuclease", "restriction endonuclease",
                                      "pd-(d/e)xk nuclease superfamily protein", "recombination endonuclease vii",
                                      "intron associated endonuclease", "exonuclease",
                                      "hnh endonuclease bacteriophage, hnh endonuclease, dna.52a"],
             # ---
             "deaminase":            ["deoxycytidylate deaminase", "dcmp deaminase", "deaminase"],
             # ---
             "thymidylate synthase": ["thymidylate synthase"],
             # ---
             "dUTPase":              ["dutpase"],
             # ---
             "kinase":               ["polynucleotide kinase", "thymidine kinase", "nucelotide kinase", "kinase",
                                      "deoxynucleoside monophosphate kinase"],
             # ---
             "phosphoesterase":      ["metallophosphoesterase", "phosphoesterase"],
             # ---
             "reductase":            ["dihydrofolate reductase", "ribonucleotide reductase",
                                      "anaerobic ribonucleoside triphosphate reductase", "reductase",
                                      "phosphoadenosine-phosphosulfate reductase"],
             # ---
             "methylase":            ["dna methylase", "adenine-specific methyltransferase",
                                      "dna adenine methylase", "cytosine specific methyltransferase",
                                      "dna n-6-adenine-methyltransferase", "methylase"],
             # ---
             "ATPase":               ["atpase"],
             # ---
             "nucleotidohydrolase":  ["nucleotidohydrolase", "deoxyuridine 5'-triphosphate nucleotidohydrolase"],
             # ---
             "transferase":          ["nucleotidyltransferase", "glycosyltransferase", "methyltransferase",
                                      "adenine specific dna methyltransferase", "transferase"],
             # ---
             "phosphohydrolase":     ["nucleoside triphosphate pyrophosphohydrolase", "phosphohydrolase"],
             # ---
             "glutaredoxin":         ["glutaredoxin"],
             # ---
             "ribonuclease":         ["endodeoxyribonuclease i", "endodeoxyribonuclease rusa", "ribonuclease",
                                      "endodeoxyribonuclease"]}
    if not only_terms:
        # get DNA sequences coding for DNA modification proteins
        get_sequences(db=BASE_DIR+"dna_modification", terms=terms, negatives=False)
    else:
        # initialize function attribute
        get_dna_modification.terms = terms
    
def get_dna_replication(only_terms: bool = False) -> None:
    """
    Retrieves DNA sequences coding for proteins belonging to the functional class of DNA replication
    proteins. If <only_terms> is set to True, the sequences are not retrieved and a function attribute is set
    to "terms".
    
    Parameters
    ----------
    only_terms: bool
        Whether to only get the terms associated to the functional class
    """
    # terms for DNA replication proteins
    terms = {"transcription factor":   ["whib family transcription factor", "transcriptional regulator",
                                        "putative transcriptional regulator", "regulatory protein",
                                        "ecf sigma factor", "transcription factor"],
             # ---
             "DNA primase-helicase":   ["dna helicase", "dna primase/helicase", "dna primase",
                                        "putative helicase", "helicase", "dnab-like replicative helicase",
                                        "replicative dna helicase", "dnab-like dsdna helicase", "dsdna helicase",
                                        "replicative helicase", "helicase of the snf2 rad54 family", "primase"],
             # ---
             "RNA polymerase":         ["dna-directed rna polymerase", "rna polymerase sigma factor",
                                        "rna polymerase", "rna dependent rna polymerase"],
             # ---
             "DNA ligase":             ["dna ligase"],
             # ---
             "DNA polymerase":         ["dna polymerase i", "dna polymerase", "dna polymerase b"],
             # ---
             "RNA ligase":             ["rna ligase"],
             # ---
             "replication initiation": ["replication initiation protein"]}
    if not only_terms:
        # get DNA sequences coding for DNA replication proteins
        get_sequences(db=BASE_DIR+"dna_replication", terms=terms, negatives=False)
    else:
        # initialize function attribute
        get_dna_replication.terms = terms
        
def get_lysis(only_terms: bool = False) -> None:
    """
    Retrieves DNA sequences coding for proteins belonging to the functional class of lysis proteins. If
    <only_terms> is set to True, the sequences are not retrieved and a function attribute is set to "terms".
    
    Parameters
    ----------
    only_terms: bool
        Whether to only get the terms associated to the functional class
    """
    # terms for lysis proteins
    terms = {"endolysin": ["n-acetylmuramoyl-l-alanine amidase", "endolysin", "peptidoglycan hydrolase",
                           "cell wall hydrolase autolysin"],
             # ---
             "holin":     ["holin", "putative holin", "holin protein"],
             # ---
             "spanin":    ["rz-like spanin", "spanin", "i-spanin", "o-spanin", "u-spanin", "rz", "rz1"]}
    if not only_terms:
        # get DNA sequences coding for lysis proteins
        get_sequences(db=BASE_DIR+"lysis", terms=terms, negatives=False)
    else:
        # initialize function attribute
        get_lysis.terms = terms

def get_lysogeny_repressor(only_terms: bool = False) -> None:
    """
    Retrieves DNA sequences coding for proteins belonging to the functional class of lysogeny/repressor
    proteins. If <only_terms> is set to True, the sequences are not retrieved and a function attribute is
    set to "terms".
    
    Parameters
    ----------
    only_terms: bool
        Whether to only get the terms associated to the functional class
    """
    # terms for lysogeny/repressor proteins
    terms = {"integrase":     ["integrase", "tyrosine integrase"],
             # ---
             "recombinase":   ["site specific recombinase xerd", "recombinase", "recombination protein",
                               "ninb protein", "rect protein", "ning recombination protein"],
             # ---
             "repressor":     ["immunity repressor", "repressor protein ci", "repressor protein",
                               "transcriptional repressor", "repressor domain protein",
                               "sos-response transcriptional repressor", "repressor"],
             # ---
             "resolvase":     ["holliday junction resolvase", "resolvase"],
             # ---
             "transposase":   ["transposase"],
             # ---
             "antirepressor": ["antirepressor protein", "antirepressor"]}
    if not only_terms:
        # get DNA sequences coding for lysogeny/repressor proteins
        get_sequences(db=BASE_DIR+"lysogeny_repressor", terms=terms, negatives=False)
    else:
        # initialize function attribute
        get_lysogeny_repressor.terms = terms
    
def get_packaging(only_terms: bool = False) -> None:
    """
    Retrieves DNA sequences coding for proteins belonging to the functional class of packaging proteins. If
    <only_terms> is set to True, the sequences are not retrieved and a function attribute is set to "terms".
    
    Parameters
    ----------
    only_terms: bool
        Whether to only get the terms associated to the functional class
    """
    # terms for packaging proteins
    terms = {"large terminase": ["terminase large subunit", "large terminase", "large subunit terminase"],
             # ---
             "small terminase": ["terminase small subunit", "putative terminase small subunit",
                                 "small terminase"]}
    if not only_terms:
        # get DNA sequences coding for packaging proteins
        get_sequences(db=BASE_DIR+"packaging", terms=terms, negatives=False)
    else:
        # initialize function attribute
        get_packaging.terms = terms
    
def get_structural(only_terms: bool = False) -> None:
    """
    Retrieves DNA sequences coding for proteins belonging to the functional class of structural proteins. If
    <only_terms> is set to True, the sequences are not retrieved and a function attribute is set to "terms".
    
    Parameters
    ----------
    only_terms: bool
        Whether to only get the terms associated to the functional class
    """
    # terms for structural proteins
    terms = {"minor tail":   ["minor tail protein", "minor tail"],
             # ---
             "major tail":   ["major tail protein", "major tail"],
             # ---
             "portal":       ["portal protein", "putative portal protein", "portal"],
             # ---
             "minor capsid": ["minor capsid protein", "minor head protein", "minor capsid component",
                              "minor capsid protein from bacteriophage", "minor capsid"],
             # ---
             "major capsid": ["major capsid protein", "major head protein", "major capsid"],
             # ---
             "head-tail":    ["head-to-tail adaptor", "head-to-tail stopper", "tail terminator",
                              "head-tail connector protein", "tail completion protein",
                              "tail tape measure protein", "tail length tape measure protein",
                              "tail length tape-measure protein", "head-tail adaptor", "tail tube protein",
                              "distal tail protein", "internal virion protein", "head fiber protein",
                              "tail completion or neck1 protein", "head-tail adaptor ad1",
                              "putative tail component", "head completion protein", "tail tape measure",
                              "tail component", "tail collar fiber protein", "tail connector protein",
                              "putative head tail adaptor", "head to tail adaptor", "tail tubular protein",
                              "head tail connector", "head tail connector protein", "head tail", "head-tail",
                              "head tail protein", "head-tail protein"],
             # ---
             "tail fiber":   ["tail fibers protein", "tail fiber protein", "putative tail protein",
                              "putative tail fiber protein", "tail protein", "43 kda tail protein",
                              "tail fiber", "tail fibre", "tail fibre protein"],
             # ---
             "tail sheath":  ["tail sheath protein", "tail sheath", "putative phage tail sheath protein",
                              "tail shaft"],
             # ---
             "baseplate":    ["baseplate protein", "baseplate wedge subunit", "baseplate hub",
                              "baseplate assembly protein", "baseplate wedge protein", "baseplate hub subunit",
                              "baseplate j like protein", "baseplate"],
             # ---
             "neck":         ["neck protein", "type i neck protein", "neck", "pre-neck appendage"],
             # ---
             "collar":       ["upper collar protein", "lower collar protein", "collar", "collar protein"],
             # ---
             "tailspike":    ["tailspike", "tail spike", "tail-spike", "tailspike protein", "tail spike protein",
                              "tail-spike protein"]}
    if not only_terms:
        # get DNA sequences coding for structural proteins
        get_sequences(db=BASE_DIR+"structural", terms=terms, negatives=False)
    else:
        # initialize function attribute
        get_structural.terms = terms
    
#def get_hypothetical(only_terms: bool = False) -> None:
#    """
#    Retrieves DNA sequences annotated as "hypothetical protein". Note that this annotation may have two distinct
#    meanings: 1. the function coded by the sequence is not yet known; 2. the function coded by the sequence is
#    already known, but, for some reason, the sequence was not properly annotated. If <only_terms> is set to True,
#    the sequences are not retrieved and a function attribute is set to "terms".
#    
#    Parameters
#    ----------
#    only_terms: bool
#        Whether to only get the terms associated to the functional class
#    """
#    # terms for hypothetical proteins
#    terms = {"hypothetical": ["hypothetical", "hypothetical protein", "conserved hypothetical protein",
#                              "protein of unknown function", "protein of unknown function duf859",
#                              "protein of unknown function", "protein of unknown function (duf551)",
#                              "unknown", "unknown function", "protein of unknown function (duf1492)",
#                              "protein of unknown function (duf4376)", "protein of unknown function (duf1351)",
#                              "protein of unknown function (duf2634)", "protein of unknown function (duf2577)",
#                              "protein of unknown function (duf1018)", "protein of unknown function (duf2612)",
#                              "protein of unknown function (duf2829)", "protein of unknown function (duf722)",
#                              "protein of unknown function (duf4969)", "protein of unknown function duf1424",
#                              "protein of unknown function (duf669)", "protein of unknown function (duf1642)",
#                              "protein of unknown function (duf1366)"]}
#    if not only_terms:
#        # get DNA sequences coding for hypothetical proteins
#        get_sequences(db=BASE_DIR+"hypothetical", terms=terms, negatives=False)
#    else:
#        # initialize function attribute
#        get_hypothetical.terms = terms
        
def get_miscellaneous() -> None:
    """
    Retrieves DNA sequences coding for miscellaneous proteins, that is, proteins whose common name is not present
    in the lists of terms initialized in "get_dna_modification", "get_dna_replication", "get_lysis",
    "get_lysogeny_repressor", "get_packaging", "get_structural" or "get_hypothetical".
    """
    # initialize tuple containing function objects ready to be called and an empty list for storing new values
    values = []
    funcs = (get_dna_modification, get_dna_replication, get_lysis, get_lysogeny_repressor, get_packaging,
             get_structural) #, get_hypothetical)
    # populate "values" by iterating through the terms stored as an attribute of each function
    for func in funcs:
        func(only_terms=True)
        new_values = [value for protein_key in func.terms.values() for value in protein_key]
        values.extend(new_values)
    # initialize terms
    terms = {"miscellaneous": values}
    # get DNA sequences coding for miscellaneous proteins
    get_sequences(db=BASE_DIR, terms=terms, negatives=True)
         
    
if __name__ == "__main__":
    
    import utils
    from pathlib import Path
    
    # get command line arguments
    args = utils.get_args(("-base_dir",),
                          ("-func_class",))
    base_dir = args.base_dir
    class_ = args.func_class
    
    # check whether <base_dir> and <class_> are valid (1)
    if base_dir is None or class_ is None:
        e = ("<base_dir> and <func_class> have no default values. Please do:\n"
             ">>> python _collect.py -base_dir <base_dir> -func_class <func_class>")
        raise Exception(e)
    
    # check whether <base_dir> is valid (2)
    if base_dir not in {"init", "cs"}:
        raise ValueError(f"{base_dir!r} is not valid for 'base_dir'. Choose one of {{'init', 'cs'}}.")
    
    # check whether <class_> is valid (2)
    options = ["dna-modification", "dna-replication", "lysis", "lysogeny-repressor", "packaging",
               "structural", "hypothetical", "miscellaneous"]
    if class_ not in options:
        raise ValueError(f"'{class_}' is not a valid class. Choose one of {{{', '.join(options)}}}.")
        
    # initialize BASE_DIR and create directory (if not already created)
    BASE_DIR = "../sequences/" if base_dir == "init" else "../_miscellaneous"
    Path(BASE_DIR).mkdir(exist_ok=True)
            
    # initialize function dictionary
    funcs = {"dna-modification": get_dna_modification,
             "dna-replication": get_dna_replication,
             "lysis": get_lysis,
             "lysogeny-repressor": get_lysogeny_repressor,
             "packaging": get_packaging,
             "structural": get_structural,
             "miscellaneous": get_miscellaneous}
    
    # call appropriate function
    get_func_class = funcs[class_]
    get_func_class()
        