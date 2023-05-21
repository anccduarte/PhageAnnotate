# -*- coding: utf-8 -*-

from pathlib import Path
from collect_sequences import CollectSequences

BASE_DIR = "../sequences/"
Path(BASE_DIR).mkdir(exist_ok=True)

def get_sequences(db: str, terms: dict) -> None:
    """
    Wrapper function for CollectSequences.
    
    Parameters
    ----------
    db: str
        The name of the directory where the retrieved sequences are stored
    terms: list[str]
        A list of search terms
    """
    # collect sequences
    print("---")
    for protein in terms:
        print(f"Collecting {protein!r} sequences...")
        CollectSequences(db=db,
                         cname=protein,
                         terms=terms[protein]).get_sequences(taxid="2731619")
        print("---")
        
def get_dna_modification() -> None:
    """
    Retrieves DNA sequences coding for proteins belonging to the functional class of DNA modification
    proteins.
    """
    # terms for DNA modification proteins
    terms = {"nuclease": ["hnh endonuclease", "exonuclease", "hnh homing endonuclease",
                          "cas4 family exonuclease", "putative exonuclease", "putative hnh endonuclease",
                          "endonuclease", "homing endonuclease", "nuclease", "endonuclease vii",
                          "putative homing endonuclease", "restriction endonuclease",
                          "pd-(d/e)xk nuclease superfamily protein", "recombination endonuclease vii",
                          "intron associated endonuclease", "exonuclease",
                          "hnh endonuclease bacteriophage, hnh endonuclease, dna.52a"],
             "deaminase": ["deoxycytidylate deaminase", "dcmp deaminase", "deaminase"],
             "thymidylate synthase": ["thymidylate synthase"],
             "dUTPase": ["dutpase"],
             "kinase": ["polynucleotide kinase", "thymidine kinase", "deoxynucleoside monophosphate kinase",
                        "nucelotide kinase", "kinase"],
             "phosphoesterase": ["metallophosphoesterase", "phosphoesterase"],
             "reductase": ["dihydrofolate reductase", "ribonucleotide reductase",
                           "anaerobic ribonucleoside triphosphate reductase",
                           "phosphoadenosine-phosphosulfate reductase", "reductase"],
             "methylase": ["dna methylase", "adenine-specific methyltransferase", "dna adenine methylase",
                           "cytosine specific methyltransferase", "dna n-6-adenine-methyltransferase",
                           "methylase"],
             "ATPase": ["atpase"],
             "nucleotidohydrolase": ["nucleotidohydrolase", "deoxyuridine 5'-triphosphate nucleotidohydrolase"],
             "transferase": ["nucleotidyltransferase", "glycosyltransferase", "methyltransferase",
                             "adenine specific dna methyltransferase", "transferase"],
             "phosphohydrolase": ["nucleoside triphosphate pyrophosphohydrolase", "phosphohydrolase"],
             "glutaredoxin": ["glutaredoxin"],
             "chromatin remodeling complex atpase": ["chromatin remodeling complex atpase"],
             "ribonuclease": ["endodeoxyribonuclease i", "endodeoxyribonuclease rusa", "endodeoxyribonuclease",
                              "ribonuclease"],
             "restriction alleviation protein": ["restriction alleviation protein"]}
    # get DNA sequences coding for DNA modification proteins
    get_sequences(db=BASE_DIR+"dna_modification", terms=terms)
    
def get_dna_replication() -> None:
    """
    Retrieves DNA sequences coding for proteins belonging to the functional class of DNA replication
    proteins.
    """
    # terms for DNA replication proteins
    terms = {"transcription factor": ["whib family transcription factor", "transcriptional regulator",
                                      "putative transcriptional regulator", "regulatory protein",
                                      "ecf sigma factor", "transcription factor"],
             "DNA primase-helicase": ["dna helicase", "dna primase/helicase", "dna primase",
                                      "putative helicase", "helicase", "dnab-like replicative helicase",
                                      "replicative dna helicase", "dnab-like dsdna helicase", "dsdna helicase",
                                      "replicative helicase", "helicase of the snf2 rad54 family", "primase"],
             "RNA polymerase": ["dna-directed rna polymerase", "rna polymerase sigma factor", "rna polymerase",
                                "rna dependent rna polymerase"],
             "DNA ligase": ["dna ligase"],
             "DNA polymerase": ["dna polymerase i", "dna polymerase", "dna polymerase b"],
             "RNA ligase": ["rna ligase"],
             "replication initiation protein": ["replication initiation protein"],
             "replisome organizer": ["replisome organizer"],
             "ParB protein": ["parb protein"],
             "chromosome partitioning protein": ["chromosome partitioning protein"]}
    # get DNA sequences coding for DNA replication proteins
    get_sequences(db=BASE_DIR+"dna_replication", terms=terms)
        
def get_lysis() -> None:
    """
    Retrieves DNA sequences coding for proteins belonging to the functional class of lysis proteins.
    """
    # terms for lysis proteins
    terms = {"endolysin": ["n-acetylmuramoyl-l-alanine amidase", "endolysin", "peptidoglycan hydrolase",
                           "cell wall hydrolase autolysin"],
             "holin": ["holin", "putative holin", "holin protein"],
             "spanin": ["rz-like spanin", "spanin", "i-spanin", "o-spanin", "u-spanin", "Rz", "Rz1"]}
    # get DNA sequences coding for lysis proteins
    get_sequences(db=BASE_DIR+"lysis", terms=terms)

def get_lysogeny_repressor() -> None:
    """
    Retrieves DNA sequences coding for proteins belonging to the functional class of lysogeny/repressor
    proteins.
    """
    # terms for lysogeny/repressor proteins
    terms = {"integrase": ["integrase", "tyrosine integrase"],
             "recombinase": ["site specific recombinase xerd", "recombinase", "recombination protein",
                             "ninb protein", "rect protein", "ning recombination protein"],
             "repressor": ["immunity repressor", "repressor protein ci", "repressor protein",
                           "transcriptional repressor", "repressor domain protein",
                           "sos-response transcriptional repressor", "repressor"],
             "resolvase": ["holliday junction resolvase", "resolvase"],
             "transposase": ["transposase"],
             "antirepressor": ["antirepressor protein", "antirepressor"],
             "excisionase": ["excisionase"]}
    # get DNA sequences coding for lysogeny/repressor proteins
    get_sequences(db=BASE_DIR+"lysogeny_repressor", terms=terms)
    
def get_packaging() -> None:
    """
    Retrieves DNA sequences coding for proteins belonging to the functional class of packaging proteins.
    """
    # terms for packaging proteins
    terms = {"large terminase": ["terminase large subunit", "large terminase", "large subunit terminase"],
             "small terminase": ["terminase small subunit", "putative terminase small subunit",
                                 "small terminase"]}
    # get DNA sequences coding for packaging proteins
    get_sequences(db=BASE_DIR+"packaging", terms=terms)
    
def get_structural() -> None:
    """
    Retrieves DNA sequences coding for proteins belonging to the functional class of structural proteins.
    """
    # terms for structural proteins
    terms = {"minor tail": ["minor tail protein", "minor tail"],
             "major tail": ["major tail protein", "major tail"],
             "portal": ["portal protein", "putative portal protein", "portal"],
             "minor capsid": ["minor capsid protein", "minor head protein", "minor capsid component",
                              "minor capsid protein from bacteriophage", "minor capsid"],
             "major capsid": ["major capsid protein", "major head protein", "major capsid"],
             "head-tail": ["head-to-tail adaptor", "head-to-tail stopper", "tail terminator",
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
             "tail fiber": ["tail fibers protein", "tail fiber protein", "putative tail protein",
                            "putative tail fiber protein", "tail protein", "43 kda tail protein",
                            "tail fiber", "tail fibre", "tail fibre protein"],
             "tail sheath": ["tail sheath protein", "tail sheath", "putative phage tail sheath protein",
                             "tail shaft"],
             "baseplate": ["baseplate protein", "baseplate wedge subunit", "baseplate hub",
                           "baseplate assembly protein", "baseplate wedge protein", "baseplate hub subunit",
                           "baseplate j like protein", "baseplate"],
             "neck": ["neck protein", "type i neck protein", "neck", "pre-neck appendage"],
             "collar": ["upper collar protein", "lower collar protein", "collar", "collar protein"],
             "tailspike": ["tailspike", "tail spike", "tail-spike", "tailspike protein", "tail spike protein",
                           "tail-spike protein"]}
    # get DNA sequences coding for structural proteins
    get_sequences(db=BASE_DIR+"structural", terms=terms)

def get_other() -> None:
    """
    Retrieves DNA sequences coding for proteins belonging to the functional class of "other".
    """
    # terms for "other" proteins
    terms = {"scaffolding protein": ["scaffolding protein", "capsid and scaffold protein",
                                     "head scaffolding protein", "capsid scaffolding protein"],
             "assembly protein": ["tail assembly chaperone", "tail assembly protein", "tail fiber assembly protein",
                                  "tail assembly chaperone protein"],
             "thioredoxin": ["thioredoxin"],
             "prohead protease": ["capsid maturation protease", "prohead protease", "head maturation protease",
                                  "prohead serine protease"],
             "virion morphogenesis protein": ["virion morphogenesis protein"],
             "AAA domain protein": ["aaa domain protein"],
             "helix-destabilizing protein": ["helix-destabilizing protein", "dna helix destabilizing protein"],
             "antitoxin": ["antitoxin"],
             "hemolysin": ["hemolysin"],
             "YopX protein": ["yopx protein"],
             "oxygenase": ["oxygenase", "2og-fe(ii) oxygenase"],
             "amidoligase": ["amidoligase", "putative amidoligase enzyme"],
             "CLP protease": ["clp protease", "putative atp dependent clp protease"],
             "pyocin activator": ["pyocin activator", "pyocin activator protein prtn"]}
    # get DNA sequences coding for "other" proteins
    get_sequences(db=BASE_DIR+"other", terms=terms)
        
    
if __name__ == "__main__":
    
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-func_class")
    args = parser.parse_args()
    class_ = args.func_class
    
    if class_ is None:
        raise Exception("'func_class' has no default value. Please, do:\n"
                        ">>> python _collect.py -func_class <func_class>")
    
    options = ["dna-modification", "dna-replication", "lysis", "lysogeny-repressor", "packaging",
               "structural", "other"]
    if class_ not in options:
        raise ValueError(f"'{class_}' is not a valid class. Choose one of {{{', '.join(options)}}}.")
    
    funcs = {"dna-modification": get_dna_modification,
             "dna-replication": get_dna_replication,
             "lysis": get_lysis,
             "lysogeny-repressor": get_lysogeny_repressor,
             "packaging": get_packaging,
             "structural": get_structural,
             "other": get_other}
    
    get_func_class = funcs[class_]
    get_func_class()
        