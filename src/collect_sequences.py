# -*- coding: utf-8 -*-

import os
import warnings
from Bio import SeqIO, SeqRecord, SeqFeature
from tqdm.auto import tqdm

# ignore Byopython warnings (SeqFeature related)
warnings.filterwarnings("ignore")

class CollectSequences:
    
    """
    Allows for the collection of DNA sequences in NCBI's Nucleotide database according
    to a set of search terms provided by the user.
    """
    
    def __init__(self, db: str, cname: str, terms: list[str]) -> None:
        """
        Initializes an instance of SeqMining.
        
        Parameters
        ----------
        db: str
            The name of the directory where the retrieved sequences are to be stored
        cname: str
            The common name of the target gene product
        terms: list[str]
            A list of search terms
        """
        self.db = CollectSequences._check_db(db)
        self.cname = "_".join("_".join(cname.split("-")).split(" "))
        self.terms = [" ".join(term.split()).lower() for term in terms]
       
    def __repr__(self) -> str:
        """
        Returns the string representation of the object.
        """
        class_ = self.__class__.__name__
        return f"{class_}({self.db!r}, {self.cname!r}, {self.terms!r})"
    
    @staticmethod
    def _check_db(db: str) -> str:
        """
        Checks whether <db> exists. If not, a new directory called <db> is created.
        Returns the string it takes as argument.
        
        Parameters
        ----------
        db: str
            The name of the directory where the retrieved sequences are stored
        """
        if not os.path.exists(db): os.mkdir(db)
        return db
    
    def _validate_product(self, product: str) -> bool:
        """
        Verifies the validity of <product>. <product> is not valid if:
        - all terms provided by the user are distinct from "product"
        - OR some term in ["not", "non"] is a substring of "product"
        
        Parameters
        ----------
        product: str
            A protein product coded by a given DNA sequence
        """
        # check whether any term is self.terms is equal to "product"
        if all(term != product for term in self.terms):
            return False
        # check whether "not" or "non" is in "product"
        if any(term in product for term in ["not", "non"]):
            return False
        # checks passed -> "product" is valid
        return True
    
    @staticmethod
    def _get_sequence(record: SeqRecord, feature: SeqFeature) -> SeqRecord:
        """
        Returns the DNA sequence to be written to the file ("temp.txt").
        
        Parameters
        ----------
        record: SeqRecord
            A SeqRecord object containing a genome
        feature: SeqFeature
            A feature of <record>
        """
        # get DNA sequence from its location in the genome
        fr, to = feature.location.start, feature.location.end
        seq = record.seq[fr:to]
        # if the sequence is located in strand -1, compute its reverse complement
        if feature.location.strand == -1:
            seq = seq.reverse_complement()
        # return sequence (or its reverse complement)
        return seq
    
    def _get_fname(self, taxid: str) -> str:
        """
        Computes and returns the name of the .fasta file where the non-redundant DNA
        sequences are to be stored.
        
        Parameters
        ----------
        taxid: str
            The txid identifier of the sought-after taxa (as in NCBI's Taxonomy database)
        """
        # construct the name of the file
        fname = f"txid{taxid}_{self.cname}"
        # check whether <fname> is the name of a file in pwd (avoid collisions)
        i = 0
        while os.path.exists(f"{self.db}/{fname}{f'({i})'*bool(i)}.fasta"):
            i += 1
        # return the name of the new .fasta file
        return fname + f"({i})"*bool(i)
    
    def _filter_fasta(self, fname: str) -> tuple:
        """
        Filters the original .fasta file ("temp.fasta") and creates a new only containing
        unique sequences. Returns a tuple consisting of:
        - nf: the number of sequences contained in the original fasta file
        - f: the number of sequences contained in the new fasta file
        
        Parameters
        ----------
        fname: str
            The name to be given to the new .fasta file
        """
        # open new file for writing non-redundant DNA sequences
        with open(f"{self.db}/{fname}.fasta", "w") as file:
            # read record in "temp.fasta"
            records = SeqIO.parse(f"{self.db}/temp.fasta", format="fasta")
            # initialize other local variables
            filt = set()
            num_nf, num_f = 0, 0
            # main loop -> add non-redundant sequences to the new file
            for record in records:
                num_nf += 1
                description, seq = record.description, record.seq
                if seq not in filt:
                    num_f += 1
                    filt.add(seq)
                    file.write(f">{description}\n{seq}\n\n")
        # return number of sequences contained in both .fasta files
        return num_nf, num_f
    
    def get_sequences(self, taxid: str) -> None:
        """
        Creates a .fasta file ("temp.fasta") containing DNA sequences related to the search
        terms provided by the user. A second .fasta file only containing non-redundant
        sequences is created by applying the auxiliary method '_filter_fasta'. Upon filtering,
        "temp.fasta" is deleted.
        
        Parameters
        ----------
        taxid: str
            The txid identifier of the sought-after taxa (as in NCBI's Taxonomy database;
            only needed to name the .fasta file)
        """
        # create and open fasta file to store DNA sequences
        with open(f"{self.db}/temp.fasta", "w") as out_file:
            # inspect the previously collected IDs
            for in_file in tqdm(os.scandir("../records")):
                # guard clause 1 -> ignore files not ending in .gb
                if not in_file.name.endswith(".gb"):
                    continue
                # read record whose ID is <id_>
                record = SeqIO.read(f"../records/{in_file.name}", "gb")
                for feature in record.features:
                    # guard clause 2 -> ignore features whose type is not "CDS"
                    if feature.type != "CDS":
                        continue
                    product = feature.qualifiers.get("product", ["not"])[0]
                    # guard clause 3 -> ignore feature if "product" is not valid
                    if not self._validate_product(product.lower()):
                        continue
                    # finally, add DNA sequence to the fasta file (all guard clauses avoided)
                    seq = CollectSequences._get_sequence(record=record, feature=feature)
                    id_, source = record.id, record.annotations['source']
                    out_file.write(f"> {id_} | {source} | {product}\n{seq}\n\n")
        # create new .fasta file containing non-redundant DNA sequences
        fname = self._get_fname(taxid=taxid)
        nf, f = self._filter_fasta(fname=fname)
        # delete original file
        os.remove(f"{self.db}/temp.fasta")
        # provide information on the collected sequences
        print(f"Filtering removed {nf-f} sequences ({f} sequences remaining in '{fname}.fasta')")
    