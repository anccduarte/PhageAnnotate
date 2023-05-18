# -*- coding: utf-8 -*-

import os
import utils
import warnings
from Bio import SeqIO, SeqRecord, SeqFeature
from tqdm.auto import tqdm

# ignore Byopython warnings (SeqFeature related)
warnings.filterwarnings("ignore")

class SeqMining:
    
    """
    Allows for the collection of DNA sequences in NCBI's Nucleotide database according to a set of
    search terms provided by the user.
    """
    
    def __init__(self, db: str, taxid: str, cname: str, terms: list[str], num_ids: int) -> None:
        """
        Initializes an instance of SeqMining.
        
        Parameters
        ----------
        db: str
            The name of the directory where the retrieved sequences are stored
        taxid: str
            The txid identifier of the sought-after taxa (as in NCBI's Taxonomy database)
        cname: str
            The common name of the target gene product
        terms: list[str]
            A list of search terms
        num_ids: int
            The number of IDs to be inspected in the database
        """
        self.db = SeqMining._check_db(db)
        self.taxid = taxid
        self.cname = "_".join("_".join(cname.split("-")).split(" "))
        self.terms = [" ".join(term.split()).lower() for term in terms]
        self.num_ids = num_ids
        
    def __repr__(self) -> str:
        """
        Returns the string representation of the object.
        """
        class_ = self.__class__.__name__
        return f"{class_}({self.db!r}, {self.taxid!r}, {self.cname!r}, {self.terms!r}, {self.num_ids})"
            
    @staticmethod
    def _check_db(db: str) -> str:
        """
        Checks whether <db> is a subdirectory of the current directory. If not, a new subdirectory
        named <db> is created. Returns the string it takes as argument.
        
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
        - all terms provided by the user cannot be found in "product"
        - OR some term in ["not", "non"] is a substring of "product" ("putative" and "possible"?)
        
        Parameters
        ----------
        product: str
            A protein product coded by a given DNA sequence
        """
        # check whether any term is self.terms is in "product"
        if all(term not in product for term in self.terms):
            return False
        # check whether "not" or "non" is in "product" (include "putative" and "possible"?)
        if any(term in product for term in ["not", "non"]):
            return False
        # "product" is valid
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
    
    def _get_fname(self) -> str:
        """
        Computes and returns the name of the .fasta file where the non-redundant DNA sequences
        are to be stored.
        """
        # construct the name of the file
        fname = f"txid{self.taxid}_{self.cname}_{self.num_ids}"
        # check whether <fname> is the name of a file in pwd (avoid collisions)
        i = 0
        while os.path.exists(f"{self.db}/{fname}{f'({i})'*bool(i)}.fasta"):
            i += 1
        # return the name of the new .fasta file
        return fname + f"({i})"*bool(i)
    
    def _filter_fasta(self, fname: str) -> tuple:
        """
        Filters the original .fasta file (generated during the search in NCBI) and creates a new
        file not containing redundant sequences. Returns a tuple consisting of:
        - nf: the number of sequences contained in the original fasta file
        - f: the number of sequences contained in the new fasta file (attained after filtering)
        
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
    
    def get_sequences(self) -> None:
        """
        Creates a .fasta file ("temp.fasta") containing DNA sequences related to the search terms
        provided by the user. A second .fasta file only containing non-redundant sequences is created
        by applying the auxiliary method '_filter_fasta'. Upon filtering, "temp.fasta" is deleted.
        """
        # retrieve IDs from NCBI's Nucleotide
        search = f"txid{self.taxid}[ORGN] ({' OR '.join(self.terms)})"
        id_list = utils.get_ids(search=search, num_ids=self.num_ids)
        # create and open fasta file to store DNA sequences
        with open(f"{self.db}/temp.fasta", "w") as file:
            # inspect the previously collected IDs
            for id_ in tqdm(id_list):
                # read record whose ID is <id_>
                record = utils.read_record(id_=id_, max_tries=5)
                for feature in record.features:
                    # guard clause 1 -> ignore features whose type is not "CDS"
                    if feature.type != "CDS":
                        continue
                    product = feature.qualifiers.get("product", ["not"])[0]
                    # guard clause 2 -> ignore feature if "product" is not valid
                    if not self._validate_product(product.lower()):
                        continue
                    # finally, add DNA sequence to the fasta file (all guard clauses avoided)
                    seq = SeqMining._get_sequence(record=record, feature=feature)
                    file.write(f"> {record.id} | {record.annotations['source']} | {product}\n{seq}\n\n")
        # create new .fasta file containing non-redundant DNA sequences
        fname = self._get_fname()
        nf, f = self._filter_fasta(fname)
        # delete original file
        os.remove(f"{self.db}/temp.fasta")
        # provide information on the collected sequences
        print(f"Filtering removed {nf-f} sequences ({f} sequences remaining in '{fname}.fasta')")
