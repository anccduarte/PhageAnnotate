# -*- coding: utf-8 -*-

import utils
import warnings
from tqdm.auto import tqdm

# ignore Byopython warnings (SeqFeature related)
warnings.filterwarnings("ignore")

class GeneProducts:
    
    """
    Allows for the construction of a .txt document containing all the names of gene products
    associated to the taxid <taxid> specified by the user, as well as the respective counts.
    """
    
    def __init__(self, taxid: str, fname: str, num_ids: int) -> None:
        """
        Initializes an instance of GeneProducts.
        
        Parameters
        ----------
        taxid: str
            The txid identifier of the sought-after taxa (as in NCBI's Taxonomy database)
        fname: str
            The name to be given to the file containing the names of the gene products
        num_ids: int
            The number of IDs to be inspected in the database
        """
        self.taxid = taxid
        self.fname = fname
        self.num_ids = num_ids
    
    def __repr__(self) -> str:
        """
        Returns the string representation of the object.
        """
        class_ = self.__class__.__name__
        return f"{class_}({self.taxid!r}, {self.fname!r}, {self.num_ids})"
        
    def _get_gene_products(self) -> dict:
        """
        Builds and returns a dictionary of unique gene products (and the respective counts)
        associated to the taxid <taxid> provided by the user.
        """
        # retrieve IDs from NCBI's Nucleotide
        id_list = utils.get_ids(search=f"txid{self.taxid}[ORGN]", num_ids=self.num_ids)
        # initialize dictionary that will contain gene product names and counts
        gene_products = {}
        # main loop -> iterate through ids
        for id_ in tqdm(id_list):
            record = utils.read_record(id_=id_, max_tries=5)
            for feature in record.features:
                # guard clause 1 -> ignore non-CDS features
                if feature.type != "CDS":
                    continue
                # guard clause 2 -> ignore CDS not containing a product
                if "product" not in feature.qualifiers:
                    continue
                # update gene products dictionary
                product = feature.qualifiers["product"][0].lower()
                if product not in gene_products:
                    gene_products[product] = 1
                else:
                    gene_products[product] += 1
        # return dictionary of gene product names
        return gene_products

    def build_txt(self) -> int:
        """
        Builds a .txt file from the dictionary returned by _get_gene_products(). Returns
        the number of entries in the file.
        """
        # get gene products
        gene_products = self._get_gene_products()
        # open new file for writing product names
        with open(f"../gene_products/{self.fname}.txt", "w") as f:
            # main loop -> iterate through gene products
            for i, gp in enumerate(gene_products):
                f.write(f"- {gp}: {gene_products[gp]}\n")
        # return number of gene product names
        return i+1
    