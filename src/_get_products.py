# -*- coding: utf-8 -*-

import os
from Bio import SeqIO, SeqRecord
from tqdm.auto import tqdm

def _update_gene_products(record: SeqRecord, gene_products: dict) -> dict:
    """
    Updates the dictionary <gene_products> by inspecting the gene products
    present in the genome record <record>.
    
    Parameters
    ----------
    record: SeqRecord
        A genome record
    gene_products: dict
        A dictionary of gene products
    """
    # iterate through features of <record>
    for feature in record.features:
        # guard clause -> ignore non-coding sequences
        if feature.type != "CDS":
            continue
        # get gene product (if non-existent -> "hypothetical protein")
        product = feature.qualifiers.get("product",
                                         ["hypothetical protein"])[0].lower()
        # update dictionary of gene products (add/update entry)
        if product not in gene_products:
            gene_products[product] = 1
        else:
            gene_products[product] += 1
    # return updated dictionary
    return gene_products

def _rank_gene_products(gene_products: dict) -> list:
    """
    Creates a list object by sorting the entries of the dictionary it takes as
    input in descending order of its values.
    
    Parameters
    ----------
    gene_products: dict
        A dictionary of gene products and the respective counts
    """
    # sort and return entries of dictionary
    items = gene_products.items()
    return sorted(items, key=lambda x: -x[1])

def _get_ranked_gene_products(path: str) -> list:
    """
    Returns a list of tuples (gene_product, count) containing unique names of gene
    products sorted by descending order of number of occurrences, constructed from
    a set of .gb files contained in a directory pointed to by <path>.
    
    Parameters
    ----------
    path: str
        The path to the directory containing the .gb files
    """
    # initialize dictionary for storing gene products
    gene_products = {}
    # iterate through files in the directory pointed to by <path>
    for file in tqdm(os.scandir(path)):
        # guard clause -> ignore files not ending in ".gb"
        if not file.name.endswith(".gb"):
            continue
        # read record with name <file.name>
        record = SeqIO.read(f"{path}/{file.name}", "gb")
        # update dictionary of gene products
        gene_products = _update_gene_products(record, gene_products)
    # rank gene products by the number of occurrences (descending order)
    ranked_gp = _rank_gene_products(gene_products=gene_products)
    # return list of ranked gene products
    return ranked_gp

def build_gene_products_txt(fname: str) -> int:
    """
    Constructs a list of gene products by calling 'get_gene_products' and builds a
    .txt file with these products and the respective counts. Returns the number of
    entries in the list of gene products.
    
    Parameters
    ----------
    fname: str
        The name tobe given to the .txt file
    """
    # get gene products
    ranked_gp = _get_ranked_gene_products(path="../records")
    # write gene products and counts to .txt file
    with open(f"../gene_products/{fname}.txt", "w") as fout:
        for i, (gp, count) in enumerate(ranked_gp):
            fout.write(f"- {gp}: {count}\n")
    # return number of gene products
    return i+1


if __name__ == "__main__":
    
    import utils
    from pathlib import Path
    
    Path("../gene_products").mkdir(exist_ok=True)
    
    args = utils.get_args(("-fname",))
    fname = args.fname
    
    if fname is None:
        raise Exception("<fname> has no default value. Please do:\n"
                        ">>> python _gene_products.py -fname <fname>")
       
    num = build_gene_products_txt(fname=fname)
    print(f"Number of distinct names: {num}")
    