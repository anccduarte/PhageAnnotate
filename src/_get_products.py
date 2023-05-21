# -*- coding: utf-8 -*-

import os
from Bio import SeqIO
from tqdm.auto import tqdm

def get_gene_products(path: str) -> dict:
    """
    Returns a dictionary of gene products containing unique names of gene products
    and the respective counts, constructed from a set of .gb files contained in a
    directory pointed to by <path>.
    
    Parameters
    ----------
    path: str
        The path to the directory containing the .gb files
    """
    # initialize dictionary for storing gene products
    gene_products = {}
    # iterate through files in the directory pointed to by <path>
    for file in tqdm(os.scandir(path)):
        # guard clause 1 -> ignore files not ending in ".gb"
        if not file.name.endswith(".gb"):
            continue
        # read record with name <file.name>
        record = SeqIO.read(f"{path}/{file.name}", "gb")
        # iterate through features of <record>
        for feature in record.features:
            # guard clause 2 -> ignore non-coding sequences
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
    # return dictionary of gene products
    return gene_products

def build_txt(fname: str) -> int:
    """
    Gets a dictionary of gene products by calling 'get_gene_products' and builds a
    .txt file with these products and the respective counts. Returns the number of
    entries in the dictionary of gene products.
    
    Parameters
    ----------
    fname: str
        The name tobe given to the .txt file
    """
    # get gene products
    gene_products = get_gene_products(path="../records")
    # write gene products and counts to .txt file
    with open(f"../gene_products/{fname}.txt", "w") as f:
        for i, gp in enumerate(gene_products):
            f.write(f"- {gp}: {gene_products[gp]}\n")
    # return number of gene products
    return i+1


if __name__ == "__main__":
    
    import argparse
    from pathlib import Path
    
    Path("../gene_products").mkdir(exist_ok=True)
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-fname", default="test")
    args = parser.parse_args()
    fname = args.fname
        
    num = build_txt(fname=fname)
    print(f"Number of distinct names: {num}")
    