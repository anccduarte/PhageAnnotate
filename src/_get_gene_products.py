# -*- coding: utf-8 -*-

if __name__ == "__main__":
    
    import argparse
    import os
    from gene_products import GeneProducts
    from pathlib import Path
    
    # make "gene_products" directory if it doesn't exist already
    Path("../gene_products").mkdir(exist_ok=True)

    # initialize parser
    parser = argparse.ArgumentParser()

    # add arguments
    parser.add_argument("-taxid")
    parser.add_argument("-fname", default="test")
    parser.add_argument("-num_ids", default="10")

    # read arguments
    args = parser.parse_args()
    taxid = args.taxid
    fname = args.fname
    num_ids = args.num_ids

    # check whether <taxid> was provided to sys.argv
    # for phages use either "txid28883[ORGN]" or "txid2731619[ORGN]"
    if taxid is None:
        raise Exception("'taxid' has no default value. Please, do:\n"
                        ">>> python gene_products.py -taxid <taxid>")

    # check validity of "taxid"
    if not taxid.isdigit():
        raise ValueError(f"'{taxid}' is not a valid number.")

    # check validity of "num_ids"
    if not num_ids.isdigit():
        raise ValueError(f"'{num_ids}' is not a valid number.")

    # avoid overwriting existing file
    i = 0
    while os.path.exists(f"../gene_products/{fname}{('_'+str(i))*bool(i)}.txt"):
        i += 1
    fname = fname + ("_" + str(i)) * bool(i)

    # initialize instance of GeneProducts and build .txt file
    gp = GeneProducts(taxid=taxid, fname=fname, num_ids=int(num_ids))
    num_gps = gp.build_txt()
    print(f"Number of gene products in '{fname}.txt' ({num_ids} IDs): {num_gps}")
