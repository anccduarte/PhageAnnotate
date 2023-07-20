# -*- coding: utf-8 -*-

def _get_count(entry: str) -> int:
    """
    Gets and returns the count associated to entry it takes as input.
    
    Parameters
    ----------
    entry: str
        An entry consisting of a gene product and the respective count
    """
    *_, last = entry.split(":")
    return int(last[1:-1])

def filter_file(nfile: str, tol: int) -> tuple:
    """
    Filters the .txt file it takes as input and creates a new .txt file
    only containing the names of the gene products whose count is, at
    least, equal to <min_occurs>.
    
    Parameters
    ----------
    nfile: str
        The name of the .txt file containing the gene products' names
    tol: int
        The minimum number of occurrences allowed
    """
    # read the contents of the file
    with open(f"../../gene_products/{nfile}") as f:
        entries = f.readlines()
    # open new file for writing the desired entries
    name, _ = nfile.split(".")
    with open(f"../../gene_products/{name}_tol{tol}.txt", "w") as fout:
        # initialize counter
        num_good = 0
        # iterate through the entries
        for entry in entries:
            count = _get_count(entry)
            if count < tol:
                break
            fout.write(entry)
            num_good += 1
    # return tuple containing the number of entries in each file
    return (len(entries), num_good)
    

if __name__ == "__main__":
    
    import sys
    sys.path.append("..")
    import utils
    
    args = utils.get_args(("-nfile",),
                          ("-tol",))
    
    nfile = args.nfile
    tol = args.tol
    
    if nfile is None or tol is None:
        e = ("'nfile' and 'tol' have no default values. Please, do:\n"
             ">>> python _filter_products.py -nfile <nfile> -tol <tol>")
        raise Exception(e)
    
    ninit, nfilt = filter_file(nfile=nfile, tol=int(tol))
    
    print(f"Number of entries in original file: {ninit}\n"
          f"Number of entries after filtering: {nfilt}")
