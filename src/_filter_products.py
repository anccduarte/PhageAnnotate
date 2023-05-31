# -*- coding: utf-8 -*-

def filter_file(nfile: str, min_occurs: int) -> tuple:
    """
    Filters the .txt file it takes as input and creates a new .txt file
    only containing the names of the gene products whose count is, at
    least, equal to <min_occurs>.
    
    Parameters
    ----------
    nfile: str
        The name of the .txt file containing the gene products' names
    min_occurs: int
        The minimum number of occurrences allowed
    """
    # read the contents of the file
    with open(f"../gene_products/{nfile}") as f:
        entries = f.readlines()
    # open new file for writing the desired entries
    name, _ = nfile.split(".")
    with open(f"../gene_products/{name}_min_{min_occurs}.txt", "w") as f:
        # initialize counter
        num_good = 0
        # iterate through the entries
        for i, entry in enumerate(entries):
            *_, num = entry.split(":")
            count = int(num[1:-1])
            if count < min_occurs:
                continue
            f.write(entry)
            num_good += 1
    # return tuple containing the number of entries in each file
    return (i+1, num_good)
    

if __name__ == "__main__":
    
    import utils
    
    args = utils.get_args(("-nfile",),
                          ("-min_occurs",))
    
    nfile = args.nfile
    min_occurs = args.min_occurs
    
    if nfile is None or min_occurs is None:
        e = "'nfile' and 'min_occurs' have no default values. Please, do:\n"
            ">>> python _filter_products.py -nfile <nfile> -min_occurs <min_occurs>"
        raise Exception(e)
    
    ninit, nfilt = filter_file(nfile=nfile, min_occurs=int(min_occurs))
    
    new_name = f"{nfile.split('.')[0]}_min_{min_occurs}.txt"
    print(f"Number of entries in {nfile!r}: {ninit}\n"
          f"Number of entries in {new_name!r}: {nfilt}")
