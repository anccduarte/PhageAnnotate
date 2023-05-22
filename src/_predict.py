# -*- coding: utf-8 -*-

if __name__ == "__main__":
    
    import argparse
    from pathlib import Path
    from predict_function import PredictFunction
    
    # create "results" directory (if not already created)
    Path("../results").mkdir(exist_ok=True)
    
    # initialize parser
    parser = argparse.ArgumentParser()
    
    # add argument
    parser.add_argument("-path")
    
    # read arguments
    args = parser.parse_args()
    path = args.path
    
    # remaining arguments
    ttable = "11"
    icodons = ("TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG")
    
    # check whether <path> was provided to sys.argv
    if path is None:
        raise Exception("'path' has no default value. Please, do:\n"
                        ">>> python gene_products.py -path <path>")
    
    # get name of the file
    name = path.split("/")[-1].split(".")[0]
    
    # get predictions
    PredictFunction(path=path,
                    ttable=ttable,
                    icodons=icodons).save_results(name=name)
