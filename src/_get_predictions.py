# -*- coding: utf-8 -*-

if __name__ == "__main__":
    
    import argparse
    from pathlib import Path
    from predict_function import PredictFunction
    
    # initialize parser
    parser = argparse.ArgumentParser()
    
    # add argument
    parser.add_argument("-path")
    
    # read arguments
    args = parser.parse_args()
    path = args.path
    
    # check whether <path> was provided to sys.argv
    if path is None:
        raise Exception("'path' has no default value. Please, do:\n"
                        ">>> python gene_products.py -path <path>")
    
    # get name of the file
    name = path.split("/")[-1].split(".")[0]
    
    # get predictions
    Path("../results").mkdir(exist_ok=True)
    PredictFunction(path=path).save_results(name=name)
