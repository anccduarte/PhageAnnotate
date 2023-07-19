# -*- coding: utf-8 -*-

if __name__ == "__main__":
    
    import utils
    from pathlib import Path
    from predict_function import PredictFunction
    
    # create "results" directory (if not already created)
    Path("../results").mkdir(exist_ok=True)
    
    # read arguments (sys.argv)
    args = utils.get_args(("-path",),
                          ("-models_dir",))
    path = args.path
    models_dir = args.models_dir
    
    # check whether <path> and <models_dir> were provided to sys.argv
    if path is None or models_dir is None:
        e = ("'path' and 'models_dir' have no default values. Please, do:\n"
             ">>> python gene_products.py -path <path> -models_dir <models_dir>")
        raise Exception(e)
    
    # remaining arguments
    ttable = "11"
    icodons = ("TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG")
    
    # get name of the file
    name = path.split("/")[-1].split(".")[0]
    
    # get predictions
    PredictFunction(path=path,
                    models_dir=models_dir,
                    ttable=ttable,
                    icodons=icodons).save_results(name=name)
