# -*- coding: utf-8 -*-

if __name__ == "__main__":
    
    import sys
    sys.path.append("src")
    import utils
    from pathlib import Path
    from predict_function import PredictFunction
    
    # create "results" directory (if not already created)
    Path("results").mkdir(exist_ok=True)
    
    # read arguments (sys.argv)
    args = utils.get_args(("-file",),
                          ("-thresh_class", "0.65"),
                          ("-thresh_role", "0.65"))
    file = args.file
    try:
        thresh_class = float(args.thresh_class)
        thresh_role = float(args.thresh_role)
    except:
        tc, tr = "thresh_class", "thresh_role"
        raise ValueError(f"{tc!r} and {tr!r} must be numerals")
    
    # check whether <file> was provided to sys.argv
    if file is None:
        raise ValueError(("'file' has no default value. Please, do:\n"
                          "% python predict.py -file <path_to_file>"))
    
    # check validity of thresholds
    thresholds = [thresh_class, thresh_role]
    if any((t<0 or t>1) for t in thresholds):
        raise ValueError("Each provided threshold must be in [0, 1].")
        
    # remaining arguments
    ttable = "11"
    icodons = ("TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG")
    
    # get name of the file
    name = file.split("/")[-1].split(".")[0]
    
    # get predictions
    PredictFunction(path=file,
                    models_dir="models",
                    thresholds=thresholds,
                    ttable=ttable,
                    icodons=icodons).predict(name=name)
    