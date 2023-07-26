# -*- coding: utf-8 -*-

if __name__ == "__main__":
    
    # command line arguments
    # ---
    # - final_model: -final_model (options: "yes" or "no")
    # - init: -init (options: "yes" or "no")
    
    # deduced from command line arguments
    # ---
    # - models_dir: deduced from <init> ("../models" or "../models_cs")
    # - db_name: deduced from <init> ("../database" or "../database_cs")
    # - test_size: deduced from <final_model> (0.1 or 0.2)
    # - behavior of the loop datasets/algorithms determined by <final_model>
    
    import utils
    from ml_model import MLModel
    from pathlib import Path
    
    # get command line arguments
    # ---
    args = utils.get_args(("-final_model",),
                          ("-init",))
    final_model = args.final_model
    init = args.init
    
    # check validity of command line arguments (1)
    # ---
    if final_model is None or init is None:
        raise Exception("'final_model' and 'init' have no default values. Please do:\n"
                        ">>> python -final_model <final_model> -init <init>")
        
    # check validity of <final_model> (2)
    # ---
    if final_model not in {"yes", "no"}:
        raise ValueError(f"{final_model!r} is not valid for 'final_model'."
                         "Choose one of {{'yes', 'no'}}.")
        
    # check validity of <init> (3)
    # ---
    if init not in {"yes", "no"}:
        raise ValueError(f"{init!r} is not valid for 'init'."
                         "Choose one of {{'yes', 'no'}}.")
        
    # reassign <final_model> and <init>
    # ---
    final_model = True if final_model == "yes" else False
    init = True if init == "yes" else False
    
    # initialize <models_dir>
    # ---
    models_dir = "../models" if init else "../models_cs"
    Path(models_dir).mkdir(exist_ok=True)
    
    # initialize <db_name> and <test_size>
    # ---
    db_name = "database" if init else "database_cs"
    test_size = 0.1 if final_model else 0.2
        
    # initialize <datasets>
    # ---
    datasets = ["all", "dna_modification", "dna_replication", "lysis",
                "lysogeny_repressor", "packaging", "structural"]
    
    # initilaize <algorithms>
    # ---
    algorithms = ["gradient-boosting", "random-forest", "naive-bayes",
                  "k-nearest-neighbors", "neural-network", "support-vector-machine",
                  "decision-tree"]
    
    # logic for model construction
    # ---
    # - if <final_model> is set to True, the combinations dataset-algorithm are
    #   defined -> ask for list of algorithms that match each dataset and zip the
    #   datasets and the algorithms
    # - if <final_model> is set to False, initialize list with all algorithms and
    #   test all combinations dataset-algorithm
    # ---
    if final_model:
        # get list of algorithms
        print(f"\nType the list of ML algorithms matching the datasets:\n{datasets}")
        algorithms_in = eval(input("Algorithms: "))
        # check validity of the input
        if any(algo not in algorithms for algo in algorithms_in):
            raise ValueError("The list of algorithms is not valid. Please choose "
                             f"one of the following:\n{{{', '.join(algorithms)}}}")
        # zip datasets and algorithms and train final models
        for data, algo in zip(datasets, algorithms_in):
            MLModel(data_path=f"../{db_name}/{data}.csv",
                    models_dir=models_dir,
                    algorithm=algo,
                    test_size=test_size,
                    init=init,
                    final_model=final_model).build_model()
            print("------------")
    else:
        # test all combinations dataset-algorithm
        for data in datasets:
            for algo in algorithms:
                MLModel(data_path=f"../{db_name}/{data}.csv",
                        models_dir=models_dir,
                        algorithm=algo,
                        test_size=test_size,
                        init=init,
                        final_model=final_model).build_model()
                print("------------")
