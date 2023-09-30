# -*- coding: utf-8 -*-

if __name__ == "__main__":
    
    """
    The idea is to build the model for "all" and the models for all other
    datasets in parallel. To do that, one has to open two terminal windows
    and write the following commands:
    ---
    [Terminal 1] % python _build_models.py -data all
    [Terminal 2] % python _build_models.py -data all-other
    ---
    The result are the final "init" models that will eventually be used to
    make predictions on "miscellaneous" proteins. According to these
    predictions the sequences are placed in the appropriate files and new
    models are trained based on the expanded database. These models will
    only be used if they yield better results than the original ones. Hence,
    the models trained here, may possibly be the ones comprising the tool.
    """
    
    import pathlib
    import utils
    from ml_model import MLModel
    
    # create "models" directory (if not already created)
    # ---
    pathlib.Path("../models").mkdir(exist_ok=True)
    
    # get command line argument "data"
    # ---
    args = utils.get_args(("-data",),
                          ("-init",))
    data, init = args.data, args.init
    
    # validate "-data" and "-init"
    # ---
    if data is None or init is None:
        raise ValueError("'data' and 'init' have no default values.")
    
    # validate "-data"
    # ---
    if data not in {"all", "all-other"}:
        raise ValueError("'data' must be in {'all', 'all-other'}.")
        
    # validate "-init"
    # ---
    if init not in {"yes", "no"}:
        raise ValueError("'init' must in {'yes', 'no'}.")
        
    # reassign "-init"
    # ---
    init = (init == "yes")
        
    # initialize "database" and "models"
    # ---
    if init:
        database, models = "../database", "../models"
    else:
        database, models = "../database_cs", "../models_cs"
        
    # mapping between datasets ("all-other") and ML algorithms
    # ---
    data_dict = {"dna-modification": "gradient-boosting",
                 "dna-replication": "gradient-boosting",
                 "lysis": "support-vector-machine",
                 "lysogeny-repressor": "support-vector-machine",
                 "packaging": "support-vector-machine",
                 "structural": "gradient-boosting"}
        
    # build model(s)
    # ---
    if data == "all":
        # build model for dataset "all"
        MLModel(data_path=f"{database}/all.csv",
                models_dir=models,
                algorithm="k-nearest-neighbors",
                test_size=0.2,
                init=init,
                final_model=True).build_model()
    else:
        # build models for all other datasets
        for dataset in data_dict:
            MLModel(data_path=f"{database}/{dataset.replace('-','_')}.csv",
                    models_dir=models,
                    algorithm=data_dict[dataset],
                    test_size=0.2,
                    init=init,
                    final_model=True).build_model()
    