# -*- coding: utf-8 -*-

if __name__ == "__main__":
    
    import os
    import pandas as pd
    from ml_dataset import MLDataset
    from pathlib import Path

    # initialize database (will contain the datasets)
    database = "../database"
    Path(database).mkdir(exist_ok=True)

    # directories where the DNA sequences are stored
    directories = {"structural": "structural",
                   "lysis": "lysis",
                   "modification_replication": "modification-replication",
                   "packaging": "packaging"}

    # initialize df for all functional classes
    df_all = pd.DataFrame()

    # initialize column to add to "df_all"
    to_add = []

    # main loop -> iterate through "directories"
    for dir_ in directories:
        # get functional class
        func_class = directories[dir_]
        # log "state" of the program to console
        print(f"Building dataset of {func_class!r} proteins...", end=" ")
        # change current working directory to "dir_"
        os.chdir(f"../sequences/{dir_}")
        # initialize df for specific functional class
        df_class = pd.DataFrame()
        # counter for number of proteins in cwd's files
        counter = 0
        # iterate through files in "dir_"
        for file in os.scandir():
            # get name of the file
            n_file = file.name
            # ignore files which do not contain DNA sequences 
            if not n_file.endswith("fasta"):
                continue
            # get name of the protein coded by the DNA sequences in "file"
            _, *temp, _ = n_file.split("_")
            prot_name = "-".join(temp)
            # featurize DNA sequences in "file"
            df = MLDataset(n_file, prot_name).build_dataset()
            # update sequence couter for functional class
            counter += df.shape[0]
            # concat "df" to "df_class" and "df_all"
            df_class = pd.concat([df_class, df])
            df_all = pd.concat([df_all, df])
        # go back to the parent directory ("sequences")
        os.chdir("..")
        # save "df_class" in the directory "database"
        df_class.to_csv(f"{database}/{dir_}.csv")
        # update "to_add" with the name of the functional class
        to_add += [func_class] * counter
        # log to console when construction of the dataset is finished
        print("DONE")

    # delete column with protein names from "df_all"
    del df_all["Function"]

    # add new column to "df_all" (functional classes)
    df_all["Func-Class"] = to_add

    # save "df_all" in the directory "database"
    df_all.to_csv(f"{database}/all.csv")
    