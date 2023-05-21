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
    directories = {"dna_modification": "dna-modification",
                   "dna_replication": "dna-replication",
                   "lysis": "lysis",
                   "lysogeny_repressor": "lysogeny-repressor",
                   "packaging": "packaging",
                   "structural": "structural",
                   "other": "other"}

    # initialize df for all functional classes
    df_all = pd.DataFrame()

    # initialize column to add to "df_all"
    to_add = []

    # main loop -> iterate through "directories"
    for i, dir_ in enumerate(directories):
        # get functional class
        func_class = directories[dir_]
        print(f"Building dataset of {func_class!r} proteins...")
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
            _, *temp = n_file.split("_")
            temp[-1] = temp[-1].split(".")[0]
            prot_name = "-".join(temp)
            # featurize DNA sequences in "file"
            df = MLDataset(n_file, prot_name).build_dataset()
            # update sequence couter for functional class
            counter += df.shape[0]
            # concat "df" to "df_class" and "df_all"
            df_class = pd.concat([df_class, df])
            df_all = pd.concat([df_all, df])
            print(f"- Added {prot_name!r} sequences")
        # go back to the parent directory ("sequences")
        os.chdir("..")
        # save "df_class" in the directory "database"
        df_class.to_csv(f"{database}/{dir_}.csv")
        # update "to_add" with the name of the functional class
        to_add += [func_class] * counter
        # separate logs of different functional classes
        if i < len(directories)-1:
            print("---")

    # delete column with protein names from "df_all"
    del df_all["Function"]

    # add new column to "df_all" (functional classes)
    df_all["Func-Class"] = to_add

    # save "df_all" in the directory "database"
    df_all.to_csv(f"{database}/all.csv")
    