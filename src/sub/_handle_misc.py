# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import os
import sys
sys.path.append("..")
import utils
from Bio import SeqIO
from pathlib import Path

def filter_miscellaneous() -> None:
    """
    Filters miscellaneous sequences to only include entries whose
    description does not include ambiguous terms. It may be pertinent
    to include other terms; nonetheless, the more problematic ones
    (i.e., occur more reccurrently) are "hypothetical", "unknown",
    "no-product" and "putative".
    """
    # ---
    def _valid(description: str) -> bool:
        descrip = description.lower()
        return "hypothetical" not in descrip and \
               "unknown" not in descrip and \
               "no-product" not in descrip and \
               "putative" not in descrip
    # ---
    orig, filt = 0, 0
    with open(f"{BASE}/miscellaneous_filt.fasta", "w") as fout:
        nfile = f"{BASE}/txid2731619_miscellaneous.fasta"
        records = SeqIO.parse(nfile, "fasta")
        for rec in records:
            orig += 1
            if not _valid(rec.description):
                continue
            filt += 1
            fout.write(f">{rec.description}\n{rec.seq}\n\n")
    # ---
    print(f"# sequences before filtering: {orig}\n"
          f"# sequences after filtering: {filt}")
    
def sample_miscellaneous() -> None:
    """
    Samples 10% of the sequences from the collection of filtered
    miscellaneous sequences. Such step is imperative given the huge
    amount of time needed to make predictions on the whole set of
    filtered miscellaneous sequences.
    """
    # ---
    def count_seqs(nfile: str) -> int:
        count = 0
        records = SeqIO.parse(nfile, "fasta")
        for _ in records:
            count += 1
        return count
    # ---
    seqs_filt = f"{BASE}/miscellaneous_filt.fasta"
    num_seqs_filt = count_seqs(seqs_filt)
    sample_size = num_seqs_filt // 10
    with open(f"{BASE}/misc_filt_sample.fasta", "w") as fout:
        records = SeqIO.parse(seqs_filt, "fasta")
        samples = np.random.permutation(num_seqs_filt)[:sample_size]
        for i, rec in enumerate(records):
            if i not in samples:
                continue
            fout.write(f">{rec.description}\n{rec.seq}\n\n")
            
def predict_misc_filt_sample() -> None:
    """
    Makes predictions on the set of filtered miscellaneous sequences.
    Encompasses an iterative process in which predictions are made on
    the same collection of sequences, but the set of thresholds is
    variable.
    """
    # establish tresholds
    thresholds = ("[0.50, 0.50]", "[0.55, 0.55]", "[0.60, 0.60]",
                  "[0.65, 0.65]", "[0.70, 0.70]", "[0.75, 0.75]",
                  "[0.80, 0.80]", "[0.85, 0.85]", "[0.90, 0.90]")
    # predict for all thresholds
    with utils.change_dir(".."):
        for thresh in thresholds:
            # log <thresh>
            print(f"{thresh = }")
            # execute predictions for <thresh>
            cmd = "python _predict.py " + \
                  "-path ../_miscellaneous/misc_filt_sample.fasta " + \
                  "-models_dir ../models " + \
                  f"-thresholds {thresh!r}"
            os.system(cmd)
            # rename .csv file (containing predictions)
            t = thresh[3:5]
            old_name = "../results/misc_filt_sample.csv"
            new_name = f"../results/misc_filt_sample_{t}.csv"
            os.system(f"mv {old_name} {new_name}")
    
def sample_misc_preds(nfile: str, rand: int) -> None:
    """
    Samples a total of 100 informative predictions from the .csv file
    pointed to by <nfile>. A prediction is regarded as informative if
    two conditions are met: 1. "Functional-Class" does not correspond
    to the string "other"; 2. "Function" is not "other-function".
    
    Parameters
    ----------
    nfile: str
        The name of the .csv file containing the predictions
    rand: int
        Controls the sampling of "informative" predictions
    """
    # ---
    def get_true_names(descrip_col: pd.Series) -> list:
        names = []
        for descrip in descrip_col.tolist():
            names.append(descrip.split("|")[-1])
        return names
    # ---
    # create directory for storing the samples
    Path(f"{BASE}/samples").mkdir(exist_ok=True)
    # read data to memory
    preds = pd.read_csv(nfile, index_col=False)
    # remove rows whose "Function" is either "" or "other-function"
    informative = preds.loc[(preds["Functional-Class"] != "other") &
                            (preds["Function"] != "other-function")]
    # sample a total of 100 informative predictions
    n_info = informative.shape[0]
    samples = np.random.RandomState(rand).permutation(n_info)[:100]
    informative_sample = informative.iloc[samples]
    # get true names associated to the samples and replace Description
    true_names = get_true_names(informative_sample.Description)
    out_df = informative_sample.iloc[:, 2:].reset_index(drop=True)
    out_df["True-Label"] = true_names
    # save "informative" samples to .csv file
    fname = nfile.split("/")[-1]
    out_df.to_csv(f"{BASE}/samples/{fname}", index=False)
    # log useful information
    print(f"# of informative predictions: {n_info}\n"
          f"# of sampled informative predictions: 100")
            
    
if __name__ == "__main__":
    
    import sys
    sys.path.append("..")
    import utils
    
    BASE = "../../_miscellaneous"
    Path(BASE).mkdir(exist_ok=True)
    
    args = utils.get_args(("-action",))
    action = args.action
    
    if action is None:
        raise ValueError("'action' has no default value")
        
    options = {"filter": filter_miscellaneous,
               "sample": sample_miscellaneous,
               "predict": predict_misc_filt_sample,
               "sample_preds": sample_misc_preds}
    
    if action not in options:
        opts = ", ".join(list(options.keys()))
        raise ValueError(f"'action' must be in {{{opts}}}")
    
    func = options[action]
    if action == "sample_preds":
        func(nfile=f"../../results/{input('nfile: ')}", rand=0)
    else:
        func()
    