# -*- coding: utf-8 -*-

import pandas as pd

# --------------------------------------------------------------------------------

# build dataframe of metrics from the accuracy, precision and recall scores
# obtained with the initial model and with the model fed on the extended datasets
# rows: models; columns: metrics (accuracy, precision, recall, f1-score)
# only modify values in "all_metrics"
# ---
def build_metrics_df() -> None:
    # ---
    names = ("all", "all-expanded",
             "dna-modification", "dna-modification-expanded",
             "dna-replication", "dna-replication-expanded",
             "lysis", "lysis-expanded",
             "lysogeny-repressor", "lysogeny-repressor-expanded",
             "packaging", "packaging-expanded",
             "structural", "structural-expanded")
    # ---
    all_metrics = ((87.44, 86.28, 82.50), (84.74, 82.96, 77.23),
                   (91.41, 90.81, 83.90), (79.45, 73.28, 66.06),
                   (91.60, 92.77, 85.80), (83.79, 84.04, 78.45),
                   (99.03, 98.62, 97.83), (96.04, 92.42, 90.54),
                   (91.94, 91.84, 84.25), (86.76, 88.52, 79.18),
                   (98.97, 98.69, 98.90), (97.33, 90.77, 82.88),
                   (83.60, 81.53, 78.54), (77.67, 74.24, 66.62))
    # ---
    metrics_dict = {"Accuracy": [], "Precision": [], "Recall": [], "F1-score": []}
    for metrics in all_metrics:
        accuracy, precision, recall = metrics
        f1_score = 2 * (precision*recall) / (precision+recall)
        metrics_dict["Accuracy"].append(accuracy)
        metrics_dict["Precision"].append(precision)
        metrics_dict["Recall"].append(recall)
        metrics_dict["F1-score"].append(round(f1_score, 2))
    # ---
    metrics_table = pd.DataFrame(metrics_dict, names)
    print(metrics_table)
    
# --------------------------------------------------------------------------------    
    
# rank gene products by their counts (more common gene products appear first in
# the new .txt file)
# ---
def rank_gene_products(tol: int) -> None:
    # ---
    def get_name(entry: str) -> str:
        first, *_ = entry.split(":")
        return first[2:]
    # ---
    def get_count(entry: str) -> int:
        *_, last = entry.split(":")
        return int(last[1:-1])
    # ---
    with open("../../gene_products/all.txt") as fin:
        entries = fin.readlines()
    # ---
    items = [(get_name(entry), c)
             for entry in entries if (c:=get_count(entry)) >= tol]
    sorted_items = sorted(items, key=lambda x: -x[1])
    # ---
    with open(f"../../gene_products/all_ranked_tol{tol}.txt", "w") as fout:
        for (name, count) in sorted_items:
            fout.write(f"- {name}: {count}\n")

