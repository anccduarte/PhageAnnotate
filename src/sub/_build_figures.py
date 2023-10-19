# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np

def _construct_figure(metrics: dict,
                      x_labels: list,
                      name_figure: str) -> None:
    """
    Constructs a figure regarding a particular evaluation metric. It is
    exploited by "hierarchical" and "hierarchical_x".
    
    Parameters
    ----------
    metrics: dict
        A dictionary of evaluation metrics
    x_labels: list
        A list of two items containing the x-axis label "ticks"
    name_figure: str
        The name to be given to the constructed figure
    """
    # initialize labels
    labels = ("ALL: GB", "DNA-mod: GB", "DNA-rep: GB", "Lysis: SVM",
              "Lys-rep: SVM", "Pack: SVM", "Struct: GB")
    # construct and assemble figures
    fig, ax = plt.subplots(2, 2)
    coords = ((0, 0), (0, 1), (1, 0), (1, 1))
    for coord, (name, vals) in zip(coords, metrics.items()):
        for i, f, l in zip(vals["init"], vals["final"], labels):
            ax[coord].plot(x_labels, [i, f], '.-', label=l)
        temp = np.linspace(min(vals["init"]+vals["final"])-2, 100, 5)
        ax[coord].set_yticks([int(val) for val in temp])
        ax[coord].set_title(name)
    # add legend and save assembly of figures
    plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    plt.savefig(f"figures/{name_figure}.png", bbox_inches="tight")
    plt.close()
    
def construct_figure(metrics: dict,
                     x_labels: list,
                     name_figure: str) -> None:
    """
    Constructs a figure regarding a particular evaluation metric. It is
    exploited by "hierarchical" and "hierarchical_x".
    
    Parameters
    ----------
    metrics: dict
        A dictionary of evaluation metrics
    x_labels: list
        A list of two items containing the x-axis labels
    name_figure: str
        The name to be given to the constructed figure
    """
    # ---
    labels = ("ALL", "DNA-mod", "DNA-rep", "Lysis", "Lys-rep",
              "Pack", "Struct")
    # ---
    fig, ax = plt.subplots(2, 2)
    fig.tight_layout(h_pad=2); fig.set_size_inches(10, 6)
    # ---
    coords = ((0, 0), (0, 1), (1, 0), (1, 1))
    for coord, (name, vals) in zip(coords, metrics.items()):
        # ---
        x = np.arange(len(labels))
        width = 0.25
        multiplier = 0
        # ---
        for x_label, (_, values) in zip(x_labels, vals.items()):
            offset = width * multiplier
            ax[coord].bar(x+offset,
                          values,
                          width,
                          align="center",
                          label=x_label)
            multiplier += 1
        # ---
        ax[coord].set_ylabel("%")
        ax[coord].set_title(name)
        ax[coord].set_xticks(x+width,
                             labels,
                             fontsize=8,
                             rotation=45)
        ax[coord].set_ylim(70, 100)
    # ---
    plt.legend(bbox_to_anchor=(0.75, -0.25))
    plt.savefig(f"figures/{name_figure}.png", bbox_inches="tight")
    plt.close()

def hierarchical() -> None:
    """
    Builds figures comparing the values attained for each evaluation
    metric regarding the initial and final Hierarchical models. These
    metrics are the accuracy score, precision score, recall score and
    f1-score.
    """
    # computed scores (accuracy, precision, recall and f1)
    a = {"init":  [88.56, 88.53, 87.91, 98.64, 89.90, 98.60, 87.44],
         "final": [90.24, 89.20, 89.21, 98.67, 90.18, 98.45, 87.63]}
    p = {"init":  [86.90, 87.36, 89.80, 97.56, 85.43, 98.27, 88.99],
         "final": [89.22, 89.05, 90.10, 97.24, 85.06, 98.12, 88.65]}
    r = {"init":  [84.02, 77.63, 77.39, 96.79, 82.20, 98.35, 81.90],
         "final": [86.14, 76.97, 78.85, 97.36, 82.59, 98.15, 81.80]}
    f = {"init":  [85.35, 81.95, 82.23, 97.17, 83.66, 98.31, 84.92],
         "final": [87.57, 82.17, 83.38, 97.30, 83.72, 98.14, 84.77]}
    # construct figures for each evaluation metric
    metrics = {"accuracy": a, "precision": p, "recall": r, "f1-score": f}
    x_labels = ["Hierarchical-init", "Hierarchical-final"]
    construct_figure(metrics, x_labels, "hierarchical")
    
def hierarchical_x() -> None:
    """
    Builds figures comparing the values attained for each evaluation
    metric regarding Hierarchical and Hierarchical-X. These metrics are
    the accuracy score, precision score, recall score and f1-score.
    """
    # computed scores (accuracy, precision, recall and f1)
    a = {"h":   [90.24, 89.20, 89.21, 98.67, 90.18, 98.45, 87.63],
         "h-x": []}
    p = {"h":   [89.22, 89.05, 90.10, 97.24, 85.06, 98.12, 88.65],
         "h-x": []}
    r = {"h":   [86.14, 76.97, 78.85, 97.36, 82.59, 98.15, 81.80],
         "h-x": []}
    f = {"h":   [87.57, 82.17, 83.38, 97.30, 83.72, 98.14, 84.77],
         "h-x": []}
    # construct figures for each evaluation metric
    metrics = {"accuracy": a, "precision": p, "recall": r, "f1-score": f}
    x_labels = ["Hierarchical", "Hierarchical-X"]
    construct_figure(metrics, x_labels, "hierarchical_x")

def miscellaneous() -> None:
    """
    Constructs a plot depicting the evolution of the error and the number
    of "informative" predictions attained via predicting the function of
    miscellaneous sequences using various sets of thresholds.
    """
    # thresholds and associated "informative" predictions and errors
    thresh = [f"0.{val}" for val in range(50, 100, 5)]
    inform = [475630, 411700, 355110, 304890, 258450, 216190, 177110,
              139420, 102040, 62470]
    errors = [62, 73, 58, 47, 51, 53, 51, 46, 40, 32]
    # construct plot thresholds-errors
    fig, ax = plt.subplots(1, 2); fig.set_size_inches(12, 4)
    left, right = "# 'informative' predictions", "Error"
    ax[0].plot(thresh, inform, '.-r'); ax[0].set_title(left)
    ax[1].plot(thresh, errors, '.-b'); ax[1].set_title(right)
    # save assembly of figures
    plt.savefig(f"figures/misc_thresholds.png")
    plt.close()
        
    
if __name__ == "__main__":
    
    import sys
    sys.path.append("..")
    import utils
    from pathlib import Path
    
    Path("figures").mkdir(exist_ok=True)
    
    args = utils.get_args(("-action",))
    action = args.action
    
    if action is None:
        raise ValueError("'action' has no default value")
        
    options = {"hierarchical": hierarchical,
               "hierarchical_x": hierarchical_x,
               "miscellaneous": miscellaneous}
    
    if action not in options:
        opt_set = "{" + ", ".join(list(options.keys())) + "}"
        raise ValueError(f"'action' must be in {opt_set}")
    
    func = options[action]
    func()
