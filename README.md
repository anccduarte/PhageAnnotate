# PhageAnnotate

## Background

The advent of bacterial strains resistant to virtually all currently available antibiotic agents is indubitably alarming. Hence, alternative strategies for efficaciously combating such multi-drug resistant bacterial strains ought to be considered. Bacteriophages (phages), the most abundant biological entities on Earth, may play a key role in such a therapeutic revolution. Phage therapy exploits the evolutionary context of phages and bacteria: millions of years of host-parasite coevolution made phages lethal bacterial predators. <br>

Nonetheless, the process of “recruiting” phages that competently combat bacterial infections is profoundly dependent on the thorough understanding of their biology and safety. Yet, under the tremendous diversity of phage’s genomes, most of their genes cannot be assigned to functions via homology-based techniques, significantly hindering such fundamental insights. <br>

## *Modus operandi*

`PhageAnnotate`, a phage annotation system powered by machine learning (*ML*), constitutes an attempt to transcend the aforementioned obstacle. The system is comprised of seven *ML* models: one model discerning functional classes (i.e., *DNA-modification*, *DNA-replication*, *lysis*, *lysogeny-repressor*, *packaging* and *structural*), and six models distinguishing functional roles pertaining a specific functional class. `PhageAnnotate`'s course of action reflects such structuring: given a new *DNA* sequence, its functional class is firstly predicted, and, only then, depending on the latter prediction, the sequence is dispatched to the appropriate model discerning functional roles. Such a hierarchical architecture allows, in a first instance, to crudely classify sequences in functional classes before performing a more fine-grained functional role categorization. The ascertained functional classes and the corresponding functional roles are depicted below. <br>


| *DNA-modification* | *DNA-replication* | *lysis* | *lysogeny-repressor* | *packaging* | *structural*
|:---:|:---:|:---:|:---:|:---:|:---:|
|*nuclease* <br> *deaminase* <br> *thymidylate synthase* <br> *dUTPase* <br> *kinase* <br> *phosphoesterase* <br> *reductase* <br> *methylase* <br> *ATPase* <br> *nucleotidohydrolase* <br> *transferase* <br> *phosphohydrolase* <br> *glutaredoxin* <br> *ribonuclease* | *transcription factor* <br> *DNA primase-helicase* <br> *RNA polymerase* <br> *DNA ligase* <br> *DNA polymerase* <br> *RNA ligase* <br> *replication initiation* | *endolysin* <br> *holin* <br> *spanin* | *integrase* <br> *recombinase* <br> *repressor* <br> *resolvase* <br> *transposase* <br> *antirepressor* | *large terminase* <br> *small terminase* | *minor tail* <br> *major tail* <br> *portal* <br> *minor capsid* <br> *major capsid* <br> *head-tail* <br> *tail fiber* <br> *tail sheath* <br> *baseplate* <br> *neck* <br> *collar* <br> *tailspike* |


`PhageAnnotate`'s predictions are guided by a set of two probability thresholds. These thresholds are a proxy for the degree of confidence in the model's predictions, and must be contained in the interval $[0.00, 1.00]$. <br>

1. The first threshold defaults to $0.65$ and determines the degree of confidence in functional class predictions. For example, if this threshold is set to $0.80$, and the model predicts a given *DNA* sequence to belong to the functional class *DNA-modification* with $0.70$ confidence, the prediction is disregarded and the sequence is labeled with *"other"*. In such cases, the prediction endeavor ends at this first classification step, as attributing a functional role to a *DNA* sequence whose functional class is *"other"* is nonsensical.

2. The second probability threshold (also defaulting to $0.65$) specifies the degree of confidence on functional role predictions, and, as hinted in the previous point, it is only employed whenever the functional class categorization is distinct from *"other"*. For example, assume that a given *DNA* sequence was assigned to the functional class *DNA-modification*, and this threshold is set to $0.80$. Now, assume that the model discerning functional roles pertaining to the functional class *DNA-modification* categorizes the sequence as *nuclease*. If the degree of confidence in such prediction is greater than o equal to $0.80$, the sequence is bounded to *nuclease*; otherwise, it is categorized as *"other-function"*.

***Currently, `PhageAnnotate` solely provides a command-line interface. However, it does not require extensive experience with such tooling. To install and use `PhageAnnotate` a handful of commands must be executed. These commands are succintely described below.***


## Installation

1. Install `miniconda`. The software is easy to install, only involving the execution of a few straightforward steps. The [miniconda's installation guide](https://docs.conda.io/projects/miniconda/en/latest/index.html#quick-command-line-install) clearly details such steps.

2. Download `PhageAnnotate`'s repository to the directory of your choice. <br>
`cd <directory>` <br>
`git clone https://github.com/anccduarte/PhageAnnotate.git`

3. Access `PhageAnnotate`'s directory. <br>
`cd PhageAnnotate`

4. Create a `conda` environment and install `PhageAnnotate`'s dependencies. <br>
`conda env create -f environment.yml`


## Usage

1. If not already in `PhageAnnotate`'s directory, access it. <br>
`cd <path to PhageAnnotate>`

2. Create a folder `genomes`. <br>
`mkdir genomes`

3. Place the *FASTA* file(s) of phage *DNA* sequences in `genomes`.

4. Activate the `conda` environment created upon `PhageAnnotate`'s installation. <br>
`conda activate phageannotate_env`

5. Run the tool for some *FASTA* file stored in `genomes`. In the first example, the default set of thresholds (i.e., $[0.65, 0.65]$) is exploited to guide the system's predictions. In the second and third examples, manually specified sets of probability thresholds are utilized by the tool (i.e., $[0.80, 0.65]$ and $[0.85, 0.75]$, respectively). <br>
`python predict.py -file <FASTA file>` <br>
`python predict.py -file <FASTA file> -thresh_class 0.80` <br>
`python predict.py -file <FASTA file> -thresh_class 0.85 -thresh_role 0.75`

6. Check the results in `results`.
