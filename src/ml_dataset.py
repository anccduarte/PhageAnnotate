# -*- coding: utf-8 -*-

import os
import pandas as pd
import warnings
from Bio import Entrez, SeqIO, Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from propy import PyPro

# Ignore Biopyhon warnings (related to translation)
warnings.filterwarnings("ignore")

# initialize Entrez.email (receive warnings in case of excessive usage of the E-utilities)
Entrez.email = "pg45464@alunos.uminho.pt"

class MLDataset:
    
    """
    Constructs a tabular dataset by computing a set of features/descriptors associated to
    the DNA sequences present in the input .fasta file.
    """
    
    def __init__(self, file: str, protein_name: str) -> None:
        """
        Initializes an instance of MLDataset.
        
        Parameters
        ----------
        file: str
            The name of the fasta file containing the DNA sequences
        protein_name: str
            The common name of the protein coded by the DNA sequences in <file>
            
        Attributes
        ----------
        _translation_table: str
            The identifier of the translation table used by the taxid in <file>
        """
        # parameters
        self.file = file
        self.protein_name = protein_name
        # attributes
        try:
            self._translation_table = self._get_translation_table()
        except:
            raise ValueError(f"The taxid identifier provided in {file!r} cannot be found in "
                             "NCBI's Taxonomy database.")
        
    def __repr__(self) -> str:
        """
        Returns the string representation of the object.
        """
        class_ = self.__class__.__name__
        return f"{class_}({self.file!r}, {self.protein_name!r})"
    
    def _get_translation_table(self) -> str:
        """
        Returns the identifier of the translation table to be used according to the 'taxid'
        present in the string <self.file>.
        """
        taxid = self.file.split("/")[-1].split("_")[0] # accounts for dir names
        with Entrez.efetch(db='taxonomy', id=taxid[4:], retmode='xml') as handle:
            record = Entrez.read(handle)
        return record[0]["GeneticCode"]["GCId"]
    
    def _translate_seq(self, seq: Seq.Seq) -> Seq.Seq:
        """
        Translates the DNA sequence it takes as input (according to <_translation_table>)
        and returns the corresponding sequence of aminoacids.
        
        Parameters
        ----------
        seq: Seq.Seq
            The DNA sequence
        """
        while len(seq) % 3 != 0: seq = f"N{seq}"
        tseq = Seq.Seq(seq).translate(table=self._translation_table)[:-1]
        return tseq
        
    @staticmethod
    def _nuc_composition(seq: Seq.Seq) -> list:
        """
        Returns a list containing the relative frequencies of each nucleotide in a DNA
        sequence.

        Parameters
        ----------
        seq: Seq.Seq
            The DNA sequence
        """
        comp = {f"{nuc}-Nuc": round(seq.count(nuc)/len(seq), 4) for nuc in "ACGT"}
        return comp

    @staticmethod
    def _amino_composition(seq: Seq.Seq) -> list:
        """
        Returns a list containing the relative frequencies of each aminoacid in a protein
        sequence.

        Parameters
        ----------
        seq: Seq.Seq
            The protein sequence
        """
        alphabet = "ACDEFGHIKLMNPQRSTVWY"
        comp = {f"{a}-Amino": round(seq.count(a)/len(seq), 4) for a in alphabet}
        return comp
    
    @staticmethod
    def _dipeptide_composition(seq: Seq.Seq) -> list:
        """
        Returns the dipeptide composition of the protein sequence.
        
        Parameters
        ----------
        seq: Seq.Seq
            The protein sequence
        """
        DesObject = PyPro.GetProDes(seq)
        return DesObject.GetDPComp()
    
    @staticmethod
    def _aromaticity(seq: Seq.Seq) -> list:
        """
        Returns a list containing the value of the protein's aromaticity (sum of the
        relative frequencies of the aminoacids phenylalanine (F), tyrosine (Y) and
        tryptophan (W)).
        
        Parameters
        ----------
        seq: Seq.Seq
            The protein sequence
        """
        arom = sum([(seq.count(amino)/len(seq)) for amino in "FYW"])
        return {"Aromaticity": round(arom, 4)}
    
    @staticmethod
    def _isoelectric_point(seq: Seq.Seq) -> list:
        """
        Returns a list containing the value of the protein's isoelectric point.
        
        Parameters
        ----------
        seq: Seq.Seq
            The protein sequence
        """
        X = ProteinAnalysis(str(seq))
        return {"Isoelectric": round(X.isoelectric_point(), 4)}
    
    @staticmethod
    def _secondary_structure_fraction(seq: Seq.Seq) -> list:
        """
        Returns a list containing the fraction of amino acids which tend to be in the
        secondary structures:
        - alpha-helix: Val (V), Ile (I), Tyr (Y), Phe (F), Trp (W), Leu (L)
        - beta-turn: Asn (N), Pro (P), Gly (G), Ser (S)
        - beta-sheet: Glu (E), Met (M), Ala (A), Leu (L)
        
        Parameters
        ----------
        seq: Seq.Seq
            The protein sequence
        """
        helix = round(sum([(seq.count(amino)/len(seq)) for amino in "VIYFWL"]), 4)
        turn = round(sum([(seq.count(amino)/len(seq)) for amino in "NPGS"]), 4)
        sheet = round(sum([(seq.count(amino)/len(seq)) for amino in "EMAL"]), 4)
        ssf = {"A-Helix": helix, "B-turn": turn, "B-sheet": sheet}
        return ssf
    
    @staticmethod
    def _ctd_descriptors(seq: Seq.Seq) -> list:
        """
        Returns a list containing the CTD descriptors (composition, transition, distribution)
        of the protein sequence.
        
        Parameters
        ----------
        seq: Seq.Seq
            The protein sequence
        """
        DesObject = PyPro.GetProDes(seq)
        return DesObject.GetCTD()

    def build_dataset(self) -> pd.DataFrame:
        """
        Returns a tabular dataset (represented by a pd.DataFrame) containing the featurized
        sequences (present in <file>).
        ---
        Feature mapping:
        - sequence length (1 feature), nucleotide composition (4 features), dipeptide
        composition (400 features), aminoacid composition (20 features), aromaticity (1
        feature), isoelectric point (1 feature), secondary structure fraction (3 features),
        ctd descriptors (147 features)
        ---
        Total number of features: 577
        """
        # initialize dataset
        dataset = []
        # read sequences present in <file>
        sequences = SeqIO.parse(self.file, format="fasta")
        # iterate through sequences and build dataset
        for i, seq in enumerate(sequences):
            # translate DNA sequence
            tseq = self._translate_seq(seq.seq)
            # get sequence descriptors
            len_ = {"Len-Prot": len(tseq)}
            nuc = MLDataset._nuc_composition(seq.seq)
            amino = MLDataset._amino_composition(tseq)
            dipep = MLDataset._dipeptide_composition(tseq)
            arom = MLDataset._aromaticity(tseq)
            isoel = MLDataset._isoelectric_point(tseq)
            ssf = MLDataset._secondary_structure_fraction(tseq)
            ctd = MLDataset._ctd_descriptors(tseq)
            # build entry and append it to "features"
            entry = {**len_, **nuc, **amino, **dipep, **arom, **isoel, **ssf, **ctd}
            dataset.append(entry)
        # build pd.DataFrame from "dataset"
        df = pd.DataFrame(dataset)
        # add label column to the dataframe
        label_col = [self.protein_name] * (i+1)
        df["Function"] = label_col
        # return df containing the featurized records
        return df
        