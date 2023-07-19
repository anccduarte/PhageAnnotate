# -*- coding: utf-8 -*-

import os
import pandas as pd
import time
import warnings
from Bio import Entrez, SeqIO, Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from propy import PyPro

# Ignore Biopyhon warnings (related to sequence translation)
warnings.filterwarnings("ignore")

class MLDataset:
    
    """
    Constructs a tabular dataset by computing a set of features/descriptors associated to
    the DNA sequences present in the input .fasta file.
    """
    
    def __init__(self, file: str, prot_name: str, ttable: str, icodons: tuple) -> None:
        """
        Initializes an instance of MLDataset.
        
        Parameters
        ----------
        file: str
            The name of the .fasta file containing the DNA sequences
        prot_name: str
            The common name of the protein coded by the DNA sequences in <file>
        ttable: str
            The identifier of the translation table to be used
        icodons: tuple
            A tuple of possible initiation codons given the table <ttable>
        """
        # parameters
        self.file = file
        self.prot_name = prot_name
        self.ttable = ttable
        self.icodons = icodons
        
    def __repr__(self) -> str:
        """
        Returns the string representation of the object.
        """
        c = self.__class__.__name__
        r = f"{c}({self.file!r}, {self.prot_name!r}, {self.ttable!r}, {self.icodons!r})"
        return r
    
    def _translate_seq(self, seq: Seq.Seq) -> Seq.Seq:
        """
        Translates the DNA sequence it takes as input (according to <self.ttable>) and
        returns the corresponding sequence of aminoacids.
        
        Parameters
        ----------
        seq: Seq.Seq
            The DNA sequence
        """
        valid = lambda s: len(s) % 3 == 0
        if len(seq) % 3 != 0:
            icodon = seq[:3]
            if icodon in self.icodons:
                while not valid(seq): seq = f"{seq}N"
            else:
                while not valid(seq): seq = f"N{seq}"
        tseq = Seq.Seq(seq).translate(table=self.ttable)[:-1]
        return tseq
    
    @staticmethod
    def _get_valid_protein(seq: Seq.Seq) -> Seq.Seq:
        """
        Returns a valid protein sequence from the sequence it takes as input (substitutes
        invalid characters for the empty string).
        
        Parameters
        ----------
        seq: Seq.Seq
            The protein sequence to be validated
        """
        invalid_chars = "*BJOUWXZ"
        for c in invalid_chars:
            seq = seq.replace(c, "")
        return seq
        
    @staticmethod
    def _nuc_composition(seq: Seq.Seq) -> dict:
        """
        Returns a dictionary containing the relative frequencies of each nucleotide in
        a DNA sequence.

        Parameters
        ----------
        seq: Seq.Seq
            The DNA sequence
        """
        return {f"{nuc}-Nuc": seq.count(nuc)/len(seq) for nuc in "ACGT"}

    @staticmethod
    def _amino_composition(seq: Seq.Seq) -> dict:
        """
        Returns a dictionary containing the relative frequencies of each aminoacid in
        a protein sequence.

        Parameters
        ----------
        seq: Seq.Seq
            The protein sequence
        """
        amino_comp = PyPro.GetProDes(seq).GetAAComp()
        return {f"{a}-Amino": val for (a, val) in amino_comp.items()}
    
    @staticmethod
    def _dipeptide_composition(seq: Seq.Seq) -> dict:
        """
        Returns a dictionary containing the dipeptide composition of the protein sequence
        it takes as input.
        
        Parameters
        ----------
        seq: Seq.Seq
            The protein sequence
        """
        dipept_comp = PyPro.GetProDes(seq).GetDPComp()
        return {f"{di}-Dipept": val for (di, val) in dipept_comp.items()}
    
    @staticmethod
    def _molecular_weight(seq: Seq.Seq) -> dict:
        """
        Returns a dictionary containing the value of the protein's molecular weight.
        
        Parameters
        ----------
        seq: Seq.Seq
            The protein sequence
        """
        return {"Molecular-Weight": ProteinAnalysis(seq).molecular_weight()}
    
    @staticmethod
    def _aromaticity(seq: Seq.Seq) -> dict:
        """
        Returns a dictionary containing the value of the protein's aromaticity (sum of
        the relative frequencies of the aminoacids phenylalanine (F), tyrosine (Y) and
        tryptophan (W)).
        
        Parameters
        ----------
        seq: Seq.Seq
            The protein sequence
        """
        return {"Aromaticity": ProteinAnalysis(seq).aromaticity()}
    
    @staticmethod
    def _instability_index(seq: Seq.Seq) -> dict:
        """
        Returns a dictionary containing the value of the protein's instability index.
        
        Parameters
        ----------
        seq: Seq.Seq
            The protein sequence
        """
        return {"Instability-Index": ProteinAnalysis(seq).instability_index()}
    
    @staticmethod
    def _isoelectric_point(seq: Seq.Seq) -> dict:
        """
        Returns a dictionary containing the value of the protein's isoelectric point.
        
        Parameters
        ----------
        seq: Seq.Seq
            The protein sequence
        """
        return {"Isoelectric-Point": ProteinAnalysis(seq).isoelectric_point()}
    
    @staticmethod
    def _molar_extinction_coefficient(seq: Seq.Seq) -> dict:
        """
        Returns a dictionary containing the value of the protein's molar extinction
        coefficient (reduced and oxidized).
        
        Parameters
        ----------
        seq: Seq.Seq
            The protein sequence
        """
        reduced, oxidized = ProteinAnalysis(seq).molar_extinction_coefficient()
        return {"Reduced-MEC": reduced, "Oxidized-MEC": oxidized}
    
    @staticmethod
    def _secondary_structure_fraction(seq: Seq.Seq) -> dict:
        """
        Returns a dictionary containing the fraction of amino acids which tend to be
        in the secondary structures:
        - alpha-helix: Val (V), Ile (I), Tyr (Y), Phe (F), Trp (W), Leu (L)
        - beta-turn: Asn (N), Pro (P), Gly (G), Ser (S)
        - beta-sheet: Glu (E), Met (M), Ala (A), Leu (L)
        
        Parameters
        ----------
        seq: Seq.Seq
            The protein sequence
        """
        helix, turn, sheet = ProteinAnalysis(seq).secondary_structure_fraction()
        return {"Helix-SSF": helix, "Turn-SSF": turn, "Sheet-SSf": sheet}
    
    @staticmethod
    def _ctd_descriptors(seq: Seq.Seq) -> dict:
        """
        Returns a dictionary containing the CTD descriptors (composition, transition,
        distribution) of the protein sequence.
        
        Parameters
        ----------
        seq: Seq.Seq
            The protein sequence
        """
        return PyPro.GetProDes(seq).GetCTD()
    
    # maybe do not consider these descriptors (to long to compute)
    # ---
    @staticmethod
    def _geary_autocorrelation(seq: Seq.Seq) -> dict:
        """
        Returns a dictionary containing the Geary autocorrelation descriptors of the
        protein sequence it takes as input.
        
        Parameters
        ----------
        seq: Seq.Seq
            The protein sequence
        """
        return PyPro.GetProDes(seq).GetGearyAuto()
    
    # maybe do not consider these descriptors (to long to compute)
    # ---
    @staticmethod
    def _moran_autocorrelation(seq: Seq.Seq) -> dict:
        """
        Returns a dictionary containing the Moran autocorrelation descriptors of the
        protein sequence it takes as input.
        
        Parameters
        ----------
        seq: Seq.Seq
            The protein sequence
        """
        return PyPro.GetProDes(seq).GetMoranAuto()
    
    # maybe do not consider these descriptors (to long to compute)
    # ---
    @staticmethod
    def _moreau_broto_autocorrelation(seq: Seq.Seq) -> dict:
        """
        Returns a dictionary containing the Moreau-Broto autocorrelation descriptors of
        the protein sequence it takes as input.
        
        Parameters
        ----------
        seq: Seq.Seq
            The protein sequence
        """
        return PyPro.GetProDes(seq).GetMoreauBrotoAuto()

    def build_dataset(self) -> pd.DataFrame:
        """
        Returns a tabular dataset (represented by a pd.DataFrame) containing the featurized
        sequences (present in <file>).
        ---
        Feature mapping:
        - sequence length (1 feature), nucleotide composition (4 features), aminoacid
        composition (20 features), dipeptide composition (400 features), molecular weight
        (1 feature), aromaticity (1 feature), instability index (1 feature), isoelectric
        point (1 feature), molar extinctioncoefficient (2 features), secondary structure
        fraction (3 features), ctd descriptors (147 features), geary autocorrelation
        descriptors (240 features), moran autocorrelation descriptors (240 features),
        moreau-broto autocorrelation descriptors (240 features)
        ---
        Total number of features: 581 (1301 including autocorrelation descriptors)
        """
        # initialize dataset
        dataset = []
        # read sequences present in <file>
        sequences = SeqIO.parse(self.file, format="fasta")
        # iterate through sequences and build dataset
        for i, seq in enumerate(sequences):
            # translate DNA sequence
            tseq = self._translate_seq(seq.seq)
            # probably temporary -> get rid of unwanted symbols in the sequence
            #                       it may be kept (there are few problematic sequences)
            tseq = MLDataset._get_valid_protein(tseq)
            # build entry by computing sequence descriptors
            entry = {"Len-Protein": len(tseq),
                     **MLDataset._nuc_composition(seq.seq),
                     **MLDataset._amino_composition(tseq),
                     **MLDataset._dipeptide_composition(tseq),
                     **MLDataset._molecular_weight(tseq),
                     **MLDataset._aromaticity(tseq),
                     **MLDataset._instability_index(tseq),
                     **MLDataset._isoelectric_point(tseq),
                     **MLDataset._molar_extinction_coefficient(tseq),
                     **MLDataset._secondary_structure_fraction(tseq),
                     **MLDataset._ctd_descriptors(tseq)} #,
                     #**MLDataset._geary_autocorrelation(tseq),
                     #**MLDataset._moran_autocorrelation(tseq),
                     #**MLDataset._moreau_broto_autocorrelation(tseq)}
            # add "entry" to "dataset"
            dataset.append(entry)
        # build pd.DataFrame from "dataset"
        df_out = pd.DataFrame(dataset)
        # add label column to the dataframe
        df_out["Function"] = [self.prot_name] * (i+1)
        # return df containing the featurized records
        return df_out
        