# -*- coding: utf-8 -*-

import numpy as np
import os
import pandas as pd
import time
import warnings
from Bio import Entrez, SeqIO, Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from propy import PyPro
from tqdm.auto import tqdm

# ignore Biopyhon warnings (related to sequence translation)
warnings.filterwarnings("ignore")

class MLDataset:
    
    """
    Constructs a tabular dataset by computing a set of features/descriptors associated to
    the DNA sequences present in the input .fasta file.
    """
    
    # INIT & REPR
    # ---
    
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
            
        Attributes
        ----------
        _bad_inds: list
            A list object storing indices of unfeaturized sequence records
        """
        # parameters
        self.file = file
        self.prot_name = prot_name
        self.ttable = ttable
        self.icodons = icodons
        # attributes
        self._bad_inds = []
        
    def __repr__(self) -> str:
        """
        Returns the string representation of the object.
        """
        c = self.__class__.__name__
        r = f"{c}({self.file!r}, {self.prot_name!r}, {self.ttable!r}, {self.icodons!r})"
        return r
    
    # TRANSLATE AND VALIDATE SEQUENCE
    # ---
    
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
    def _get_valid_protein(tseq: Seq.Seq) -> Seq.Seq:
        """
        Returns a valid protein sequence from the sequence it takes as input (substitutes
        invalid characters for the empty string).
        
        Parameters
        ----------
        seq: Seq.Seq
            The protein sequence to be validated
        """
        non_valid_chars = "*BJOUXZ"
        for c in non_valid_chars:
            tseq = tseq.replace(c, "")
        return tseq
    
    # COMPUTE SEQUENCE DESCRIPTORS
    # ---
            
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
    
    # BUILD DATASET ENTRY FROM SEQUENCE
    # ---
        
    @staticmethod
    def _build_entry(seq: Seq.Seq, tseq: Seq.Seq) -> dict:
        """
        Given a sequence of aminoacids, it builds and returns an entry to be added to the
        dataset of featurized sequences.
        
        Parameters
        ----------
        seq: Seq.Seq
            The DNA sequence which translates to the protein sequence
        tseq: Seq.Seq
            The sequence of aminoacids
        """
        return {"Len-Protein": len(tseq),
                **MLDataset._nuc_composition(seq.seq),
                **MLDataset._amino_composition(tseq),
                **MLDataset._dipeptide_composition(tseq),
                **MLDataset._molecular_weight(tseq),
                **MLDataset._aromaticity(tseq),
                **MLDataset._instability_index(tseq),
                **MLDataset._isoelectric_point(tseq),
                **MLDataset._molar_extinction_coefficient(tseq),
                **MLDataset._secondary_structure_fraction(tseq),
                **MLDataset._ctd_descriptors(tseq)}
        
    # ADD ENTRY TO DATASET AND BUILD FINAL DATAFRAME
    # ---
        
    @staticmethod
    def _add_entry(dataset: pd.DataFrame, entry: dict, id_seq: int) -> pd.DataFrame:
        """
        Adds an entry <entry> to the dataset <dataset> at index <id_seq>. The pd.DataFrame
        parameter "dtype" is set to np.float32 to reduce RAM usage.
        
        Parameters
        ----------
        dataset: pd.DataFrame
            A dataset of featurized sequences
        entry: dict
            The entry to be added to <dataset>
        id_seq: int
            The index at which the new entry is to be added
        """
        idx = range(id_seq, id_seq+1)
        df_entry = pd.DataFrame(entry, index=idx, dtype=np.float32)
        new_df = pd.concat([dataset, df_entry])
        return new_df
    
    # BUILD_DATASET
    # ---
            
    def build_dataset(self) -> pd.DataFrame:
        """
        Returns a tabular dataset (represented by a pd.DataFrame) containing the featurized
        sequences (present in <file>).
        ---
        Feature mapping:
        - sequence length (1 feature), nucleotide composition (4 features), aminoacid
        composition (20 features), dipeptide composition (400 features), molecular weight
        (1 feature), aromaticity (1 feature), instability index (1 feature), isoelectric
        point (1 feature), molar extinction coefficient (2 features), secondary structure
        fraction (3 features), ctd descriptors (147 features)
        ---
        Total number of features: 581
        """
        # initialize <datasets> and <dataset>
        datasets = []
        dataset = pd.DataFrame()
        # read sequences present in <file>
        sequences = SeqIO.parse(self.file, format="fasta")
        # iterate through sequences and build dataset
        for i, seq in tqdm(enumerate(sequences)):
            # append <dataset> to <datasets> and reinitialize <dataset> every 5000 sequences
            # (speedup the concatenation of dataset and entry -> the largest the dataset,
            # the longer it takes to concatenate, since a copy of it has to be created)
            if i % 5000 == 0:
                datasets.append(dataset)
                dataset = pd.DataFrame()
            # translate DNA sequence and get rid of non valid symbols (*BJOUXZ)
            tseq = self._translate_seq(seq.seq)
            tseq = MLDataset._get_valid_protein(tseq=tseq)
            # build entry by computing sequence descriptors and add it to <dataset>
            try:
                entry = MLDataset._build_entry(seq, tseq)
            except:
                self._bad_inds.append(i)
                name_prot = " | ".join(seq.description.split(" | ")[2:])
                print(f"Featurization failed for {name_prot!r} sequence...")
            else:    
                dataset = MLDataset._add_entry(dataset, entry, i)
        # build pd.DataFrame from the list of datasets
        df_out = pd.concat(datasets + [dataset])
        # add label column to the dataframe and return it
        df_out["Function"] = [self.prot_name] * (i+1 - len(self._bad_inds))
        return df_out
        