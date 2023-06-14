import pandas as pd
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from codon import CODONS


class CDS(SeqRecord):

    def __init__(self, seq, id):
        super().__init__(seq=seq, id=id)

    def translate(self, table=1):
        """Translates the CDS sequence to a protein sequence."""
        protein = self.seq.translate(table=table)
        return protein
    
    def get_codons(self):
        """
        Returns a list of codons in the DNA sequence.
        """
        codons = []
        for i in range(0, len(self.sequence), 3):
            codons.append(self.sequence[i:i + 3])
        return codons

    def calc_codon_freq(self):
        """
        Returns a Series of codon frequencies in the CDS.
        """
        codon_freq = {codon: 0 for codon in CODONS}
        for codon in self.get_codons():
            codon_freq[codon] += 1
        return pd.Series(codon_freq)


class SeqDict(defaultdict):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __setitem__(self, key, value):
        if not isinstance(value, SeqRecord):
                raise TypeError("value must be a SeqRecord")
        super().__setitem__(key, value)


class CDSDict(defaultdict):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __setitem__(self, key, value):
        if not isinstance(value, CDS):
                raise TypeError("value must be a SeqRecord")
        super().__setitem__(key, value)
    
    def calc_codon_freq(self):
        """
        Returns a pandas.DataFrame of codon frequencies in the CDSs.
        """
        codon_freq = pd.DataFrame([self[id].calc_codon_freq() for id in self])
        return codon_freq.transpose()
