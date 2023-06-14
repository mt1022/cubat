import sys
import click
from codon import CODONS, GeneticCode
from seq import CDS, CDSDict


def read_cds(fasta):
    """
    Read and QC of CDS sequences. Return a CDSDict object.
    """