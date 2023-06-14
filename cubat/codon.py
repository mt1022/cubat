# Available Genetic Codes (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi):
# 1. The Standard Code
# 2. The Vertebrate Mitochondrial Code
# 3. The Yeast Mitochondrial Code
# 4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
# 5. The Invertebrate Mitochondrial Code
# 6. The Ciliate, Dasycladacean and Hexamita Nuclear Code
# 9. The Echinoderm and Flatworm Mitochondrial Code
# 10. The Euplotid Nuclear Code
# 11. The Bacterial, Archaeal and Plant Plastid Code
# 12. The Alternative Yeast Nuclear Code
# 13. The Ascidian Mitochondrial Code
# 14. The Alternative Flatworm Mitochondrial Code
# 16. Chlorophycean Mitochondrial Code
# 21. Trematode Mitochondrial Code
# 22. Scenedesmus obliquus Mitochondrial Code
# 23. Thraustochytrium Mitochondrial Code
# 24. Rhabdopleuridae Mitochondrial Code
# 25. Candidate Division SR1 and Gracilibacteria Code
# 26. Pachysolen tannophilus Nuclear Code
# 27. Karyorelict Nuclear Code
# 28. Condylostoma Nuclear Code
# 29. Mesodinium Nuclear Code
# 30. Peritrich Nuclear Code
# 31. Blastocrithidia Nuclear Code
# 33. Cephalodiscidae Mitochondrial UAA-Tyr Code

from Bio.Data import CodonTable

# codons
CODONS = [i + j + k for i in 'ACGT' for j in 'ACGT' for k in 'ACGT']


def get_codon_table(table=1):
    """
    Get Genetic Code by ID

    see available tables: `CodonTable.generic_by_id`
    """
    codon_table = CodonTable.generic_by_id[table]

    # insert stop codons
    c2a = codon_table.forward_table
    for codon in codon_table.stop_codons:
        c2a[codon] = '*'
    
    return {codon: c2a[codon] for codon in CODONS}


def read_custom_table(path):
    """
    Read custom codon table from a csv file
    """
    codon_table = dict()
    with open(path, 'rt') as fh:
        for line in fh:
            ary = line.rstrip().split(',')
            codon_table[ary[0]] = ary[1]
    for codon in CODONS:
        if codon not in codon_table:
            exit(f'Codon Table is not complete: {codon}')
    return codon_table


class GeneticCode():

    def __init__(self, type, table=None, path=None):
        if type == 'ncbi':
            self.codon_table = get_codon_table(table)
        elif type == 'custom':
            self.codon_table = read_custom_table(path)
        else:
            raise ValueError(f'type must be one of "ncbi" or "custom"')
        
        self.codon_weights = pd.DataFrame()
    
    def aa_to_codon(self):
        """
        Return a Dict of AA to codons mappings
        """
        a2c = dict()
        for codon, aa in self.codon_table.items():
            if aa not in a2c:
                a2c[aa] = []
            a2c[aa].append(codon)
        return a2c
