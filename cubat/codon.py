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


def get_codon_table(id=1):
    """
    Get Genetic Code by ID

    see available tables: `CodonTable.generic_by_id`
    """
    codon_table = CodonTable.generic_by_id[id]

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
