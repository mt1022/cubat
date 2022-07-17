import re
import time
import numpy as np
import pandas as pd
from Bio import SeqIO
# from Bio.Seq import Seq
import genecode_data
from Bio.Data import CodonTable
from Bio.Seq import Seq

_genecode = 1
_reference_file = None
memories = {}

pd.options.display.max_columns = None

# reference file
ref_heg = [544814.0, 655524.0, 234253.0, 406803.0, 477158.0, 567853.0, 383463.0, 143585.0, 385678.0, 495846.0,
           327792.0, 403140.0, 413535.0, 632996.0, 224870.0, 1295575.0, 558133.0, 645805.0, 538723.0, 226405.0,
           340239.0, 483801.0, 382585.0, 1109029.0, 148414.0, 345088.0, 200253.0, 376554.0, 502403.0, 677714.0,
           414756.0, 615448.0, 476341.0, 197169.0, 532711.0, 615538.0, 764918.0, 1035170.0, 384496.0, 626669.0,
           373118.0, 375881.0, 349016.0, 470297.0, 223912.0, 919492.0, 595935.0, 912370.0, 510981.0, 242826.0,
           702386.0, 820851.0, 922616.0, 1291190.0, 345134.0, 727678.0, 524128.0, 531888.0, 17361.0, 20594.0,
           36263.0]
ref_bg = [687589.0, 512749.0, 345215.0, 295841.0, 459081.0, 342347.0, 414992.0, 355639.0, 504963.0, 376561.0,
          418699.0, 312233.0, 749622.0, 559009.0, 677630.0, 580714.0, 575017.0, 428803.0, 519794.0, 445452.0,
          472034.0, 352006.0, 803248.0, 688366.0, 312558.0, 233081.0, 282540.0, 242131.0, 532065.0, 396773.0,
          497528.0, 371017.0, 449746.0, 385423.0, 657750.0, 490499.0, 969364.0, 830724.0, 579225.0, 431940.0,
          403343.0, 345656.0, 573163.0, 427420.0, 518118.0, 444016.0, 660594.0, 492619.0, 597152.0, 511746.0,
          872554.0, 650683.0, 1192155.0, 1021651.0, 621672.0, 463594.0, 561968.0, 481594.0, 23436.0, 27347.0,
          23436.0]
fop_opt = np.array([2, 5, 6, 8, 11, 12, 17, 18, 21, 28, 29, 34, 36, 39, 47, 50, 53, 57])
cbi_opt = np.array(
    [0, 2, 5, 6, 8, 11, 12, 15, 17, 18, 21, 26, 28, 29, 32, 34, 36, 37, 39, 44, 46, 47, 50, 51, 53, 55, 57])
tGCN = [12, 8, 0, 6, 6, 9, 1, 6, 5, 5, 3, 4, 6, 8, 4, 7, 15, 8, 0, 13, 5, 5, 7, 4, 20, 13, 7, 9, 5, 4, 4, 4, 0, 0, 1,
         0, 0, 0, 0, 7, 15, 10, 0, 9, 9, 26, 9, 9, 25, 13, 13, 9, 8, 14, 29, 0, 3, 0, 10, 0, 0, 0, 0, 0]
reference_set=[510981, 912370, 242826, 595935, 403140, 327792, 820851, 702386, 922616, 1291190, 655524, 544814, 524128, 727678, 531888, 345134, 483801, 340239, 229688, 677714, 502403, 764918, 1035170, 224870, 632996, 1295575, 413535, 234253, 406803, 710465, 615538, 532711, 538723, 645805, 226405, 558133, 382585, 1109029, 373118, 375881, 200253, 345088, 376554, 148414, 626669, 384496, 383463, 567853, 143585, 477158, 476341, 615448, 197169, 414756, 223912, 470297, 919492, 349016, 417781, 495846, 385678, 20594, 17361, 36263]
reference_order=['GCA', 'GCC', 'GCG', 'GCT', 'TGC', 'TGT', 'GAC', 'GAT', 'GAA', 'GAG', 'TTC', 'TTT', 'GGA', 'GGC', 'GGG', 'GGT', 'CAC', 'CAT', 'ATA', 'ATC', 'ATT', 'AAA', 'AAG', 'CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG', 'ATG', 'AAC', 'AAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CAA', 'CAG', 'AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT', 'AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT', 'ACA', 'ACC', 'ACG', 'ACT', 'GTA', 'GTC', 'GTG', 'GTT', 'TGG', 'TAC', 'TAT', 'TAA', 'TAG', 'TGA']



def set_genecode(genecode):
    _genecode = genecode
    return print('Genecode has been set')


def input_reference_file(file_path):
    _reference_file = file_path
    return print('Reference file has been set')


def generate_reference_file():
    reference_dataframe=pd.DataFrame(index=['TTT', 'TTC', 'TTA', 'TTG',
                                                'TCT', 'TCC', 'TCA', 'TCG',
                                                'TAT', 'TAC', 'TAA', 'TAG',
                                                'TGT', 'TGC', 'TGA', 'TGG',
                                                'CTT', 'CTC', 'CTA', 'CTG',
                                                'CCT', 'CCC', 'CCA', 'CCG',
                                                'CAT', 'CAC', 'CAA', 'CAG',
                                                'CGT', 'CGC', 'CGA', 'CGG',
                                                'ATT', 'ATC', 'ATA', 'ATG',
                                                'ACT', 'ACC', 'ACA', 'ACG',
                                                'AAT', 'AAC', 'AAA', 'AAG',
                                                'AGT', 'AGC', 'AGA', 'AGG',
                                                'GTT', 'GTC', 'GTA', 'GTG',
                                                'GCT', 'GCC', 'GCA', 'GCG',
                                                'GAT', 'GAC', 'GAA', 'GAG',
                                                'GGT', 'GGC', 'GGA', 'GGG'])
    reference_dataframe['']


def info_preprocess(parse_data):
    data = list(parse_data)

    id_list = []
    sequence_length_list = []
    frequency_check = []
    start_codon_check = []
    seq_list = []

    for record in data:
        id_list.append(record.name)
        sequence_length_list.append(len(record.seq))
        start_codon_check.append(str(record.seq[0:3]))
        seq_list.append(str(record.seq))
        if len(record.seq) % 3 == 0:
            frequency_check.append("True")
        else:
            frequency_check.append("False")

    info_tuple = (id_list, seq_list, sequence_length_list, frequency_check, start_codon_check)

    return info_tuple


def count_codon(sequence):
    codon_num = {'TTT': 0, 'TTC': 0, 'TTA': 0, 'TTG': 0, 'TCT': 0, 'TCC': 0, 'TCA': 0, 'TCG': 0, 'TAT': 0, 'TAC': 0,
                 'TAA': 0, 'TAG': 0, 'TGT': 0,
                 'TGC': 0, 'TGA': 0, 'TGG': 0, 'CTT': 0, 'CTC': 0, 'CTA': 0, 'CTG': 0, 'CCT': 0, 'CCC': 0, 'CCA': 0,
                 'CCG': 0,
                 'CAT': 0,
                 'CAC': 0, 'CAA': 0, 'CAG': 0, 'CGT': 0, 'CGC': 0, 'CGA': 0, 'CGG': 0, 'ATT': 0, 'ATC': 0, 'ATA': 0,
                 'ATG': 0,
                 'ACT': 0, 'ACC': 0, 'ACA': 0, 'ACG': 0, 'AAT': 0, 'AAC': 0, 'AAA': 0, 'AAG': 0, 'AGT': 0, 'AGC': 0,
                 'AGA': 0,
                 'AGG': 0, 'GTT': 0, 'GTC': 0, 'GTA': 0, 'GTG': 0, 'GCT': 0, 'GCC': 0, 'GCA': 0, 'GCG': 0, 'GAT': 0,
                 'GAC': 0,
                 'GAA': 0, 'GAG': 0, 'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0}

    depart = len(sequence) / 8
    times = int(depart)
    rest = int((depart - times) * 8)
    for i in range(times):
        codon_num[sequence.pop()] += 1
        codon_num[sequence.pop()] += 1
        codon_num[sequence.pop()] += 1
        codon_num[sequence.pop()] += 1
        codon_num[sequence.pop()] += 1
        codon_num[sequence.pop()] += 1
        codon_num[sequence.pop()] += 1
        codon_num[sequence.pop()] += 1
    for i in range(rest):
        codon_num[sequence.pop()] += 1

    codon_array = np.array(list(codon_num.values()))

    return codon_array


def frequency_count(seq_data):
    if type(seq_data) == tuple:
        seq_data = list(seq_data)
    times = len(seq_data)
    result_array = np.zeros((times, 64), dtype=int)
    for t in range(0, times):
        seq_str = seq_data[t]
        seq_list = re.findall('.{3}', seq_str)
        result_array[t] = count_codon(seq_list)

    return result_array


def rscu_compute(codon_dataframe,genecode):
    stop_codon = []
    rscu_dataframe = pd.DataFrame(columns=codon_dataframe.columns, index=codon_dataframe.index)
    for codon in codon_dataframe.columns:
        try:
            amino_acid_value = CodonTable.unambiguous_dna_by_id[genecode].forward_table[codon]
            rscu_codons = \
                [key for key, value in CodonTable.unambiguous_dna_by_id[genecode].forward_table.items() if
                 value == amino_acid_value]
            rscu_dataframe[codon] = codon_dataframe[codon] / (
                        codon_dataframe[rscu_codons].sum(axis=1) / len(rscu_codons))

        except KeyError:
            stop_codon.append(codon)
    return rscu_dataframe


def enc_compute(frequency_array):
    # call up the database
    familys = genecode_data.all_24families[_genecode]
    codon_location = genecode_data.codon_family24[_genecode]
    amount_of_eachfamily = genecode_data.number_of_eachfamily[_genecode]
    amount_of_single = genecode_data.number_of_single

    # judge that whether the frequency_array was a tuple.
    # tips:look that of the "return" of def process_data. the result is in a form of tuple
    if type(frequency_array) == tuple:
        frequency_array = frequency_array[0]

    number_of_sequence = frequency_array.shape[0]
    process_arr = np.delete(frequency_array, genecode_data.delete_location[_genecode], axis=1)
    FcF_array, an_array = np.empty((21, number_of_sequence), dtype=float), np.empty((21, number_of_sequence),
                                                                                    dtype=float)

    # do the computing
    number_of_families = len(familys)
    num1xa = process_arr
    for i in range(0, number_of_families):
        F_array = num1xa[:, codon_location[i]:codon_location[i + 1]]
        an = np.sum(F_array, axis=1)
        an_array[i] = an
        FcF_array[i] = np.sum(np.square(F_array + 1), axis=1) / np.square(an + familys[i])

    family2, family3, family4 = (familys == 2), (familys == 3), (familys == 4)
    family2sum = amount_of_eachfamily[0] * np.sum(an_array[family2], axis=0) / np.sum(
        np.multiply(FcF_array[family2], an_array[family2]), axis=0)
    family3_prevent0 = np.sum(np.multiply(FcF_array[family3], an_array[family3]), axis=0)
    family3_prevent0[family3_prevent0 == 0] = 1
    family3sum = amount_of_eachfamily[1] * np.sum(an_array[family3], axis=0) / family3_prevent0
    family4sum = amount_of_eachfamily[2] * np.sum(an_array[family4], axis=0) / np.sum(
        np.multiply(FcF_array[family4], an_array[family4]), axis=0)

    enc_array = amount_of_single[_genecode] + family2sum + family3sum + family4sum

    return enc_array


def ite_compute(frequency_array):
    times = frequency_array.shape[0]
    ite_array = np.empty(times, dtype=float)
    Wte = np.array(np.zeros((29, 2), dtype=float))
    for i in range(0, 58, 2):
        Wte[int(i / 2), 0] = ref_heg[i] / ref_bg[i]
        Wte[int(i / 2), 1] = ref_heg[i + 1] / ref_bg[i + 1]
        Wte[int(i / 2)] *= 1 / max([ref_heg[i] / ref_bg[i], ref_heg[i + 1] / ref_bg[i + 1]])
    Wte = Wte.flatten()
    for t in range(0, times):
        num = frequency_array[t]
        num2 = np.delete(num, [10, 11, 14, 15, 35, 36])
        mut = np.multiply(np.log(Wte), num2)
        sum1 = np.sum(mut)
        ITE = np.exp(sum1 / np.sum(num))
        ite_array[t] = ITE
    return ite_array


def chi_squared(processed_array):
    codon_familyami = genecode_data.codon_familyami[_genecode]
    all_amifamilies = genecode_data.all_amifamilies[_genecode]

    times = processed_array.shape[0]
    all_ami = np.zeros((18, times), dtype=float)
    for i in range(0, len(all_amifamilies)):
        codon_count = processed_array[:, codon_familyami[i]:codon_familyami[i + 1]]
        ami_count = np.sum(codon_count, axis=1)
        ami_count[ami_count == 0] = 1
        every_codon = np.square((codon_count.T / ami_count) - (1 / all_amifamilies[i]))
        all_codon = every_codon * ami_count * all_amifamilies[i]
        all_ami[i] = np.sum(all_codon, axis=0)
    sum_arr = np.sum(all_ami, axis=0)
    chi2_array = sum_arr / np.sum(processed_array, axis=1)

    return chi2_array


def cbi_compute(processed_array, optimized_codon=cbi_opt):
    codon_familyami = genecode_data.codon_familyami[_genecode]
    all_amifamilies = genecode_data.all_amifamilies[_genecode]

    all_codon = np.sum(processed_array, axis=1)
    opt_codon = np.sum(processed_array[:, optimized_codon], axis=1)
    rand_number = np.empty(processed_array.shape[0], dtype=float)
    cbi_list = list(optimized_codon)

    for i in range(0, len(codon_familyami) - 1):
        opt_ami = 0
        one_ami = processed_array[:, codon_familyami[i]:codon_familyami[i + 1]]
        ami_count = np.sum(one_ami, axis=1)
        while cbi_list[0] < codon_familyami[i + 1]:
            cbi_list.pop(0)
            opt_ami += 1
            if len(cbi_list) == 0:
                break
        weight = opt_ami / all_amifamilies[i]
        rand_number += ami_count * weight

    cbi_array = (opt_codon - rand_number) / (all_codon - rand_number)
    return cbi_array


def chi_squared(processed_array):
    codon_location_ami = genecode_data.codon_familyami[_genecode]
    family = genecode_data.all_amifamilies[_genecode]
    times = processed_array.shape[0]
    all_ami = np.zeros((18, times), dtype=float)

    for i in range(0, len(family)):
        codon_count = processed_array[:, codon_location_ami[i]:codon_location_ami[i + 1]]
        ami_count = np.sum(codon_count, axis=1)
        ami_count[ami_count == 0] = 1
        every_codon = np.square((codon_count.T / ami_count) - (1 / family[i]))
        all_codon = every_codon * ami_count * family[i]
        all_ami[i] = np.sum(all_codon, axis=0)
    sum_arr = np.sum(all_ami, axis=0)
    X2_array = sum_arr / np.sum(processed_array, axis=1)

    return X2_array


def fop_compute(processed_array, optimized_codon=fop_opt):
    opt_array = processed_array[:, optimized_codon]
    fop_array = np.sum(opt_array, axis=1) / np.sum(processed_array, axis=1)

    return fop_array

def tai_compute(frequency_array, SI_C=0.28, SI_A=0.9999, SG_U=0.41, SU_G = 0.68):
    times=frequency_array.shape[0]
    tai_array = np.empty(times, dtype=float)
    list2 = []
    tRNA=['TTT', 'TTC', 'TTA', 'TTG', 'TCT', 'TCC', 'TCA', 'TCG', 'TAT', 'TAC', 'TAG', 'TAA', 'TGT', 'TGC', 'TGA', 'TGG', 'CTT', 'CTC', 'CTA', 'CTG', 'CCT', 'CCC', 'CCA', 'CCG', 'CAT', 'CAC', 'CAA', 'CAG', 'CGT', 'CGC', 'CGA', 'CGG', 'ATT', 'ATC', 'ATA', 'ATG', 'ACT', 'ACC', 'ACA', 'ACG', 'AAT', 'AAC', 'AAA', 'AAG', 'AGT', 'AGC', 'AGA', 'AGG', 'GTT', 'GTC', 'GTA', 'GTG', 'GCT', 'GCC', 'GCA', 'GCG', 'GAT', 'GAC', 'GAA', 'GAG', 'GGT', 'GGC', 'GGA', 'GGG']
    for i in tRNA:
        list2.append(tGCN[tRNA.index(Seq(i).reverse_complement())])
    Ws = np.array(np.zeros((1, 64), dtype=float))
    for i in range(0, 61, 4):
        Ws[0, i] = 1 * list2[i] + (1 - SG_U) * list2[i + 1]
        Ws[0, i + 1] = 1 * list2[i + 1] + (1 - SI_C) * list2[i]
        Ws[0, i + 2] = 1 * list2[i + 2] + (1 - SI_A) * list2[i]
        Ws[0, i + 3] = 1 * list2[i + 3] + (1 - SU_G) * list2[i + 2]
    Ws = np.delete(Ws, [10, 11, 14, 35], axis=1)
    Ws[Ws == 0] = np.mean(Ws != 0)
    Ws = Ws / Ws.max()
    for m in range(0,times):
        num = frequency_array[m]
        num = np.delete(num, [10, 11, 14, 35])
        mut = np.multiply(np.log(Ws), num)
        sum1 = np.sum(mut, axis=1)
        sum0 = np.sum(sum1, axis=0)
        tAI = np.exp(sum0 / sum(num))
        tai_array[m]=tAI
    return tai_array

def codon_table_completion(genecode):
    codon_table=CodonTable.unambiguous_dna_by_id[genecode].forward_table
    codons = ["TTT", "TTC", "TTA", "TTG",
              "TCT", "TCC", "TCA", "TCG",
              "TAT", "TAC", "TAA", "TAG",
              "TGT", "TGC", "TGA", "TGG",
              "CTT", "CTC", "CTA", "CTG",
              "CCT", "CCC", "CCA", "CCG",
              "CAT", "CAC", "CAA", "CAG",
              "CGT", "CGC", "CGA", "CGG",
              "ATT", "ATC", "ATA", "ATG",
              "ACT", "ACC", "ACA", "ACG",
              "AAT", "AAC", "AAA", "AAG",
              "AGT", "AGC", "AGA", "AGG",
              "GTT", "GTC", "GTA", "GTG",
              "GCT", "GCC", "GCA", "GCG",
              "GAT", "GAC", "GAA", "GAG",
              "GGT", "GGC", "GGA", "GGG"]
    complete_codon_table = codon_table
    for codon in codons:
        try:
            codon_table[codon]
        except KeyError:
            complete_codon_table.update({codon: "*"})
    return complete_codon_table


def generate_wij_dataframe(codon_dataframe):
    wij_dict = {}
    for codon in codon_dataframe.columns:
        amino_acid_value = codon_table_completion(_genecode)[codon]
        wij_codons = \
            [key for key, value in codon_table_completion(_genecode).items() if
             value == amino_acid_value]
        wij = \
            codon_dataframe[wij_codons].sum()[codon]/codon_dataframe[wij_codons].sum().sort_values(ascending=False)[0]
        wij_dict.update({codon: wij})
    return pd.DataFrame(pd.Series(wij_dict))


def cai_compute(codon_dataframe, reference):
    reset_reference_dataframe = reference.loc[codon_dataframe.columns]
    reset_reference_dataframe.columns = ["wij"]

    codon_dataframe.loc["All sequences"] = codon_dataframe.apply(lambda x: x.sum())
    cai_up_dataframe = pd.DataFrame(columns=codon_dataframe.columns, index=codon_dataframe.index)
    cai2_up_dataframe = pd.DataFrame(columns=codon_dataframe.columns, index=codon_dataframe.index)

    for codon in codon_dataframe.columns:
        cai_up_dataframe[codon] = codon_dataframe[codon] * np.log(reset_reference_dataframe.loc[codon, "wij"])
        cai2_up_dataframe[codon] = codon_dataframe[codon] * reset_reference_dataframe.loc[codon, "wij"]

    cai_up_dataframe["Total"] = cai_up_dataframe.apply(lambda x: x.sum(), axis=1)
    cai2_up_dataframe["Total"] = cai2_up_dataframe.apply(lambda x: x.sum(), axis=1)
    codon_dataframe["Total frequency"] = codon_dataframe.apply(lambda x: x.sum(), axis=1)
    cai_value = cai_up_dataframe["Total"] / codon_dataframe["Total frequency"]
    cai2_value = cai2_up_dataframe["Total"] / codon_dataframe["Total frequency"]

    for key in cai_value.keys():
        cai_value[key] = np.exp(cai_value[key])

    return (cai_value, cai2_value)

def csc_compute(codon_dataframe, mrna_hl_location):
    mrna_hl = pd.read_excel(mrna_hl_location)  # Half-life of mrna
    array_mrna = mrna_hl[["Total Half-life"]].values.T
    csc_dataframe = pd.DataFrame(columns=codon_dataframe.columns, index=["csc"])
    for codon in csc_dataframe:
        csc_dataframe[codon] = float(
            np.corrcoef(array_mrna, codon_dataframe[codon].values.T)[[0], [1]])
    return csc_dataframe


def codon_bias_analyse(file_path, file_format='fasta', quality_control=True, rscu=True, enc=True, ite=True, X2=True,
                       fop=True, cbi=True,tai=True,cai=True,csc=False):
    delete_location = genecode_data.delete_location[_genecode]
    sorting = genecode_data.sorting[_genecode]

    file = SeqIO.parse(file_path, file_format)
    information = info_preprocess(file)
    frequency_array = frequency_count(information[1])
    frequency_dataframe = pd.DataFrame(frequency_array,
                                       columns=['TTT', 'TTC', 'TTA', 'TTG',
                                                'TCT', 'TCC', 'TCA', 'TCG',
                                                'TAT', 'TAC', 'TAA', 'TAG',
                                                'TGT', 'TGC', 'TGA', 'TGG',
                                                'CTT', 'CTC', 'CTA', 'CTG',
                                                'CCT', 'CCC', 'CCA', 'CCG',
                                                'CAT', 'CAC', 'CAA', 'CAG',
                                                'CGT', 'CGC', 'CGA', 'CGG',
                                                'ATT', 'ATC', 'ATA', 'ATG',
                                                'ACT', 'ACC', 'ACA', 'ACG',
                                                'AAT', 'AAC', 'AAA', 'AAG',
                                                'AGT', 'AGC', 'AGA', 'AGG',
                                                'GTT', 'GTC', 'GTA', 'GTG',
                                                'GCT', 'GCC', 'GCA', 'GCG',
                                                'GAT', 'GAC', 'GAA', 'GAG',
                                                'GGT', 'GGC', 'GGA', 'GGG'],
                                       index=np.array(information[0]))
    print(frequency_dataframe)
    index_dataframe = pd.DataFrame(index=np.array(information[0]))
    index_dataframe['sequnece_length'] = information[2]

    if quality_control:
        index_dataframe['mult_of_3'] = information[3]
        index_dataframe['start_codon'] = information[4]

    if enc:
        enc_array = enc_compute(frequency_array)
        index_dataframe['enc'] = enc_array

    if ite:
        ite_array = ite_compute(frequency_array)
        index_dataframe['ite'] = ite_array

    if tai:
        tai_array=tai_compute(frequency_array)
        index_dataframe['tai']=tai_array

    processed_array = np.delete(frequency_array, delete_location, axis=1)[:, sorting]

    if X2:
        X2_array = chi_squared(processed_array)
        index_dataframe['X2'] = X2_array

    if fop:
        fop_array = fop_compute(processed_array)
        index_dataframe["fop"] = fop_array

    if cbi:
        cbi_array = cbi_compute(processed_array)
        index_dataframe['cbi'] = cbi_array


    if cai:
        reference_dataframe = pd.DataFrame(columns=reference_order)
        reference_dataframe.loc[len(reference_dataframe)] = np.array(reference_set)
        reference_weight=generate_wij_dataframe(reference_dataframe)
        # reference_dataframe=pd.read_excel('Homo_sapiens.xlsx')
        cai_dataframe,cai2_dataframe=cai_compute(frequency_dataframe,reference_weight)
        index_dataframe['cai']=cai_dataframe
        index_dataframe['cai2']=cai2_dataframe

    print(index_dataframe)


    if rscu:
        rscu_dataframe = rscu_compute(frequency_dataframe,_genecode)
        rscu_dataframe=rscu_dataframe.drop(columns='Total frequency')
        print(rscu_dataframe)

    if csc:
        csc_dataframe = csc_compute(frequency_dataframe,mrna_hl_location=None)
        print(csc_dataframe)

    return


set_genecode(2)
lll = codon_bias_analyse('Homo_sapiens.GRCh38.cds2.FAS')
