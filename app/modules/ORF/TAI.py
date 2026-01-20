from collections import Counter
import numpy as np
from modules.shared_functions_and_vars import *


from tqdm import tqdm

def geometric_mean(iterable):
    a = np.log(iterable)
    return np.exp(a.mean())

def antiCodon_L(tgcn_dict):
    """This function get tRNA list of microorganism and return countered tGCN (tRNA Gene Copy Number)"""

    antiCodon_list = []
    for row in tqdm(tgcn_dict.keys()):
        antiCodon = row[-7:-4]
        antiCodon_list += [antiCodon]
    tGCN = Counter(antiCodon_list)
    return tGCN

def print_counter(c):
    for t in c:
        print(str(t) + '\t' + str(c[t]))
    return


def weight_cal(Sij, tGCN):
    """ get the sequence with Sij table and tGCN (from prev. func) and return list of codons' weight."""


    W_dict = {}
    codon_list = ['ATA', 'ATC', 'ATT', 'ATG', 'ACA', 'ACC', 'ACG', 'ACT', 'AAC', 'AAT', 'AAA', 'AAG', 'AGC', 'AGT', 'AGA', 'AGG', 'CTA', 'CTC', 'CTG', 'CTT', 'CCA', 'CCC', 'CCG', 'CCT', 'CAC', 'CAT', 'CAA', 'CAG', 'CGA', 'CGC', 'CGG', 'CGT', 'GTA', 'GTC', 'GTG', 'GTT', 'GCA', 'GCC', 'GCG', 'GCT', 'GAC', 'GAT', 'GAA', 'GAG', 'GGA', 'GGC', 'GGG', 'GGT', 'TCA', 'TCC', 'TCG', 'TCT', 'TTC', 'TTT', 'TTA', 'TTG', 'TAC', 'TAT', 'TGC', 'TGT', 'TGG', 'TAG', 'TGA', 'TAA']
    for codon in codon_list:
        W = 0
        for key in Sij.keys():
            if key[2] == codon[-1]:
                S = Sij[key]
                i_codon = reverse_complement(codon)
                anticodon_pos1 = key[0]
                anti_codon = anticodon_pos1 + i_codon[1:3]
                tGC_curr = tGCN.get(str(anti_codon))

                if not tGC_curr:
                    continue

                else:
                    W += (1-S)*tGC_curr

        W_dict[codon] = W

    factor = max(W_dict.values())
    for k in W_dict:
        W_dict[k] = W_dict[k]/factor

    return W_dict


def TAI_cal(W_dict):

    w_norm_np = np.array(list(W_dict.values()))
    #w_norm_np = w_np/np.max(w_np)

    W_positive = w_norm_np[w_norm_np > 0]

    TAI_positive = geometric_mean(W_positive)
    W_norm_switched = W_dict.copy()
    for k, key in enumerate(W_norm_switched.keys()):
        if w_norm_np[k] == 0:
            W_norm_switched[key] = TAI_positive

    return W_norm_switched

# ---------------------------------------------------


class TAI(object):
    def __init__(self,  tgcn):
        """
        :param cds_path: cds file (fasta)
        """
        self.Sij = Sij

        W_dict = weight_cal(self.Sij, tgcn)
        TAI_weight = TAI_cal(W_dict)
        self.index = TAI_weight


# from yeast!!
Sij = {'A:T': 0,
       'G:C': 0,
       'T:A': 0,
       'C:G': 0.41,
       'G:T': 0.28,
       'A:C': 0.9999,
       'A:A': 0.68,
       'T:G': 0.89}

