from not_used_code.running_modules_functions import *
import json
from Bio.SeqIO import read
from Bio import pairwise2
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import time
import numpy as np
import pandas as pd
from modules.main import run_orf_optimization



# ribosomal_msa_dict = {}
# for seq in parse('16s_sequences_msa_clustalo.fasta', 'fasta'):
#     ribosomal_msa_dict[seq.id.replace('_', ' ')] = str(seq.seq)
# with open('data_for_analysis/org_name_to_16s_msa.json', 'w') as fp:
#     json.dump(ribosomal_msa_dict, fp)

def diff_letters(a,b):
    return sum ( a[i] != b[i] for i in range(len(a)) )


def intersect(lst1, lst2):
    return set(lst1).intersection(lst2)

with open('data_for_analysis/org_name_to_dict.json', 'r') as fp:
    org_dict = json.load(fp)



with open('data_for_analysis/org_name_to_16s_msa.json', 'r') as fp:
    ribosomal_msa_dict = json.load(fp)

with open('data_for_analysis/org_name_to_16s.json', 'r') as fp:
    ribosomal_dict = json.load(fp)



not_really_small_genome = ['Agromyces allii', 'Arthrobacter crystallopoietes', 'Arthrobacter luteolus', 'Arthrobacter pascens',
              'Arthrobacter subterraneus', 'Arthrobacter tumbae', 'Brevibacterium frigoritolerans', 'Janibacter limosus',
              'Knoellia subterranea', 'Mycolicibacterium smegmatis', 'Nocardioides daejeonensis',
              'Nocardioides oleivorans', 'Nocardioides sediminis', 'Nocardioides terrigena',
              'Paenarthrobacter nitroguajacolicus', 'Paenibacillus aceris',
              'Paenibacillus oryzisoli', 'Paenibacillus prosopidis', 'Paenibacillus qinlingensis', 'Pedococcus badiiscoriae',
              'Pedococcus bigeumensis', 'Pedococcus dokdonensis', 'Peribacillus muralis', 'Peribacillus simplex',
              'Phycicoccus duodecadis', 'Priestia flexa', 'Pseudarthrobacter phenanthrenivorans',
              'Rhodanobacter denitrificans', 'Terrabacter aerolatus', 'Terrabacter tumescens',
              'Yonghaparkia alkaliphila'] #deleted all 3 arabidopsis_microbiome that are smaller than 1kb



contain_more_than_10_ribosomal_genes = ['Arthrobacter crystallopoietes', 'Arthrobacter luteolus',
              'Arthrobacter subterraneus', 'Arthrobacter tumbae', 'Brevibacterium frigoritolerans', 'Janibacter limosus',
              'Knoellia subterranea', 'Mycolicibacterium smegmatis', 'Nocardioides daejeonensis',
              'Nocardioides oleivorans', 'Nocardioides sediminis',
              'Paenarthrobacter nitroguajacolicus', 'Paenibacillus aceris',
              'Paenibacillus qinlingensis', 'Pedococcus badiiscoriae',
              'Pedococcus bigeumensis', 'Peribacillus muralis', 'Peribacillus simplex',
              'Phycicoccus duodecadis', 'Priestia flexa',
              'Rhodanobacter denitrificans', 'Terrabacter tumescens']



gene = read('zorA anti-phage defense.fasta', 'fasta')
cds = str(gene.seq)
print(len(cds)/3)
# aln_scores = []

final_tested_org = not_really_small_genome
# func_options = ['single_codon_global', 'single_codon_local', 'zscore_hill_climbing_average', 'zscore_hill_climbing_weakest_link']
func_options = ['single_codon_global', 'zscore_hill_climbing_average']

spearman_dict = {}
for translation_function in func_options:
    opt_scores = []
    msa_scores = []

    tic = time.time()
    for org1 in final_tested_org:
        for org2 in final_tested_org:
            # if org1.split(' ')[0] == org2.split(' ')[0]:
            #     continue
            if org1 == org2:
                continue
            org_dict[org1]['optimized'] = True
            org_dict[org1]['tai_profile'] = {}
            org_dict[org1]['tai_std'] = {}
            org_dict[org1]['tai_avg'] = {}
            org_dict[org2]['optimized'] = False
            org_dict[org2]['tai_profile'] = {}
            org_dict[org2]['tai_std'] = {}
            org_dict[org2]['tai_avg'] = {}
            software_dict = {
                'sequence': cds,
                'tuning_param': 0.5,
                'organisms': {},
                }
            software_dict['organisms'][org1] = org_dict[org1]
            software_dict['organisms'][org2] = org_dict[org2]
            inner_tic = time.time()
            final_cds, optimization_index, weakest_score = run_orf_optimization(software_dict)
            print('TIME: ', time.time()-inner_tic)
            # alignment_score = pairwise2.align.globalxx(
            #     ribosomal_dict[org1], ribosomal_dict[org2], score_only=True)
            # aln_scores.append(alignment_score)
            msa_score = diff_letters(ribosomal_msa_dict[org1], ribosomal_msa_dict[org2])
            msa_scores.append(msa_score)
            opt_scores.append(optimization_index)
    spearman_dict[translation_function] = spearmanr(msa_scores, opt_scores, nan_policy='omit')
    print(spearmanr(msa_scores, opt_scores, nan_policy='omit'))
    trans_func_name = translation_function.replace('_', ' ')
    plt.scatter(msa_scores, opt_scores, s=0.1)

            # print(org1, org2, msa_scores, optimization_index)
    toc = time.time()
    print(toc-tic)

print(spearman_dict)
plt.legend(['single codon optimization', 'hill climbing optimization'], loc ="upper right")
plt.title(f'Evolutionary distance and model performance')
plt.xlabel('# Different aligned characters')
plt.ylabel('Optimization score')
plt.show()


