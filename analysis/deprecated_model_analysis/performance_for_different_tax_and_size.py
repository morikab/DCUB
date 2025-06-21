from tax_splitting_func import format_output
import json
from Bio.SeqIO import read
from modules.main import run_orf_optimization
import matplotlib.pyplot as  plt
import time
import numpy as np

# open organism dict
with open('data_for_analysis/org_name_to_dict.json') as f:
    org_dict = json.load(f)

gene = read('zorA anti-phage defense.fasta', 'fasta')
initial_seq = str(gene.seq)
runtimes_for_tax = 10
translation_function = 'zscore_hill_climbing_average'


model_preferences = {'RE': False, #todo: test restcition enzymes
                     'translation': True,
                     'transcription': False,
                     'translation_function': translation_function#, 'single_codon_global', 'single_codon_local’, 'zscore_hill_climbing_average', 'zscore_hill_climbing_weakest_link'
}

tax_to_color = {'Phylum':'red', 'Family':'blue', 'Genus':'green'}
all_data = {}
for tax_con, color in tax_to_color.items():
    print(tax_con, color)
    np_runs_for_tax = np.empty(shape=(len(org_dict)-2, runtimes_for_tax))
    for i in range(runtimes_for_tax):
        tic = time.time()
        scores_for_x = []
        for x in range(2,len(org_dict)):
            inner_tic = time.time()
            usr_dict = format_output(initial_seq, org_dict, tax_con, x, run_times=50)
            final_cds, optimization_index, weakest_score = run_orf_optimization(usr_dict, model_preferences)
            scores_for_x.append(optimization_index)
            print('TIME; ', time.time()-inner_tic)
        np_runs_for_tax[:,i]=scores_for_x
    opt_scores = np.nanmean(np_runs_for_tax, axis=1)
    sub_microbiome_size = np.array(list(range(2,len(org_dict))))
    plt.plot(sub_microbiome_size, opt_scores, color=color)
    all_data[tax_con] = sub_microbiome_size.tolist()


translation_function_name = translation_function.replace('_', ' ')
plt.title(f'Scale-up analysis for {translation_function_name} optimization of translation')
plt.xlabel('Size of sub-microbiome')
plt.ylabel('Optimization score')
plt.legend(["Phylum", "Family", "Genus"], loc ="upper right")
plt.show()

with open(f'{translation_function} optimization per tax.json', 'w') as fp:
    json.dump(all_data, fp)