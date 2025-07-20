from scipy.stats import spearmanr
import numpy as np
from modules import models
from sklearn.cluster import AgglomerativeClustering, KMeans


# TODO: try different clustering methods, and cluster eval methods- DBI, Silhouette- find min val
# TODO: write efficiently

def dict_to_cluster_np_array(user_input: models.ModuleInput): #todo: make this work on tai as well
    clustering_mat = []
    opt_org_list = []
    for organism in user_input.organisms:
        if organism.is_optimized:
            cai_profile = organism.cai_profile
            clustering_mat.append(list(cai_profile.values()))
            opt_org_list.append(organism.name)
    clustering_mat=np.array(clustering_mat)
    return clustering_mat, opt_org_list


def make_distance_matrix(clustering_mat):
    n_samples = clustering_mat.shape
    n_samples = n_samples[0]
    distance_matrix = np.zeros(shape=(n_samples, n_samples))
    for i in range(n_samples):
        sample_i = clustering_mat[i,:]
        for k in range(n_samples):
            sample_k = clustering_mat[k,:]
            distance_matrix[i,k] = 1-spearmanr(sample_k, sample_i)[0] # using -1 for distance instead of similarity
    return distance_matrix


def create_n_clusters(clustering_mat, n_clus):
    dist_metric = 'precomputed'
    distance_mat = make_distance_matrix(clustering_mat)
    clustering = AgglomerativeClustering(n_clusters=n_clus,
                                         affinity=dist_metric,
                                         linkage='average').fit(distance_mat, )

    return clustering.labels_
#
# def find_best_clustering(clustering_mat, max_clus_num, c_index='dbi', c_method = 'alggomerative' ):
#     dist_metric = 'precomputed'
#     scores = []
#     clusters = []
#
#     distance_mat = make_distance_matrix(clustering_mat)
#     n_samples = distance_mat.shape
#     n_samples = n_samples[0]
#     for n_clus in range(2, min(n_samples, max_clus_num)):
#         ##### clustering options ######
#         # if c_method == 'kmeans':
#         #     clustering = KMeans(n_clusters=n_clus).fit(clustering_mat)
#         # else:
#         clustering = AgglomerativeClustering(n_clusters=n_clus,
#                                              affinity= dist_metric,
#                                              linkage='average').fit(distance_mat, )
#
#         labels = clustering.labels_
#
#         ##### cluster eval indexes ####
#         if c_index == 'dbi':
#             score = davies_bouldin_score(clustering_mat, labels)
#         else:
#             score = silhouette_score(clustering_mat, labels, metric= dist_metric)
#
#         scores.append(score)
#         clusters.append(labels)
#
#     best_score = min(scores)
#     best_clusturing = clusters[scores.index(best_score)]
#     return best_clusturing, best_score


def return_list_of_sub_microbiomes(best_clusturing: list, module_input: models.ModuleInput):
    opt_org_list = [org.name for org in module_input.organisms if org.is_optimized]
    deopt_org_list = [org.name for org in module_input.organisms if not org.is_optimized]
    c_assignment_dict = dict(zip(opt_org_list, best_clusturing))
    opt_org_clusters = list({n: [k for k in c_assignment_dict.keys() if c_assignment_dict[k] == n]
         for n in set(c_assignment_dict.values())}.values())

    inp_obj_list = []
    for c_opt_org_list in opt_org_clusters:
        opt_and_deopt = c_opt_org_list+deopt_org_list
        new_module_input = models.ModuleInput(
            organisms=[],
            sequence=module_input.sequence,
            tuning_parameter=module_input.tuning_parameter,
            clusters_count=module_input.clusters_count,
        )

        new_module_input.organisms = [module_input.organisms[i] for i in range(len(opt_org_list + deopt_org_list))
                                    if module_input.organisms[i].name in opt_and_deopt]
        inp_obj_list.append(new_module_input)
    return inp_obj_list
