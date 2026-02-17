from modules.sequence_family.clustering_optimization import *


class SequenceFamilyModule(object):
    @staticmethod
    def run_module(user_input: models.ModuleInput):
        """
        in this part, module input is split into different inputs according to the sequence family theory
        """
        object_list = []
        opt_org_list = [org.name for org in user_input.organisms if org.is_optimized]

        if user_input.clusters_count > len(opt_org_list):
            raise ValueError('Number of clusters must be smaller (or equal) to number of optimized organisms')

        elif user_input.clusters_count>1:
            clustering_mat, opt_org_list = dict_to_cluster_np_array(user_input)
            #todo: decide minimal val to devide into 2 sequences
            best_clustering = create_n_clusters(clustering_mat, n_clus=user_input.clusters_count)
            object_list = return_list_of_sub_microbiomes(best_clustering, user_input)

        elif user_input.clusters_count == 1:
            object_list = [user_input]

        else:
            raise ValueError('Number of clusters should be a positive integer larger than 1 ')

        return object_list
