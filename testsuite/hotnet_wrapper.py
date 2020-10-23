from testsuite.algorithm_wrapper import AlgorithmWrapper
import networkx as nx
import numpy as np
import subprocess
import multiprocessing as mp
import itertools as itt


class HotNetWrapper(AlgorithmWrapper):

    @staticmethod
    def permute_scores(path_scores, path_bins, path_permuted_scores):
        hotnet = 'cd ../algorithms/hierarchical_hotnet/src; python permute_scores.py'
        command = f'{hotnet} -i ../../{path_scores} -bf ../../{path_bins} -o ../../{path_permuted_scores}'
        subprocess.call(command, shell=True, stdout=subprocess.PIPE)

    @staticmethod
    def construct_hierarchy(path_sim_matrix, path_map, path_scores, path_h_ggi, path_h_map):
        hotnet = 'cd ../algorithms/hierarchical_hotnet/src; python construct_hierarchy.py'
        command = f'{hotnet} -smf ../../{path_sim_matrix} -igf ../../{path_map} -gsf ../../{path_scores} -helf ../../{path_h_ggi} -higf ../../{path_h_map}'
        subprocess.call(command, shell=True, stdout=subprocess.PIPE)

    def run_algorithm(self, ggi_network, expression_data, phenotypes, seed_genes, p_values, indicator_matrix, prefix):
        """Runs the algorithm.

        Parameters
        ----------
        ggi_network : nx.Graph
            Possibly randomized GGI network.
        expression_data : pd.DataFrame
            Expression data (indices are sample IDs, column names are gene IDs).
        phenotypes : np.array, shape (n_samples,)
            Phenotype data (indices are sample IDs).
        seed_genes : list of str
            Seed genes (entries are gene IDs).
        p_values : dict of str: float
            P-values for all genes (keys are gene IDs).
        indicator_matrix : pd.DataFrame
            Indicator matrix obtained from expression data (indices are sample IDs, column names are gene IDs).
        prefix : str
            Prefix to be used for temporary files and directories.

        Returns
        -------
        result_genes : list of str
            Set of genes computed by the algorithm.
        mean_degree : float
            Mean degree of the result genes.
        """

        # Write GGI network, mapping, and scores in format required by Hierarchical HotNet.
        path_ggi = f'../temp/{prefix}_hotnet_ggi.tsv'
        path_map = f'../temp/{prefix}_hotnet_map.tsv'
        path_scores = f'../temp/{prefix}_hotnet_scores.tsv'
        with open(path_ggi, 'w') as ggi_file, open(path_map, 'w') as map_file, open(path_scores, 'w') as scores_file:
            first_non_isolated_node = -1
            last_non_isolated_node = -1
            for node in ggi_network.nodes():
                if ggi_network.degree[node] > 0:
                    last_non_isolated_node = node
                    if first_non_isolated_node == -1:
                        first_non_isolated_node = node
            gene_ids = nx.get_node_attributes(ggi_network, 'GeneID')
            for node in range(first_non_isolated_node, last_non_isolated_node + 1):
                gene_index = node - first_non_isolated_node + 1
                gene_id = gene_ids[node]
                map_file.write(f'{gene_index}\t{gene_id}\n')
                scores_file.write(f'{gene_index}\t{-np.log10(p_values[gene_id])}\n')
            for u, v in ggi_network.edges():
                ggi_file.write(f'{u - first_non_isolated_node + 1}\t{v - first_non_isolated_node + 1}\n')

        # path_ggi = f'../temp/network_1_edge_list.tsv'
        # path_map = f'../temp/network_1_index_gene.tsv'
        # path_scores = f'../temp/scores_1.tsv'

        # Compile Fortran module and set number of cores.
        command = 'cd ../algorithms/hierarchical_hotnet/src; f2py -c fortran_module.f95 -m fortran_module > /dev/null'
        subprocess.call(command, shell=True, stdout=subprocess.PIPE)
        num_cores = 34

        # Construct similarity matrix.
        print('Construct similarity matrix.')
        hotnet = 'cd ../algorithms/hierarchical_hotnet/src; python construct_similarity_matrix.py'
        path_sim_matrix = f'../temp/{prefix}_hotnet_sim_matrix.h5'
        command = f'{hotnet} -i ../../{path_ggi} -o ../../{path_sim_matrix}'
        subprocess.call(command, shell=True, stdout=subprocess.PIPE)

        # Find permutation bins.
        hotnet = 'cd ../algorithms/hierarchical_hotnet/src; python find_permutation_bins.py'
        path_bins = f'../temp/{prefix}_hotnet_bins.tsv'
        command = f'{hotnet} -elf ../../{path_ggi} -igf ../../{path_map} -gsf ../../{path_scores} -o ../../{path_bins}'
        subprocess.call(command, shell=True, stdout=subprocess.PIPE)

        # Permute the scores in parallel.
        print('Permute the scores in parallel.')
        pool = mp.Pool(num_cores)
        paths_permuted_scores = [f'../temp/{prefix}_hotnet_scores_{i}.tsv' for i in range(100)]
        args = list(itt.product([path_scores], [path_bins], paths_permuted_scores))
        pool.starmap(self.permute_scores, args)

        # Construct the hierarchies in parallel.
        print('Construct the hierarchies in parallel.')
        path_h_ggi = f'../temp/{prefix}_hotnet_h_ggi.tsv'
        path_h_map = f'../temp/{prefix}_hotnet_h_map.tsv'
        paths_permuted_h_ggis = [f'../temp/{prefix}_hotnet_h_ggi_{i}.tsv' for i in range(100)]
        paths_permuted_h_maps = [f'../temp/{prefix}_hotnet_h_map_{i}.tsv' for i in range(100)]
        args = [(path_sim_matrix, path_map, paths_permuted_scores[i], paths_permuted_h_ggis[i], paths_permuted_h_maps[i])
                for i in range(100)]
        args.append((path_sim_matrix, path_map, path_scores, path_h_ggi, path_h_map))
        pool.starmap(self.construct_hierarchy, args)

        # Process the hierarchies.
        print('Process the hierarchies.')
        hotnet = 'cd ../algorithms/hierarchical_hotnet/src; python process_hierarchies.py'
        path_output = f'../temp/{prefix}_hotnet_results.tsv'
        pelf_str = ' '.join([f'../../{path}' for path in paths_permuted_h_ggis])
        pigf_str = ' '.join([f'../../{path}' for path in paths_permuted_h_maps])
        command = f'{hotnet} -oelf ../../{path_h_ggi} -oigf ../../{path_h_map} -pelf {pelf_str} -pigf {pigf_str} -cf ../../{path_output} -nc {num_cores}'
        subprocess.call(command, shell=True, stdout=subprocess.PIPE)

        # Read the results.
        print('Read the results.')
        result_genes = []
        with open(path_output, 'r') as hotnet_results:
            for line in hotnet_results:
                if not line.startswith('#'):
                    result_genes = line.strip().split('\t')
                    break
        print(result_genes)

        # Delete temporary data.
        print('Delete temporary data.')
        subprocess.call(f'rm ../temp/{prefix}_hotnet_*', shell=True)

        # Return the results.
        return result_genes, AlgorithmWrapper.mean_degree(ggi_network, result_genes)




