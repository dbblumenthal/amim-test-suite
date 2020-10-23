from testsuite.algorithm_wrapper import AlgorithmWrapper
import networkx as nx
import subprocess
import pandas as pd
import numpy as np
import os
import sys
#os.chdir("algorithms/NetCore")
sys.path.append('algorithms/NetCore/netcore/')
from permutations_test import make_network_permutations


class NetCoreWrapper(AlgorithmWrapper):

    def run_algorithm(self, ggi_network, expression_data, phenotypes, seed_genes, gene_scores, indicator_matrix,
                      prefix):
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
        gene_scores : dict of str: float
            Scores for all genes (keys are gene IDs).
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

        # Write GGI network in format required by netcore (no header)
        path_to_network = f'../temp/{prefix}_netcore.txt'
        with open(path_to_network, 'w') as edge_list_file:
            gene_ids = nx.get_node_attributes(ggi_network, 'GeneID')
            for u, v in ggi_network.edges():
                edge_list_file.write(f'{gene_ids[u]}\t{gene_ids[v]}\n')

        path_to_pvals = f'../temp/{prefix}_netcore_pvals.txt'
        gene_scores = [(k,-np.log10(v)) for (k,v) in gene_scores.items()]
        gene_scores = pd.DataFrame(gene_scores)
        gene_scores.to_csv(path_to_pvals,index = False, sep = "\t")

        permutations_num = 100
        make_network_permutations(net_file=path_to_network,
                                  net_name=f'{prefix}_netcore',
                                  output_path="../temp/",
                                  num_perm=permutations_num,
                                  swap_factor=100,
                                  num_cores=34)


        # Run netcore
        edge_list_path = f'../../../temp/{prefix}_netcore.txt'
        weights_path = f'../../../temp/{prefix}_netcore_pvals.txt'
        permutations_path = f'../../../temp/{prefix}_netcore_edge_permutations/'
        results_path = f'../../../temp/netcore_results/'
        command = f'cd ../algorithms/NetCore/netcore; python netcore.py -e {edge_list_path} -w {weights_path} -pd {permutations_path} -o {results_path} -np {permutations_num-1}'
        subprocess.call(command, shell = True)

        # Read the results.
        path_to_output = '../temp/netcore_results/core_norm_subnetworks.txt'
        result_genes = pd.read_csv(path_to_output,sep = ',')
        result_genes = list(result_genes[result_genes.sum_weights == max(result_genes.sum_weights)].nodes)[0].split("\t")

        # Delete temporary data.
        subprocess.call(f'rm -r ../temp/{prefix}_netcore*', shell=True)
        subprocess.call(f'rm -r ../temp/netcore_results', shell=True)


        # Return the results.
        return result_genes, AlgorithmWrapper.mean_degree(ggi_network, result_genes)
