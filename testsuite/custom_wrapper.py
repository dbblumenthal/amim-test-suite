from testsuite.algorithm_wrapper import AlgorithmWrapper
import networkx as nx
import subprocess
import pandas as pd

# import os
# from testsuite.utils import *
#
# os.chdir('/home/olga/Dropbox/testing-onfah/testsuite')
# expression_data = load_expression_data("GSE3790")
# phenotypes = load_phenotypes("GSE3790")
# ggi_network = load_ggi_network("HPRD", expression_data)
# prefix = 'HPRD_ORIGINAL'

class CustomWrapper(AlgorithmWrapper):

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

        # 1. Write GGI network in format required by your method
        path_to_network = f'../temp/{prefix}_custom_ggi.txt'

        with open(path_to_network, 'w') as edge_list_file:
            gene_ids = nx.get_node_attributes(ggi_network, 'GeneID')
            for u, v in ggi_network.edges():
                edge_list_file.write(f'{gene_ids[u]}\t{gene_ids[v]}\n')

        # 2. Write expression data and phenotype in the format required by your method
        path_to_expression = f'../temp/{prefix}_custom_expr.csv'
        expression_data.T.to_csv(path_to_expression)

        # 3. Insert the command to run your method, direct the output to path_to_output
        path_to_output = f'../temp/{prefix}_custom_output.txt'
        command = f'cd ../../; python bicon_test.py {prefix}'
        subprocess.call(command, shell = True, stdout=subprocess.PIPE)

        # 4. Process results such that they are formatted as a list of strings (entez IDs)
        result_genes = []
        with open(path_to_output, 'r') as results:
            for line in results:
                result_genes.append(line.strip())

        # Delete temporary data.
        subprocess.call(f'rm ../temp/{prefix}_custom_*', shell=True)

        return result_genes, AlgorithmWrapper.mean_degree(ggi_network, result_genes)
