from testsuite.algorithm_wrapper import AlgorithmWrapper
import networkx as nx
import subprocess
import pandas as pd


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

        # 2. Write expression data in format required by your method
        path_to_expression = f'../temp/{prefix}_custom_expr.txt'

        # 3. Write phenotype data in format required by your method
        path_to_expression = f'../temp/{prefix}_custom_phenotype.txt'

        # 4. Insert the command to run your method, direct the output to path_to_output
        path_to_output = f'../temp/{prefix}_custom_output.txt'
        command = ""
        subprocess.call(command, shell=True, stdout=subprocess.PIPE)

        # 5. Process results such that they are formatted as a list of strings (entez IDs)
        result_genes = []

        # 6. Delete temporary data.
        subprocess.call(f'rm ../temp/{prefix}_custom_*', shell=True)

        return result_genes, AlgorithmWrapper.mean_degree(ggi_network, result_genes)