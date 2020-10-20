from testsuite.algorithm_wrapper import AlgorithmWrapper
import networkx as nx
import subprocess
import pandas as pd


class GrandForestWrapper(AlgorithmWrapper):

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

        # Write GGI network in format required by ClustEx2.
        path_to_network = f'../temp/{prefix}_gf_ggi.txt'
        AlgorithmWrapper.save_network_as_edge_list(ggi_network, path_to_network, '\t', 'source\ttarget')

        expression_data["phenotype"] = phenotypes
        expression_data.to_csv(f'../temp/{prefix}_gf_expr.txt', index = False)

        AlgorithmWrapper.save_network_as_edge_list(ggi_network, path_to_network, '\t', 'source\ttarget')

        # Run GF.
        command = f'cd ../algorithms/grand_forest/; Rscript grandforest.R {prefix}'
        subprocess.call(command, shell = True, stdout=subprocess.PIPE)

        # Read the results.
        result_genes = []
        path_to_output = f'../temp/{prefix}_gf_output.txt'
        with open(path_to_output, 'r') as results:
            for line in results:
                result_genes.append(line.strip())

        # Delete temporary data.
        subprocess.call(f'rm ../temp/{prefix}_gf_*', shell=True)

        # Return the results.
        return result_genes, AlgorithmWrapper.mean_degree(ggi_network, result_genes)
