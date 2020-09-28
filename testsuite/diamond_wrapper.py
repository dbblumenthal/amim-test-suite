from testsuite.algorithm_wrapper import AlgorithmWrapper
import networkx as nx
import subprocess


class DIAMOnDWrapper(AlgorithmWrapper):

    def run_algorithm(self, ggi_network, expression_data, phenotypes, seed_genes, gene_scores, indicator_matrix):
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

        Returns
        -------
        result_genes : list of str
            Set of genes computed by the algorithm.
        mean_degree : float
            Mean degree of the result genes.
        """

        # Write GGI network in format required by DIAMOnD.
        path_ggi = '../temp/diamond_ggi.txt'
        AlgorithmWrapper.save_network_as_edge_list(ggi_network, path_ggi, ',', None)

        # Write seed genes in format required by DIAMOnD.
        path_seeds = '../temp/diamond_seed_genes.txt'
        AlgorithmWrapper.save_array(seed_genes, path_seeds, '\n', None)

        # Run DIAMOnD.
        diamond = 'cd ../algorithms/diamond/; python DIAMOnD.py'
        path_output = '../temp/diamond_results.txt'
        command = f'{diamond} ../{path_ggi} ../{path_seeds} 200 1 ../{path_output}'
        subprocess.call(command, shell=True)

        # Read the results.
        result_genes = []
        with open(path_output, 'r') as diamond_results:
            for line in diamond_results:
                result_genes.append(line.strip())

        # Delete temporary data.
        subprocess.call('rm ../temp/diamond_*', shell=True)

        # Return the results.
        return result_genes, AlgorithmWrapper.mean_degree(ggi_network, result_genes)

