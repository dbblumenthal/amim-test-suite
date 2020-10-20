from testsuite.algorithm_wrapper import AlgorithmWrapper
import networkx as nx
import subprocess


class ClustEx2Wrapper(AlgorithmWrapper):

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
        path_to_network = f'../temp/{prefix}_clustex2_ggi.txt'
        AlgorithmWrapper.save_network_as_edge_list(ggi_network, path_to_network, '\t', 'source\ttarget')

        # Write seed genes in format required by ClustEx2.
        path_to_seeds = f'../temp/{prefix}_clustex2_seed_genes.txt'
        AlgorithmWrapper.save_array(seed_genes, path_to_seeds, '\n', 'gene')

        # Run ClustEx2.
        clustex2 = 'cd ../algorithms/clustex2/; ./clustex2'
        command = f'{clustex2} --gene_list ../{path_to_seeds} --network ../{path_to_network} -D -G -C -s 100 -j {prefix}_clustex2'
        subprocess.call(command, shell=True, stdout=subprocess.PIPE)

        # Read the results.
        result_genes = []
        path_to_output = f'../temp/{prefix}_clustex2_genes_100_0.1.txt'
        with open(path_to_output, 'r') as clustex2_results:
            for line in clustex2_results:
                result_genes.append(line.strip())

        # Delete temporary data.
        subprocess.call(f'rm ../temp/{prefix}_clustex2_*', shell=True)
        subprocess.call(f'rm ../algorithms/clustex2/{prefix}_clustex2_*', shell=True)

        # Return the results.
        return result_genes, AlgorithmWrapper.mean_degree(ggi_network, result_genes)
