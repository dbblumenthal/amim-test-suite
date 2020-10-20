from testsuite.algorithm_wrapper import AlgorithmWrapper
import networkx as nx
import subprocess


class PinnacleZWrapper(AlgorithmWrapper):

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

        # Write GGI network in format required by PinnacleZ.
        path_ggi = f'../temp/{prefix}_pinnaclez_ggi.sif'
        AlgorithmWrapper.save_network_as_edge_list(ggi_network, path_ggi, '\tpp\t', None)

        # Write the expression data in the format required by PinnacleZ.
        path_expr = f'../temp/{prefix}_pinnaclez_expression.txt'
        expression_data.T.to_csv(path_expr, sep='\t', index_label='genes')

        # Write the phenotypes in the format required by PinnacleZ.
        path_phen = f'../temp/{prefix}_pinnaclez_phenotypes.txt'
        AlgorithmWrapper.save_array(phenotypes, path_phen, '\n', None, True)

        # Run PinnacleZ.
        pinnaclez = 'cd ../algorithms/pinnaclez/; java -Xmx2g -jar pinnaclez.jar'
        path_output = f'../temp/{prefix}_pinnaclez_results.txt'
        command = f'{pinnaclez} ../{path_phen} ../{path_expr} ../{path_ggi} -o ../{path_output}'
        subprocess.call(command, shell=True)

        # Read the results.
        with open(path_output, 'r') as results:
            result_genes = [line for line in results if not line.startswith('#')][0].strip().split('\t')[-1].split(' ')

        # Delete temporary data.
        subprocess.call(f'rm ../temp/{prefix}_pinnaclez_*', shell=True)

        # Return the results.
        return result_genes, AlgorithmWrapper.mean_degree(ggi_network, result_genes)