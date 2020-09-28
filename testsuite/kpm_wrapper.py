from testsuite.algorithm_wrapper import AlgorithmWrapper
import networkx as nx
import subprocess


class KPMWrapper(AlgorithmWrapper):

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

        # Write GGI network in format required by KPM.
        path_ggi = '../temp/kpm_ggi.sif'
        AlgorithmWrapper.save_network_as_edge_list(ggi_network, path_ggi, '\tpp\t', None)

        # Write the indicator matrix in the format required by KPM.
        path_indicator_matrix = '../temp/kpm_indicator_martrix.txt'
        indicator_matrix.T.to_csv(path_indicator_matrix, sep='\t', header=False)

        # Run PinnacleZ.
        kpm = 'cd ../algorithms/kpm/; java -jar KPM-5.0.jar -strategy=INES -algo=GREEDY -L1=5 -K=2 -pSingleFile=true ' \
              '-maxSolutions=1'
        output_dir = '../temp/kpm/'
        command = f'{kpm} -matrix1=../{path_indicator_matrix} -graphFile=../{path_ggi} -resultsDir=../{output_dir}'
        subprocess.call(command, shell=True)

        # Read the results.
        # with open(path_output, 'r') as results:
        #    result_genes = [line for line in results if not line.startswith('#')][0].strip().split('\t')[-1].split(' ')

        # Delete temporary data.
        # subprocess.call('rm ../temp/pinnaclez_*', shell=True)

        # Return the results.
        return result_genes, AlgorithmWrapper.mean_degree(ggi_network, result_genes)