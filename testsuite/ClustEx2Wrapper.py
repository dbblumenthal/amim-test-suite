from testsuite.algorithm_wrapper import AlgorithmWrapper
import testsuite.utils as utils
import networkx as nx
import subprocess


class ClustEx2Wrapper(AlgorithmWrapper):

    def run_algorithm(self, ggi_network, expression_data, phenotypes, seed_genes, gene_scores):
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

        Returns
        -------
        result_genes : list of str
            Set of genes computed by the algorithm.
        """

        # Write GGI network in format required by ClustEx2.
        path_to_network = '../temp/clustex2_ggi.txt'
        utils.save_network_as_edge_list(ggi_network, path_to_network, '\t', 'source\ttarget')

        # Write seed genes in format required by DIAMOnD.
        path_to_seeds = '../temp/clustex2_seed_genes.txt'
        utils.save_array(seed_genes, path_to_seeds, '\n', 'gene')

        # Run ClustEx2.
        clustex2 = '../algorithms/clustex2/clustex2'
        command = f'{clustex2} --gene_list {path_to_seeds} --network {path_to_network} -D -G -C -s 100 -j clustex2'
        subprocess.call(command, shell=True)

        # Read the results.
        result_genes = []
        path_to_output = '../temp/clustex2_genes_100_0.1.txt'
        with open(path_to_output, 'r') as clustex2_results:
            for line in clustex2_results:
                result_genes.append(line.strip())

        # Delete temporary data.
        subprocess.call('rm ../temp/clustex2_*')
        subprocess.call('rm ../algorithms/clustex2/clustex2_*')

        # Return the results.
        return result_genes
