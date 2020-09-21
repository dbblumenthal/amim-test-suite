from testsuite.algorithm_wrapper import AlgorithmWrapper
import testsuite.utils as utils
import networkx as nx
import subprocess


class DIAMOnDWrapper(AlgorithmWrapper):

    def run_algorithm(self, ggi_network, expression_data, phenotypes, seed_genes):
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

        Returns
        -------
        result_genes : list of str
            Set of genes computed by the algorithm.
        """

        # Write GGI network in format required by DIAMOnD.
        path_to_network = '../temp/DIAMOnD_ggi.txt'
        with open(path_to_network, 'w') as diamond_network:
            gene_ids = nx.get_node_attributes(ggi_network, utils.gene_id_attribute_name())
            for u, v in ggi_network.edges():
                diamond_network.write(f'{gene_ids[u]},{gene_ids[v]}\n')

        # Write seed genes in format required by DIAMOnD.
        path_to_seeds = '../temp/DIAMOnD_seed_genes.txt'
        with open(path_to_seeds, 'w') as diamond_seed_genes:
            for seed_gene in seed_genes:
                diamond_seed_genes.write(f'{seed_gene}\n')

        # Run DIAMOnD.
        path_to_diamond = '../algorithms/DIAMOnD/DIAMOnD.py'
        path_to_output = '../temp/DIAMOnD_results.txt'
        subprocess.call(f'python {path_to_diamond} {path_to_network} {path_to_seeds} 200 1 {path_to_output}')

        # Read the results.
        result_genes = []
        with open(path_to_output, 'r') as diamond_results:
            for line in diamond_results:
                result_genes.append(line.strip())

        # Delete temporary data.
        subprocess.call('rm ../temp/DIAMOnD_*')

        # Return the results.
        return result_genes

