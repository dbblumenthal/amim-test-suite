from testsuite.algorithm_wrapper import AlgorithmWrapper
import testsuite.utils as utils
import networkx as nx
import subprocess


class HierarchicalHotNetWrapper(AlgorithmWrapper):

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

        # Write GGI network in format required by Hierarchical HotNet.
        path_ggi = '../temp/hotnet_ggi.tsv'
        with open(path_ggi, 'w') as network_file:
            for u, v in ggi_network.edges():
                network_file.write(f'{u + 1}\t{v + 1}\n')

        # Write index to gene file and gene score file in format required by Hierarchical HotNet.
        path_index_gene = '../temp/hotnet_index_gene.tsv'
        path_scores = '../temp/hotnet_gene_scores.tsv'
        with open(path_index_gene, 'w') as index_gene_file, open(path_scores, 'w') as gene_score_file:
            gene_ids = nx.get_node_attributes(ggi_network, utils.gene_id_attribute_name())
            for node in ggi_network.nodes():
                gene_id = gene_ids[node]
                index_gene_file.write(f'{node + 1}\t{gene_id}\n')
                gene_score_file.write(f'{gene_id}\t{gene_scores[gene_id]}\n')

        # Run Hierarchical HotNet.
        hotnet_sim_matrix = '../algorithms/hierarchical-hotnet/src/construct_similarity_matrix.py'
        hotnet_find_bins = '../algorithms/hierarchical-hotnet/src/find_permutation_bins.py'
        hotnet_permute = '../algorithms/hierarchical-hotnet/src/permute_scores.py'
        path_sim_matrix = '../temp/hotnet_similarity_matrix.h5'
        command = f'python {hotnet_sim_matrix} -i {path_ggi} -o {path_sim_matrix}'
        subprocess.call(command, shell=True)
        path_bins = '../temp/hotnet_score_bins.tsv'
        command = f'python {hotnet_find_bins} -gsf {path_scores} -igf {path_index_gene} -elf {path_ggi} -o {path_bins}'
        subprocess.call(command, shell=True)
        num_perms = 100
        paths_scores = [path_scores] + [f'../temp/hotnet_gene_scores_{i}.tsv' for i in range(num_perms)]
        for i in range(num_perms):
            command = f'python {hotnet_permute} -i {path_scores} -bf {path_bins} -o {path_scores[i + 1]}'
            # todo: continue here

        result_genes = []

        # Return the results.
        return result_genes