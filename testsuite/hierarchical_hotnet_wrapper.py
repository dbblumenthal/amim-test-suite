from testsuite.algorithm_wrapper import AlgorithmWrapper
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
        mean_degree : float
            Mean degree of the result genes.
        """

        # Write GGI network in format required by Hierarchical HotNet.
        path_ggi = '../temp/hotnet_ggi.tsv'
        with open(path_ggi, 'w') as network_file:
            for u, v in ggi_network.edges():
                network_file.write(f'{u + 1}\t{v + 1}\n')

        # Write index to gene file and gene score file in format required by Hierarchical HotNet.
        path_index_gene = '../temp/hotnet_index_gene.tsv'
        path_scores = '../temp/hotnet_gene_scores_0.tsv'
        with open(path_index_gene, 'w') as index_gene_file, open(path_scores, 'w') as gene_score_file:
            gene_ids = nx.get_node_attributes(ggi_network, 'GeneID')
            for node in ggi_network.nodes():
                gene_id = gene_ids[node]
                index_gene_file.write(f'{node + 1}\t{gene_id}\n')
                gene_score_file.write(f'{gene_id}\t{gene_scores[gene_id]}\n')

        # Run Hierarchical HotNet.
        hotnet_sim_matrix = 'cd ../algorithms/hierarchical-hotnet/src/; python construct_similarity_matrix.py'
        hotnet_find_bins = 'cd ../algorithms/hierarchical-hotnet/src/; python find_permutation_bins.py'
        hotnet_permute = 'cd ../algorithms/hierarchical-hotnet/src/; python permute_scores.py'
        hotnet_constr_hierarchy = 'cd ../algorithms/hierarchical-hotnet/src/; python construct_hierarchy.py'
        hotnet_proc_hierarchies = 'cd ../algorithms/hierarchical-hotnet/src/; python process_hierarchies.py'
        path_sim_matrix = '../../../temp/hotnet_similarity_matrix.h5'
        command = f'{hotnet_sim_matrix} -i ../../{path_ggi} -o {path_sim_matrix}'
        subprocess.call(command, shell=True)
        path_bins = '../../../temp/hotnet_score_bins.tsv'
        command = f'{hotnet_find_bins} -gsf ../../{path_scores} -igf ../../{path_index_gene} -elf ../../{path_ggi} -o {path_bins}'
        subprocess.call(command, shell=True)
        num_perms = 100
        paths_scores = [path_scores] + [f'../temp/hotnet_gene_scores_{i + 1}.tsv' for i in range(num_perms)]
        for i in range(num_perms):
            command = f'{hotnet_permute} -i ../../{path_scores} -bf {path_bins} -o ../../{paths_scores[i + 1]}'
            subprocess.call(command, shell=True)
        paths_h_ggis = [f'../../../temp/hotnet_ggi_{i}.tsv' for i in range(num_perms + 1)]
        paths_h_index_genes = [f'../../../temp/hotnet_index_gene_{i}.tsv' for i in range(num_perms + 1)]
        for i in range(num_perms + 1):
            command = f'{hotnet_constr_hierarchy} -smf {path_sim_matrix} -igf ../../{path_index_gene} ' \
                      f'-gsf ../../{paths_scores[i]} -helf {paths_h_ggis[i]} -higf {paths_h_index_genes[i]}'
            subprocess.call(command, shell=True)
        path_clusters = '../temp/hotnet_clusters.tsv'
        command = f'{hotnet_proc_hierarchies} -oelf {paths_h_ggis[0]} -oigf {paths_h_index_genes[0]} ' \
                  f'-pelf {" ".join(paths_h_ggis[1:])} -pigf {" ".join(paths_h_index_genes[1:])} -cf ../../{path_clusters}'
        subprocess.call(command, shell=True)

        # Read the results.
        with open(path_clusters, 'r') as hotnet_results:
            result_genes = [line for line in hotnet_results if not line.startswith('#')][0].strip().split('\t')

        # Return the results.
        return result_genes, AlgorithmWrapper.mean_degree(ggi_network, result_genes)
