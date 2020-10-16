from testsuite.algorithm_wrapper import AlgorithmWrapper
import networkx as nx
import subprocess


class GiGAWrapper(AlgorithmWrapper):

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

        # Write GGI network in format required by GiGA.
        path_ggi = f'../temp/{prefix}_giga_ggi.txt'
        AlgorithmWrapper.save_network_as_edge_list(ggi_network, path_ggi, '\t', None)

        # Write the sorted gene list in the format required by GiGA.
        path_gene_list = f'../temp/{prefix}_giga_sorted_genes.txt'
        sorted_gene_list = sorted(list(p_values.items()), key=lambda item: item[1])
        with open(path_gene_list, 'w') as gene_list:
            for gene, _ in sorted_gene_list:
                gene_list.write(f'{gene}\n')

        # Write the mapping file in the format required by GiGA.
        path_map = f'../temp/{prefix}_giga_mapping.txt'
        gene_ids = list(expression_data.columns)
        with open(path_map, 'w') as mapping_file:
            for gene_id in gene_ids:
                mapping_file.write(f'{gene_id}\t{gene_id}\n')

        # Run GiGA.
        giga = 'cd ../algorithms/giga/; perl GiGA.pl'
        path_output = f'../temp/{prefix}_giga_results.txt'
        command = f'{giga} -i../{path_gene_list} -n../{path_ggi} -g../{path_map} -o../{path_output} -fTXT'
        subprocess.call(command, shell=True, stdout=subprocess.PIPE)

        # Read the results.
        result_genes = []
        with open(path_output, 'r') as giga_results:
            index = 0
            for line in giga_results:
                index += 1
                if index == 1:
                    continue
                if line.startswith('-'):
                    result_genes.append(line.strip().split('\t')[1])
                else:
                    break

        # Delete temporary data.
        subprocess.call(f'rm ../temp/{prefix}_giga_*', shell=True)

        # Return the results.
        return result_genes, AlgorithmWrapper.mean_degree(ggi_network, result_genes)
