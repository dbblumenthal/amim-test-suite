from testsuite.algorithm_wrapper import AlgorithmWrapper
import networkx as nx
import subprocess


class GXNAWrapper(AlgorithmWrapper):

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

        # Write GGI network in format required by GXNA.
        path_ggi = '../temp/gxna_ggi.txt'
        AlgorithmWrapper.save_network_as_edge_list(ggi_network, path_ggi, ' ', None)

        # Write the expression data in the format required by GXNA.
        path_expr = '../temp/gxna_expression.txt'
        gene_ids = list(expression_data.columns)
        with open(path_expr, 'w') as expression_data_file:
            for gene_id in gene_ids:
                expression_data_file.write(gene_id)
                for sample_id in range(expression_data.shape[0]):
                    expression_data_file.write(f' {expression_data.loc[sample_id, gene_id]}')
                expression_data_file.write('\n')

        # Write the phenotypes in the format required by GXNA.
        path_phen = '../temp/gxna_phenotypes.txt'
        AlgorithmWrapper.save_array(phenotypes, path_phen, ' ', None)

        # Write the mapping file in the format required by GXNA.
        path_map = '../temp/gxna_mapping.txt'
        with open(path_map, 'w') as mapping_file:
            for gene_id in gene_ids:
                mapping_file.write(f'{gene_id} {gene_id}\n')

        # Run GXNA.
        gxna = 'cd ../algorithms/gxna/; ./gxna'
        command = f'{gxna} -name gxna -phenoFile ../{path_phen} -expFile ../{path_expr} -edgeFile ../{path_ggi} -mapFile ../{path_map}'
        subprocess.call(command, shell=True)

        # Read the results.
        result_genes = []
        path_to_output = '../temp/gxna_000_0.txt'
        with open(path_to_output, 'r') as gxna_results:
            for line in gxna_results:
                result_genes.append(line.strip().split()[1])

        # Delete temporary data.
        subprocess.call('rm ../temp/gxna_*', shell=True)

        # Return the results.
        return result_genes, AlgorithmWrapper.mean_degree(ggi_network, result_genes)
