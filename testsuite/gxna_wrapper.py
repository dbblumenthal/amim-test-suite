from testsuite.algorithm_wrapper import AlgorithmWrapper
import networkx as nx
import subprocess


class GXNAWrapper(AlgorithmWrapper):

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

        # Write GGI network in format required by GXNA.
        path_ggi = f'../temp/{prefix}_gxna_ggi.txt'
        AlgorithmWrapper.save_network_as_edge_list(ggi_network, path_ggi, ' ', None)

        # Write the expression data in the format required by GXNA.
        path_expr = f'../temp/{prefix}_gxna_expression.txt'
        expression_data.to_csv(path_expr, sep=' ', header=False)

        # Write the phenotypes in the format required by GXNA.
        path_phen = f'../temp/{prefix}_gxna_phenotypes.txt'
        gene_ids = list(expression_data.columns)
        AlgorithmWrapper.save_array(phenotypes, path_phen, ' ', None)

        # Write the mapping file in the format required by GXNA.
        path_map = f'../temp/{prefix}_gxna_mapping.txt'
        with open(path_map, 'w') as mapping_file:
            for gene_id in gene_ids:
                mapping_file.write(f'{gene_id} {gene_id}\n')

        # Run GXNA.
        gxna = 'cd ../algorithms/gxna/; ./gxna'
        command = f'{gxna} -name {prefix}_gxna -phenoFile ../{path_phen} -expFile ../{path_expr} -edgeFile ../{path_ggi} -mapFile ../{path_map}'
        subprocess.call(command, shell=True, stdout=subprocess.PIPE)

        # Read the results.
        result_genes = []
        path_to_output = f'../temp/{prefix}_gxna_000_0.txt'
        with open(path_to_output, 'r') as gxna_results:
            for line in gxna_results:
                result_genes.append(line.strip().split()[1])

        # Delete temporary data.
        subprocess.call(f'rm ../temp/{prefix}_gxna_*', shell=True)

        # Return the results.
        return result_genes, AlgorithmWrapper.mean_degree(ggi_network, result_genes)
