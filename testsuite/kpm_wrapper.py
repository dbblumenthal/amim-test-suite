from testsuite.algorithm_wrapper import AlgorithmWrapper
import networkx as nx
import subprocess
import numpy as np
import os


class KPMWrapper(AlgorithmWrapper):

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

        # Write GGI network in format required by KPM.
        path_ggi = f'../temp/{prefix}_kpm_ggi.sif'
        AlgorithmWrapper.save_network_as_edge_list(ggi_network, path_ggi, '\tpp\t', None)

        # Write the indicator matrix in the format required by KPM.
        path_indicator_matrix = f'../temp/{prefix}_kpm_indicator_martrix.txt'
        indicator_matrix.T.to_csv(path_indicator_matrix, sep='\t', header=False)


        # Run KPM.
        num_case_exceptions = int(np.ceil(indicator_matrix.shape[0] / 10))
        kpm = f'cd ../algorithms/kpm/; java -Xmx2g -jar keypathwayminer-standalone-5.0.jar -strategy=INES -algo=GREEDY -L1={num_case_exceptions} -K=2 -maxSolutions=1'
        output_dir = f'../temp/{prefix}_kpm/'
        command = f'{kpm} -matrix1=../{path_indicator_matrix} -graphFile=../{path_ggi} -resultsDir=../{output_dir}'
        subprocess.call(command, shell=True)

        # Read the results.
        path_output = None
        for directory in os.listdir(f'../temp/{prefix}_kpm/tables/'):
            if directory.startswith('INES_GREEDY'):
                path_output = f'../temp/{prefix}_kpm/tables/{directory}/Pathway-k-2-l-{num_case_exceptions}-NODES-KPM.txt'
                break
        with open(path_output, 'r') as results:
            result_genes = [line.strip().split('\t')[1] for line in results]

        # Delete temporary data.
        subprocess.call(f'rm -rf ../temp/{prefix}_kpm/', shell=True)
        subprocess.call(f'rm ../temp/{prefix}_kpm_*', shell=True)

        # Return the results.
        return result_genes, AlgorithmWrapper.mean_degree(ggi_network, result_genes)
