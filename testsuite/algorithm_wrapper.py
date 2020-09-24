import abc
import networkx as nx
import numpy as np


class AlgorithmWrapper(object):
    """Abstract wrapper class for the network enrichment algorithms used in the tests."""

    @abc.abstractmethod
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
        pass

    @staticmethod
    def save_network_as_edge_list(ggi_network, path_to_edge_list, sep, header):
        """Saves a GGI network as an edge list.

        Parameters
        ----------
        ggi_network : nx.Graph
            GGI network that should be saved.
        path_to_edge_list : str
            Path to the output file.
        sep : str
            Separator for source and target of an edge.
        header : str or None
            If not None, a header is written to the output file.
        """
        with open(path_to_edge_list, 'w') as edge_list_file:
            if header is not None:
                edge_list_file.write(f'{header}\n')
            gene_ids = nx.get_node_attributes(ggi_network, 'GeneID')
            for u, v in ggi_network.edges():
                edge_list_file.write(f'{gene_ids[u]}{sep}{gene_ids[v]}\n')

    @staticmethod
    def save_array(array, path_to_array, sep, header):
        """Saves an array in a text-file.

        Parameters
        ----------
        array : list or np.array, shape (n,)
            Array that should be saved.
        path_to_array : str
            Path to the output file.
        sep : str
            Separator for the values in the array.
        header : str or None
            If not None, a header is written to the output file.
        """
        with open(path_to_array, 'w') as array_file:
            if header is not None:
                array_file.write(f'{header}{sep}')
            for value in array:
                array_file.write(f'{value}{sep}')

    @staticmethod
    def mean_degree(ggi_network, result_genes):
        """Computes mean degree of the result genes.

        Parameters
        ----------
        ggi_network : nx.Graph
            GGI network.
        result_genes : str
            Gene IDs of the nodes contained in the result.

        Returns
        -------
        mean_degree : float
            Mean degree of the result genes.
        """
        gene_ids = nx.get_node_attributes(ggi_network, 'GeneID')
        degrees = [ggi_network.degree[node] for node in ggi_network.nodes() if gene_ids[node] in set(result_genes)]
        return np.mean(degrees)

