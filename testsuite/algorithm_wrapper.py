import abc


class AlgorithmWrapper(object):
    """Abstract wrapper class for the network enrichment algorithms used in the tests."""

    @abc.abstractmethod
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
        pass
