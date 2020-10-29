import testsuite.utils as utils
import testsuite.meaningfulness_scores as scores
import testsuite.network_generators as generators
import pandas as pd
import numpy as np


class TestRunner(object):
    """Runs the tests."""

    def __init__(self):
        """Constructs TestRunner object."""
        self.condition_selectors = list(utils.ConditionSelector)
        self.algorithm_selectors = list(utils.AlgorithmSelector)
        self.phenotypes = {sel: utils.load_phenotypes(sel) for sel in self.condition_selectors}
        self.pathways = {sel: utils.get_pathways(sel) for sel in self.condition_selectors}
        self.expression_data = {sel: utils.load_expression_data(sel) for sel in self.condition_selectors}
        self.indicator_matrix = {sel: utils.compute_indicator_matrix(self.expression_data[sel], self.phenotypes[sel]) for sel in self.condition_selectors}
        self.gene_p_values = {sel: utils.compute_gene_p_values(self.expression_data[sel], self.phenotypes[sel]) for sel in self.condition_selectors}
        self.seed_genes = {sel: utils.extract_seed_genes(self.gene_p_values[sel]) for sel in self.condition_selectors}
        self.algorithm_wrappers = {sel: utils.get_algorithm_wrapper(sel) for sel in self.algorithm_selectors}
        self.network_generator_names = []
        self.ggi_network_names = []
        self.condition_names = []
        self.algorithm_names = []
        self.random_seeds = []
        self.nums_genes_seeds = []
        self.lcc_ratios_seeds = []
        self.mean_shortest_distances_seeds = []
        self.mean_degrees_results = []
        self.nums_genes_results = []
        self.result_genes = []
        self.mean_mutual_informations = []
        self.neg_log_gsea_p_values = []
        self.results = None
        self.outfile = ''

    def run_on_network(self, ggi_network, seed, ggi_network_name, network_generator_name, condition_selector,
                       algorithm_selector, verbose):
        """Runs the tests for a given condition on a given network.

        Parameters
        ----------
        ggi_network : nx.Graph
            Possibly randomized GGI network.
        seed : convertible to np.uint32 or None
            Seed used by the generator.
        ggi_network_name : str
            String representation of the original GGI network.
        network_generator_name : str
            String representation of the employed network generator.
        condition_selector : ConditionSelector
            Specifies for which condition the tests should be run.
        algorithm_selector : AlgorithmSelector
            Specifies the algorithm that should be run.
        verbose : bool
            Print progress to stdout.
        """
        phenotypes = self.phenotypes[condition_selector]
        pathways = self.pathways[condition_selector]
        expression_data = self.expression_data[condition_selector]
        seed_genes = self.seed_genes[condition_selector]
        p_values = self.gene_p_values[condition_selector]
        indicator_matrix = self.indicator_matrix[condition_selector]
        lcc_ratio, mean_shortest_distance = utils.compute_seed_statistics(ggi_network, seed_genes)
        prefix = f'{ggi_network_name}_{network_generator_name}'
        if verbose:
            print(f'\t\talgorithm = {str(algorithm_selector)}')
        algorithm_wrapper = self.algorithm_wrappers[algorithm_selector]
        result_genes, mean_degree = algorithm_wrapper.run_algorithm(ggi_network, expression_data, phenotypes,
                                                                    seed_genes, p_values, indicator_matrix,
                                                                    prefix)
        mean_mutual_information = scores.compute_mean_mutual_information(expression_data, phenotypes, result_genes)
        neg_log_gsea_p_value = scores.compute_neg_log_gsea_p_value(pathways, result_genes)
        self.ggi_network_names.append(ggi_network_name)
        self.random_seeds.append(seed)
        self.network_generator_names.append(network_generator_name)
        self.condition_names.append(str(condition_selector))
        self.nums_genes_seeds.append(len(seed_genes))
        self.lcc_ratios_seeds.append(lcc_ratio)
        self.mean_shortest_distances_seeds.append(mean_shortest_distance)
        self.algorithm_names.append(str(algorithm_selector))
        self.mean_degrees_results.append(mean_degree)
        self.nums_genes_results.append(len(result_genes))
        self.result_genes.append(','.join(result_genes))
        self.mean_mutual_informations.append(mean_mutual_information)
        self.neg_log_gsea_p_values.append(neg_log_gsea_p_value)

    def run_on_original_network(self, ggi_network_selector, condition_selector, algorithm_selector, verbose):
        """Runs the tests for a given condition on a given original network.

        Parameters
        ----------
        ggi_network_selector : GGINetworkSelector
            Specifies on which GGI network the tests should be run.
        condition_selector : ConditionSelector
            Specifies for which condition the tests should be run.
        algorithm_selector : AlgorithmSelector
            Specifies the algorithm that should be run.
        verbose : bool
            Print progress to stdout.
        """
        ggi_network = utils.load_ggi_network(ggi_network_selector, self.expression_data[condition_selector])
        ggi_network_name = str(ggi_network_selector)
        network_generator_name = str(utils.NetworkGeneratorSelector.ORIGINAL)
        self.run_on_network(ggi_network, -1, ggi_network_name, network_generator_name, condition_selector, algorithm_selector, verbose)

    def run_on_random_networks(self, ggi_network_selector, condition_selector, network_generator_selector,
                               algorithm_selector, verbose):
        """Runs the tests for a given condition on randomized version of a given original network.

        Parameters
        ----------
        ggi_network_selector : GGINetworkSelector
            Specifies on which GGI network the tests should be run.
        condition_selector : ConditionSelector
            Specifies for which condition the tests should be run.
        network_generator_selector : NetworkGeneratorSelector
            Specifies which random generator should be used.
        algorithm_selector : AlgorithmSelector
            Specifies the algorithm that should be run.
        verbose : bool
            Print progress to stdout.
        """
        original_ggi_network = utils.load_ggi_network(ggi_network_selector, self.expression_data[condition_selector])
        num_randomizations = 10
        ggi_network_name = str(ggi_network_selector)
        network_generator_name = str(network_generator_selector)
        seeds = [np.random.randint(low=0, high=np.iinfo(np.uint32).max) for _ in range(num_randomizations)]
        for seed in seeds:
            ggi_network = generators.generate_network(original_ggi_network, seed, network_generator_selector)
            self.run_on_network(ggi_network, seed, ggi_network_name, network_generator_name, condition_selector, algorithm_selector, verbose)


    def clear(self):
        """Clears the results of the last previous run."""
        self.network_generator_names = []
        self.ggi_network_names = []
        self.condition_names = []
        self.algorithm_names = []
        self.random_seeds = []
        self.nums_genes_seeds = []
        self.lcc_ratios_seeds = []
        self.mean_shortest_distances_seeds = []
        self.mean_degrees_results = []
        self.nums_genes_results = []
        self.results = []
        self.mean_mutual_informations = []
        self.neg_log_gsea_p_values = []
        self.results = None
        self.outfile = ''

    def run_all(self, ggi_network_selector, network_generator_selector, algorithm_selector, condition_selectors, verbose):
        """Runs all tests.

        Parameters
        ----------
        ggi_network_selector : GGINetworkSelector
            Specifies on which network the tests should be run.
        network_generator_selector : NetworkGeneratorSelector
            Specifies which random generator should be used.
        algorithm_selector : AlgorithmSelector
            Specifies the algorithm that should be run.
        condition_selectors : list of ConditionSelector or None
            Specifies on which condition the algorithm should be run.
        verbose : bool
            Print progress to stdout.
        """
        self.clear()
        self.outfile = f'../results/{str(ggi_network_selector)}_{str(network_generator_selector)}_{str(algorithm_selector)}.csv'
        if verbose:
            print(f'GGI network = {str(ggi_network_selector)}')
        try:
            if condition_selectors[0] is None:
                condition_selectors = self.condition_selectors
        except TypeError:
            if condition_selectors is None:
                condition_selectors = self.condition_selectors

        for condition_selector in condition_selectors:
            if verbose:
                print(f'\tcondition = {str(condition_selector)}')
            if network_generator_selector == utils.NetworkGeneratorSelector.ORIGINAL:
                self.run_on_original_network(ggi_network_selector, condition_selector, algorithm_selector, verbose)
            else:
                self.run_on_random_networks(ggi_network_selector, condition_selector, network_generator_selector, algorithm_selector, verbose)
        self.results = pd.DataFrame({'network_generator_name': self.network_generator_names,
                                     'ggi_network_name': self.ggi_network_names,
                                     'condition_name': self.condition_names,
                                     'algorithm_name': self.algorithm_names,
                                     'random_seed': self.random_seeds,
                                     'num_genes_seeds': self.nums_genes_seeds,
                                     'lcc_ratio_seeds': self.lcc_ratios_seeds,
                                     'mean_shortest_distance_seeds': self.mean_shortest_distances_seeds,
                                     'mean_degree_result': self.mean_degrees_results,
                                     'num_genes_result': self.nums_genes_results,
                                     'result_genes': self.result_genes,
                                     'mean_mutual_information': self.mean_mutual_informations,
                                     'neg_log_gsea_p_value': self.neg_log_gsea_p_values})

    def get_results(self):
        """Returns the results of the previous run.

        Returns
        -------
        results : pd.DataFrame
            A dataframe containing the results.
        """
        return self.results

    def save_results(self):
        """Writes the results to CSV.

        Parameters
        ----------
        filename : str
            Name of the CSV file to which the results should be written.
        """
        self.results.to_csv(self.outfile, index=False)
