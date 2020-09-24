import testsuite.utils as utils
import testsuite.meaningfulness_scores as scores
import testsuite.network_generators as generators
import pandas as pd
import numpy as np
import argparse


class TestRunner(object):
    """Runs the tests."""

    def __init__(self):
        """Constructs TestRunner object."""
        self.ggi_network_selectors = list(utils.GGINetworkSelector)
        self.condition_selectors = list(utils.ConditionSelector)
        self.algorithm_selectors = list(utils.AlgorithmSelector)
        self.ggi_networks = {sel: utils.load_ggi_network(sel) for sel in self.ggi_network_selectors}
        self.phenotypes = {sel: utils.load_phenotypes(sel) for sel in self.condition_selectors}
        self.pathways = {sel: utils.load_pathways(sel) for sel in self.condition_selectors}
        self.expression_data = {sel: utils.load_expression_data(sel) for sel in self.condition_selectors}
        self.gene_scores = {sel: utils.compute_gene_scores(self.expression_data[sel], self.phenotypes[sel]) for sel in self.condition_selectors}
        self.seed_genes = {sel: utils.extract_seed_genes(self.gene_scores[sel]) for sel in self.condition_selectors}
        self.algorithm_wrappers = {sel: utils.get_algorithm_wrapper(sel) for sel in self.algorithm_selector}
        self.network_generator_names = []
        self.ggi_network_names = []
        self.condition_names = []
        self.algorithm_names = []
        self.random_seeds = []
        self.nums_seed_genes = []
        self.lcc_ratios = []
        self.mean_shortest_distances = []
        self.mean_mutual_informations = []
        self.neg_log_gsea_p_values = []
        self.results = None
        self.outfile = ''

    def run_on_network(self, ggi_network, seed, ggi_network_name, network_generator_name, condition_selector, num_runs, verbose):
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
        num_runs : int
            Number of runs for each algorithm.
        verbose : bool
            Print progress to stdout.
        """
        phenotypes = self.phenotypes[condition_selector]
        pathways = self.pathways[condition_selector]
        expression_data = self.expression_data[condition_selector]
        seed_genes = self.seed_genes[condition_selector]
        lcc_ratio, mean_shortest_distance = utils.compute_seed_statistics(ggi_network, seed_genes)
        for algorithm_selector in self.algorithm_selectors:
            if verbose:
                print(f'\t\talgorithm = {str(algorithm_selector)}')
            algorithm_wrapper = self.algorithm_wrappers[algorithm_selector]
            for run in range(num_runs):
                print(f'\t\t\trun {run + 1} of {num_runs}')
                result_genes = algorithm_wrapper.run_algorithm(ggi_network, expression_data, phenotypes, seed_genes)
                mean_mutual_information = scores.compute_mean_mutual_information(expression_data, phenotypes, result_genes)
                neg_log_gsea_p_value = scores.compute_neg_log_gsea_p_value(pathways, result_genes)
                self.ggi_network_names.append(ggi_network_name)
                self.random_seeds.append(seed)
                self.network_generator_names.append(network_generator_name)
                self.condition_names.append(str(condition_selector))
                self.nums_seed_genes.append(len(seed_genes))
                self.lcc_ratios.append(lcc_ratio)
                self.mean_shortest_distances.append(mean_shortest_distance)
                self.algorithm_names.append(str(algorithm_selector))
                self.mean_mutual_informations.append(mean_mutual_information)
                self.neg_log_gsea_p_values.append(neg_log_gsea_p_value)

    def run_on_original_network(self, ggi_network_selector, condition_selector, num_runs, verbose):
        """Runs the tests for a given condition on a given original network.

        Parameters
        ----------
        ggi_network_selector : GGINetworkSelector
            Specifies on which GGI network the tests should be run.
        condition_selector : ConditionSelector
            Specifies for which condition the tests should be run.
        num_runs : int
            Number of runs for each algorithm.
        verbose : bool
            Print progress to stdout.
        """
        ggi_network = self.ggi_networks[ggi_network_selector]
        ggi_network_name = str(ggi_network_selector)
        network_generator_name = str(utils.NetworkGeneratorSelector.ORIGINAL)
        self.run_on_network(ggi_network, -1, ggi_network_name, network_generator_name, condition_selector, num_runs, verbose)

    def run_on_random_networks(self, ggi_network_selector, condition_selector, network_generator_selector,
                               num_randomizations, num_runs, verbose):
        """Runs the tests for a given condition on randomized version of a given original network.

        Parameters
        ----------
        ggi_network_selector : GGINetworkSelector
            Specifies on which GGI network the tests should be run.
        condition_selector : ConditionSelector
            Specifies for which condition the tests should be run.
        network_generator_selector : NetworkGeneratorSelector
            Specifies which random generator should be used.
        num_randomizations : int
            Specifies how many randomized versions of the original GGI networks should be generated.
        verbose : bool
            Print progress to stdout.
        """
        original_ggi_network = self.ggi_networks[ggi_network_selector]
        ggi_network_name = str(ggi_network_selector)
        network_generator_name = str(network_generator_selector)
        seeds = [np.random.randint(low=0, high=np.iinfo(np.uint32).max) for _ in range(num_randomizations)]
        for seed in seeds:
            ggi_network = generators.generate_network(original_ggi_network, seed, network_generator_selector)
            self.run_on_network(ggi_network, seed, ggi_network_name, network_generator_name, condition_selector, num_runs, verbose)

    def clear(self):
        """Clears the results of the last previous run."""
        self.ggi_network_names = []
        self.network_generator_names = []
        self.condition_names = []
        self.lcc_ratios = []
        self.mean_shortest_distances = []
        self.algorithm_names = []
        self.mean_mutual_informations = []
        self.neg_log_gsea_p_values = []
        self.results = None
        self.outfile = ''

    def run_all(self, num_randomizations, num_runs, network_generator_selector, verbose):
        """Runs all tests.

        Parameters
        ----------
        num_randomizations : int
            Specifies how many randomized versions of the original GGI networks should be generated.
        num_runs : int
            Number of runs for each algorithm.
        network_generator_selector : NetworkGeneratorSelector
            Specifies which random generator should be used.
        verbose : bool
            Print progress to stdout.
        """
        self.clear()
        self.outfile = f'../results/{str(network_generator_selector)}.csv'
        for ggi_network_selector in self.ggi_network_selectors:
            if verbose:
                print(f'GGI network = {str(ggi_network_selector)}')
            for condition_selector in self.condition_selectors:
                if verbose:
                    print(f'\tcondition = {str(condition_selector)}')
                if network_generator_selector == utils.NetworkGeneratorSelector.ORIGINAL:
                    self.run_on_original_network(ggi_network_selector, condition_selector, num_runs, verbose)
                else:
                    self.run_on_random_networks(ggi_network_selector, condition_selector, network_generator_selector,
                                                num_randomizations, num_runs, verbose)
        self.results = pd.DataFrame({'network_generator_name': self.network_generator_names,
                                     'ggi_network_name': self.ggi_network_names,
                                     'condition_name': self.condition_names,
                                     'algorithm_name': self.algorithm_names,
                                     'random_seed': self.random_seeds,
                                     'num_seed_genes': self.nums_seed_genes,
                                     'lcc_ratio': self.lcc_ratios,
                                     'mean_shortest_distance': self.mean_shortest_distances,
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
        self.results.to_csv(self.outfile)


if __name__ == '__main__':
    parser = argparse.ArgumentParser('tests the one-network-fits-all hypothesis')
    parser.add_argument('generator', type=utils.NetworkGeneratorSelector, choices=list(utils.NetworkGeneratorSelector))
    parser.add_argument('--num-randomizations', type=int, default=5, help='number of network randomizations')
    parser.add_argument('--num-runs', type=int, default=1, help='number of runs for each algorithm')
    parser.add_argument('--verbose', action='store_true', help='print progress to stdout')
    args = parser.parse_args()
    test_runner = TestRunner()
    test_runner.run_all(args.num_randomizations, args.num_runs, args.generator, args.verbose)
    test_runner.save_results()
