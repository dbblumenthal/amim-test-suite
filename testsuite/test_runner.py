import testsuite.utils as utils
import testsuite.meaningfulness_scores as scores
import testsuite.network_generators as generators
import pandas as pd
import argparse


class TestRunner(object):
    """Runs the tests."""

    def __init__(self, ggi_network_selectors, condition_selectors, algorithm_selectors):
        """Constructs TestRunner object.red = 'red'

        Parameters
        ----------
        ggi_network_selectors : list of GGINetworkSelector
            Specifies on which GGI networks the tests should be run.
        condition_selectors : list of ConditionSelector
            Specifies on which conditions the tests should be run.
        algorithm_selectors : list of AlgorithmSelector
            Specifies which algorithms should be employed for the tests.
        """
        self.ggi_network_selectors = ggi_network_selectors
        self.ggi_networks = {selector: utils.load_ggi_network(selector) for selector in self.ggi_network_selectors}
        self.condition_selectors = condition_selectors
        self.phenotypes = {selector: utils.load_phenotypes(selector) for selector in self.condition_selectors}
        self.pathways = {selector: utils.load_pathways(selector) for selector in self.condition_selectors}
        self.expression_data = {selector: utils.load_expression_data(selector) for selector in self.condition_selectors}
        self.seed_genes = {selector: utils.extract_seed_genes(self.expression_data[selector], self.phenotypes[selector])
                           for selector in self.condition_selectors}
        self.algorithm_selectors = algorithm_selectors
        self.algorithm_wrappers = {selector: utils.get_algorithm_wrapper(selector) for selector in self.algorithm_selector}
        self.network_generator_selectors = [utils.NetworkGeneratorSelector.REWIRED,
                                            utils.NetworkGeneratorSelector.SHUFFLED,
                                            utils.NetworkGeneratorSelector.SCALE_FREE,
                                            utils.NetworkGeneratorSelector.UNIFORM]
        self.ggi_network_names = []
        self.network_generator_names = []
        self.condition_names = []
        self.lcc_ratios = []
        self.mean_shortest_distances = []
        self.algorithm_names = []
        self.mean_mutual_informations = []
        self.neg_log_gsea_p_values = []
        self.results = None

    def run_on_network(self, ggi_network, ggi_network_name, network_generator_name, condition_selector):
        """Runs the tests for a given condition on a given network.

        Parameters
        ----------
        ggi_network : nx.Graph
            Possibly randomized GGI network.
        ggi_network_name : str
            String representation of the original GGI network.
        network_generator_name : str
            String representation of the employed network generator.
        condition_selector : ConditionSelector
            Specifies for which condition the tests should be run.
        """
        phenotypes = self.phenotypes[condition_selector]
        pathways = self.pathways[condition_selector]
        expression_data = self.expression_data[condition_selector]
        seed_genes = self.seed_genes[condition_selector]
        lcc_ratio, mean_shortest_distance = utils.compute_seed_statistics(ggi_network, seed_genes)
        for algorithm_selector in self.algorithm_selectors:
            algorithm_wrapper = self.algorithm_wrappers[algorithm_selector]
            result_genes = algorithm_wrapper.run_algorithm(ggi_network, expression_data, phenotypes, seed_genes)
            mean_mutual_information = scores.compute_mean_mutual_information(expression_data, phenotypes, result_genes)
            neg_log_gsea_p_value = scores.compute_neg_log_gsea_p_value(pathways, result_genes)
            self.ggi_network_names.append(ggi_network_name)
            self.network_generator_names.append(network_generator_name)
            self.condition_names.append(str(condition_selector))
            self.lcc_ratios.append(lcc_ratio)
            self.mean_shortest_distances.append(mean_shortest_distance)
            self.algorithm_names.append(str(algorithm_selector))
            self.mean_mutual_informations.append(mean_mutual_information)
            self.neg_log_gsea_p_values.append(neg_log_gsea_p_value)

    def run_on_original_network(self, ggi_network_selector, condition_selector):
        """Runs the tests for a given condition on a given original network.

        Parameters
        ----------
        ggi_network_selector : GGINetworkSelector
            Specifies on which GGI network the tests should be run.
        condition_selector : ConditionSelector
            Specifies for which condition the tests should be run.
        """
        ggi_network = self.ggi_networks[ggi_network_selector]
        ggi_network_name = str(ggi_network_selector)
        network_generator_name = 'NONE'
        self.run_on_network(ggi_network, ggi_network_name, network_generator_name, condition_selector)

    def run_on_random_networks(self, ggi_network_selector, condition_selector, network_generator_selector, k):
        """Runs the tests for a given condition on randomized version of a given original network.

        Parameters
        ----------
        ggi_network_selector : GGINetworkSelector
            Specifies on which GGI network the tests should be run.
        condition_selector : ConditionSelector
            Specifies for which condition the tests should be run.
        network_generator_selector : NetworkGeneratorSelector
            Specifies which random generator should be used.
        k : int
            Specifies how many randomized versions of the original GGI networks should be generated.
        """
        original_ggi_network = self.ggi_networks[ggi_network_selector]
        ggi_network_name = str(ggi_network_selector)
        network_generator_name = str(network_generator_selector)
        for _ in range(k):
            ggi_network = generators.generate_network(original_ggi_network, network_generator_selector)
            self.run_on_network(ggi_network, ggi_network_name, network_generator_name, condition_selector)

    def clear_results(self):
        """Clears the results of the last previous run."""
        self.ggi_network_names = []
        self.network_generator_names = []
        self.condition_names = []
        self.lcc_ratios = []
        self.mean_shortest_distances = []
        self.algorithm_names = []
        self.mean_mutual_informations = []
        self.neg_log_gsea_p_values = []

    def run_all(self, k, verbose):
        """Runs all tests.

        Parameters
        ----------
        k : int
            Specifies how many randomized versions of the original GGI networks should be generated.
        verbose : bool
            Print progress to stdout.
        """
        self.clear_results()
        for ggi_network_selector in self.ggi_network_selectors:
            if verbose:
                print(f'GGI network = {ggi_network_selector}')
            for condition_selector in self.condition_selectors:
                if verbose:
                    print(f'\tcondition = {condition_selector}')
                    print('\t\tgenerator = NONE')
                self._run_on_original_network(ggi_network_selector, condition_selector)
                for network_generator_selector in self.network_generator_selectors:
                    if verbose:
                        print(f'\t\tgenerator = {network_generator_selector}')
                    self.run_on_random_networks(ggi_network_selector, condition_selector, network_generator_selector, k)
        self.results = pd.DataFrame({'ggi_network_name': self.ggi_network_names,
                                     'network_generator_name': self.network_generator_names,
                                     'condition_name': self.condition_names,
                                     'lcc_ratio': self.lcc_ratios,
                                     'mean_shortest_distance': self.mean_shortest_distances,
                                     'algorithm_name': self.algorithm_names,
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

    def save_results(self, filename):
        """Writes the results to CSV.

        Parameters
        ----------
        filename : str
            Name of the CSV file to which the results should be written.
        """
        self.results.to_csv(filename)


if __name__ == '__main__':
    parser = argparse.ArgumentParser('tests the one-network-fits-all hypothesis')
    parser.add_argument('outfilename', help='name of output file name (CSV)')
    parser.add_argument('--ggis', nargs='+', type=utils.GGINetworkSelector, choices=list(utils.GGINetworkSelector),
                        default=list(utils.GGINetworkSelector), help='employed GGI networks')
    parser.add_argument('--conditions', nargs='+', type=utils.ConditionSelector, choices=list(utils.ConditionSelector),
                        default=list(utils.ConditionSelector), help='conditions for which the tests should be run')
    parser.add_argument('--algorithms', nargs='+', type=utils.AlgorithmSelector, choices=list(utils.AlgorithmSelector),
                        default=list(utils.AlgorithmSelector), help='employed algorithms')
    parser.add_argument('--k', type=int, default=5, help='number of randomizations for each network generator')
    parser.add_argument('--verbose', action='store_true', help='print progress to stdout')
    args = parser.parse_args()
    test_runner = TestRunner(args.ggis, args.conditions, args.algorithms)
    test_runner.run_all(args.k, args.verbose)
    test_runner.save_results(args.outfilename)
