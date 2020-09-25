import testsuite.utils as utils
import pandas as pd
import scipy.stats as sps
import numpy as np
import argparse


class ResultsAnalyzer(object):

    def __init__(self):
        """Constructs a ResultsAnalyzer object."""
        self.network_generator_selectors = [selector for selector in list(utils.NetworkGeneratorSelector)
                                            if not selector == utils.NetworkGeneratorSelector.ORIGINAL]
        self.ggi_network_selectors = list(utils.GGINetworkSelector)
        self.condition_selectors = list(utils.ConditionSelector)
        self.algorithm_selectors = list(utils.AlgorithmSelector)
        self.network_generator_names = []
        self.ggi_network_names = []
        self.condition_names = []
        self.algorithm_names = []
        self.means_num_seed_genes = []
        self.deltas_mean_lcc_ratios = []
        self.deltas_mean_shortest_distances = []
        self.p_values_mi = []
        self.p_values_gsea = []
        self.analyzed_results = None

    def clear(self):
        """Clears the results of the previous analysis."""
        self.network_generator_names = []
        self.ggi_network_names = []
        self.condition_names = []
        self.algorithm_names = []
        self.means_num_seed_genes = []
        self.deltas_mean_lcc_ratios = []
        self.deltas_mean_shortest_distances = []
        self.p_values_mi = []
        self.p_values_gsea = []
        self.analyzed_results = None

    def update_meta_information(self, network_generator_name, ggi_network_name, condition_name, algorithm_name):
        """Updates the meta information in the result lists.

        Parameters
        ----------
        network_generator_name : str
            Name of the network generator used for creating the results view.
        ggi_network_name
            Name of the GGI network used for creating the results view or 'ALL'.
        condition_name
            Name of the condition used for creating the results view or 'ALL'
        algorithm_name
            Name of the algorithm used for creating the results view or 'ALL'
        """
        self.network_generator_names.append(network_generator_name)
        self.ggi_network_names.append(ggi_network_name)
        self.condition_names.append(condition_name)
        self.algorithm_names.append(algorithm_name)

    def compute_test_statistics(self, original, randomized):
        """Computes the test-statistics for the given results view.

        Parameters
        ----------
        original : pd.DataFrame
            View of the results for the original GGI networks.
        randomized : pd.DataFrame
            View of the results for the randomized GGI networks.
        """
        self.means_num_seed_genes.append(np.mean(original['num_seed_genes']))
        self.deltas_mean_lcc_ratios.append(np.mean(original['lcc_ratio']) - np.mean(randomized['lcc_ratio']))
        self.deltas_mean_shortest_distances.append(np.mean(original['mean_shortest_distance']) - np.mean(randomized['mean_shortest_distance']))
        p_value_mi = sps.mannwhitneyu(original['mean_mutual_information'], randomized['mean_mutual_information'],
                                           alternative='greater')[1]
        p_value_gsea = sps.mannwhitneyu(original['neg_log_gsea_p_value'], randomized['neg_log_gsea_p_value'],
                                             alternative='greater')[1]
        self.p_values_mi.append(p_value_mi)
        self.p_values_gsea.append(p_value_gsea)

    def analyze_results(self, verbose):
        """Analyzes the results produced by TestRunner.

        Parameters
        ----------
        results : pd.DataFrame or str
            Results produced by TestRunner or path to CSV file.
        verbose : bool
            Print progress to stdout.
        """
        self.clear()
        results_for_original_networks = pd.read_csv('../results/ORIGINAL.csv')
        for network_generator_selector in self.network_generator_selectors:
            network_generator_name = str(network_generator_selector)
            if verbose:
                print(f'generator = {network_generator_name}')
                print('\tselect all')
            results_for_randomized_networks = pd.read_csv(f'../results/{network_generator_name}.csv')
            self.update_meta_information(network_generator_name, 'ALL', 'ALL', 'ALL')
            self.compute_test_statistics(results_for_original_networks, results_for_randomized_networks)
            for ggi_network_selector in self.ggi_network_selectors:
                ggi_network_name = str(ggi_network_selector)
                if verbose:
                    print(f'\tselect network = {ggi_network_name}')
                original = results_for_original_networks.loc[results_for_original_networks['ggi_network_name'] == ggi_network_name]
                randomized = results_for_randomized_networks.loc[results_for_randomized_networks['ggi_network_name'] == ggi_network_name]
                self.update_meta_information(network_generator_name, ggi_network_name, 'ALL', 'ALL')
                self.compute_test_statistics(original, randomized)
            for condition_selector in self.condition_selectors:
                condition_name = str(condition_selector)
                if verbose:
                    print(f'\tselect condition = {condition_name}')
                original = results_for_original_networks.loc[results_for_original_networks['condition_name'] == condition_name]
                randomized = results_for_randomized_networks.loc[results_for_randomized_networks['condition_name'] == condition_name]
                self.update_meta_information(network_generator_name, 'ALL', condition_name, 'ALL')
                self.compute_test_statistics(original, randomized)
            for algorithm_selector in self.algorithm_selectors:
                algorithm_name = str(algorithm_selector)
                if verbose:
                    print(f'\tselect algorithm = {algorithm_name}')
                original = results_for_original_networks.loc[results_for_original_networks['algorithm_name'] == algorithm_name]
                randomized = results_for_randomized_networks.loc[results_for_randomized_networks['algorithm_name'] == algorithm_name]
                self.update_meta_information(network_generator_name, 'ALL', 'ALL', algorithm_name)
                self.compute_test_statistics(original, randomized)
        self.analyzed_results = pd.DataFrame({'network_generator_name': self.network_generator_names,
                                              'ggi_network_name': self.ggi_network_names,
                                              'condition_name': self.condition_names,
                                              'algorithm_name': self.algorithm_names,
                                              'mean_num_seed_genes': self.means_num_seed_genes,
                                              'delta_mean_lcc_ratios': self.deltas_mean_lcc_ratios,
                                              'delta_mean_shortest_path_distances': self.deltas_mean_shortest_distances,
                                              'p_value_mi': self.p_values_mi,
                                              'p_value_gsea': self.p_values_gsea})

    def get_analyzed_results(self):
        """Returns the analyzed results.

        Returns
        -------
        analyzed_results : pd.DataFrame
            A dataframe containing the analyzed results.
        """
        return self.analyzed_results

    def save_analyzed_results(self, filename):
        """Writes the analyzed results to CSV.

        Parameters
        ----------
        filename : str
            Name of the CSV file to which the analyzed results should be written.
        """
        self.analyzed_results.to_csv(filename)


if __name__ == '__main__':
    parser = argparse.ArgumentParser('analyzes the results')
    parser.add_argument('--verbose', action='store_true', help='print progress to stdout')
    args = parser.parse_args()
    analyzer = ResultsAnalyzer()
    analyzer.analyze_results(args.verbose)
    analyzer.save_analyzed_results(args.output)
