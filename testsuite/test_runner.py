import testsuite.utils as utils
import testsuite.meaningfulness_scores as scores
import testsuite.network_generators as generators
import pandas as pd


class TestRunner(object):

    def __init__(self, ggi_network_selectors, condition_selectors, algorithm_wrappers):
        self.ggi_network_selectors = ggi_network_selectors
        self.ggi_networks = {selector: utils.load_ggi_network(selector) for selector in self.ggi_network_selectors}
        self.condition_selectors = condition_selectors
        self.phenotypes = {selector: utils.load_phenotypes(selector) for selector in self.condition_selectors}
        self.pathways = {selector: utils.load_pathways(selector) for selector in self.condition_selectors}
        self.expression_data = {selector: utils.load_expression_data(selector) for selector in self.condition_selectors}
        self.seed_genes = {selector: utils.extract_seed_genes(self.expression_data[selector], self.phenotypes[selector])
                           for selector in self.condition_selectors}
        self.algorithm_wrappers = algorithm_wrappers
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

    def run_on_network(self, ggi_network, ggi_network_name, network_generator_name, condition_selector):
        phenotypes = self.phenotypes[condition_selector]
        pathways = self.pathways[condition_selector]
        expression_data = self.expression_data[condition_selector]
        seed_genes = self.seed_genes[condition_selector]
        condition_name = str(condition_selector).split('.')[-1]
        lcc_ratio, mean_shortest_distance = utils.compute_seed_statistics(ggi_network, seed_genes)
        for algorithm_wrapper in self.algorithm_wrappers:
            result_genes = algorithm_wrapper.run_algorithm(ggi_network, expression_data, phenotypes, seed_genes)
            mean_mutual_information = scores.compute_mean_mutual_information(expression_data, phenotypes, result_genes)
            neg_log_gsea_p_value = scores.compute_neg_log_gsea_p_value(pathways, result_genes)
            self.ggi_network_names.append(ggi_network_name)
            self.network_generator_names.append(network_generator_name)
            self.condition_names.append(condition_name)
            self.lcc_ratios.append(lcc_ratio)
            self.mean_shortest_distances.append(mean_shortest_distance)
            self.algorithm_names.append(str(algorithm_wrapper))
            self.mean_mutual_informations.append(mean_mutual_information)
            self.neg_log_gsea_p_values.append(neg_log_gsea_p_value)

    def run_on_original_network(self, ggi_network_selector, condition_selector):
        ggi_network = self.ggi_networks[ggi_network_selector]
        ggi_network_name = str(ggi_network_selector).split('.')[-1]
        network_generator_name = 'NONE'
        self.run_on_network(ggi_network, ggi_network_name, network_generator_name, condition_selector)

    def run_on_random_networks(self, ggi_network_selector, condition_selector, network_generator_selector, k):
        original_ggi_network = self.ggi_networks[ggi_network_selector]
        ggi_network_name = str(ggi_network_selector).split('.')[-1]
        network_generator_name = str(network_generator_selector).split('.')[-1]
        for _ in range(k):
            ggi_network = generators.generate_network(original_ggi_network, network_generator_selector)
            self.run_on_network(ggi_network, ggi_network_name, network_generator_name, condition_selector)

    def clear_results(self):
        self.ggi_network_names = []
        self.network_generator_names = []
        self.condition_names = []
        self.lcc_ratios = []
        self.mean_shortest_distances = []
        self.algorithm_names = []
        self.mean_mutual_informations = []
        self.neg_log_gsea_p_values = []

    def run_all(self, k):
        self.clear_results()
        for ggi_network_selector in self.ggi_network_selectors:
            for condition_selector in self.condition_selectors:
                self._run_on_original_network(ggi_network_selector, condition_selector)
                for network_generator_selector in self.network_generator_selectors:
                    self.run_on_random_networks(ggi_network_selector, condition_selector, network_generator_selector, k)

    def save_results(self, filename):
        results = pd.DataFrame({'ggi_network_name': self.ggi_network_names,
                                'network_generator_name': self.network_generator_names,
                                'condition_name': self.condition_names,
                                'lcc_ratio': self.lcc_ratios,
                                'mean_shortest_distance': self.mean_shortest_distances,
                                'algorithm_name': self.algorithm_names,
                                'mean_mutual_information': self.mean_mutual_informations,
                                'neg_log_gsea_p_value': self.neg_log_gsea_p_values})
        results.to_csv(filename)
