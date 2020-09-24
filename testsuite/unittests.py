import itertools as itt
import testsuite.utils as utils
import testsuite.network_generators as generators
import testsuite.meaningfulness_scores as scores


def algorithms():
    return [utils.AlgorithmSelector.GXNA]
    #return list(utils.AlgorithmSelector)


def ggi_networks():
    return [utils.GGINetworkSelector.HPRD]
    #return list(utils.GGINetworkSelector)


def conditions():
    return [utils.ConditionSelector.LC]
    #return list(utils.ConditionSelector)


def network_generators():
    return list(utils.NetworkGeneratorSelector)


def run_algorithms(algorithm_selector, ggi_network_selector, condition_selector, network_generator_selector):
    print(f'algorithm_selector = {algorithm_selector}')
    print(f'ggi_network_selector = {ggi_network_selector}')
    print(f'condition_selector = {condition_selector}')
    print(f'network_generator_selector = {network_generator_selector}')
    print('\tutils.load_expression_data() ...')
    expression_data = utils.load_expression_data(condition_selector)
    print('\tutils.load_ggi_network() ...')
    ggi_network = utils.load_ggi_network(ggi_network_selector, expression_data)
    print('\tutils.utils.load_phenotypes() ...')
    phenotypes = utils.load_phenotypes(condition_selector)
    print('\tutils.compute_gene_scores() ...')
    gene_scores = utils.compute_gene_scores(expression_data, phenotypes)
    print('\tutils.extract_seed_genes() ...')
    seed_genes = utils.extract_seed_genes(gene_scores)
    if network_generator_selector != utils.NetworkGeneratorSelector.ORIGINAL:
        print('\tgenerators.generate_network() ...')
        ggi_network = generators.generate_network(ggi_network, None, network_generator_selector)
    print('\tutils.get_algorithm_wrapper ...')
    algorithm_wrapper = utils.get_algorithm_wrapper(algorithm_selector)
    print('\talgorithm_wrapper.run_algorithm() ...')
    result_genes, mean_degree = algorithm_wrapper.run_algorithm(ggi_network,expression_data,phenotypes, seed_genes, gene_scores)
    print('\tscores.compute_mean_mutual_information() ...')
    mean_mutual_information = scores.compute_mean_mutual_information(expression_data, phenotypes, result_genes)
    print(f'\tmean_degree = {mean_degree}, mean_mutual_information = {mean_mutual_information}')


if __name__ == '__main__':
    configs = list(itt.product(algorithms(), ggi_networks(), conditions(), network_generators()))
    for config in configs:
        run_algorithms(*config)