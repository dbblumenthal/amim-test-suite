import itertools as itt
import pytest
import testsuite.utils as utils
import testsuite.network_generators as generators
import testsuite.meaningfulness_scores as scores


def algorithms():
    return [utils.AlgorithmSelector.HOTNET]
    # return list(utils.AlgorithmSelector)


def ggi_networks():
    return [utils.GGINetworkSelector.HPRD]
    # return list(utils.GGINetworkSelector)


def conditions():
    return [utils.ConditionSelector.CD]
    # return list(utils.ConditionSelector)


def network_generators():
    return [utils.NetworkGeneratorSelector.ORIGINAL]
    # return list(utils.NetworkGeneratorSelector)


def load_data(ggi_network_selector, condition_selector, network_generator_selector):
    print('\tutils.load_expression_data() ...')
    expression_data = utils.load_expression_data(condition_selector)
    print('\tutils.load_ggi_network() ...')
    ggi_network = utils.load_ggi_network(ggi_network_selector, expression_data)
    print('\tutils.utils.load_phenotypes() ...')
    phenotypes = utils.load_phenotypes(condition_selector)
    print('\tutils.compute_gene_scores() ...')
    gene_scores = utils.compute_gene_p_values(expression_data, phenotypes)
    print('\tutils.extract_seed_genes() ...')
    seed_genes = utils.extract_seed_genes(gene_scores)
    print('\tutils.compute_indicator_matrix() ...')
    indicator_matrix = utils.compute_indicator_matrix(expression_data, phenotypes)
    if network_generator_selector != utils.NetworkGeneratorSelector.ORIGINAL:
        print('\tgenerators.generate_network() ...')
        ggi_network = generators.generate_network(ggi_network, None, network_generator_selector)
    return ggi_network, expression_data, phenotypes, seed_genes, gene_scores, indicator_matrix


def run_algorithm(algorithm_wrapper, data, pathways, prefix):
    print('\talgorithm_wrapper.run_algorithm() ...')
    result_genes, mean_degree = algorithm_wrapper.run_algorithm(*data, prefix)
    print('\tscores.compute_mean_mutual_information() ...')
    mean_mutual_information = scores.compute_mean_mutual_information(data[1], data[2], result_genes)
    print('\tscores.compute_neg_log_gsea_p_value() ...')
    neg_log_gsea_p_value = scores.compute_neg_log_gsea_p_value(pathways, result_genes)
    print(f'\tmean_degree = {mean_degree}, mean_mutual_information = {mean_mutual_information}, '
          f'neg_log_gsea_p_value = {neg_log_gsea_p_value}')


@pytest.mark.parametrize(
    'algorithm_selector,ggi_network_selector,condition_selector,network_generator_selector',
    list(itt.product(algorithms(), ggi_networks(), conditions(), network_generators())),
)
def test_algorithm(algorithm_selector, ggi_network_selector, condition_selector, network_generator_selector):
    print(f'{str(algorithm_selector)}-{str(ggi_network_selector)}-{str(condition_selector)}-{str(network_generator_selector)}')
    data = load_data(ggi_network_selector, condition_selector, network_generator_selector)
    print('\tutils.get_algorithm_wrapper() ...')
    algorithm_wrapper = utils.get_algorithm_wrapper(algorithm_selector)
    pathways = utils.get_pathways(condition_selector)
    run_algorithm(algorithm_wrapper, data, pathways, f'{str(ggi_network_selector)}_{str(network_generator_selector)}')
