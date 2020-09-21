import networkx as nx
import numpy as np
import testsuite.utils as utils


def generate_rewired_network(ggi_network, seed):
    """Generates random GGI network where nodes keep their degrees but edges are rewired.

    Parameters
    ----------
    ggi_network : nx.Graph
        Original GGI network.
    seed : convertible to np.uint32 or None
        Seed used by the generator.

    Returns
    -------
    rewired_network : nx.Graph
        Random GGI network where nodes keep their degrees but edges are rewired.
    """
    degree_view = ggi_network.degree()
    degree_sequence = [degree_view[node] for node in ggi_network.nodes()]
    rewired_network = nx.random_degree_sequence_graph(degree_sequence, seed=seed, tries=1000)
    gene_ids = nx.get_node_attributes(ggi_network, utils.gene_id_attribute_name())
    nx.set_node_attributes(rewired_network, gene_ids, utils.gene_id_attribute_name())
    return rewired_network


def generate_shuffled_network(ggi_network, seed):
    """Generates random GGI network by shuffling the node labels (gene IDs).

    Parameters
    ----------
    ggi_network : nx.Graph
        Original GGI network.
    seed : convertible to np.uint32 or None
        Seed used by the generator.

    Returns
    -------
    shuffled_network : nx.Graph
        Random GGI network where the nodes labels are shuffled.
    """
    shuffled_network = nx.Graph(ggi_network)
    shuffled_gene_ids = list(nx.get_node_attributes(shuffled_network, utils.gene_id_attribute_name()).values())
    np.random.seed(seed=seed)
    np.random.shuffle(shuffled_gene_ids)
    shuffled_gene_ids = {node: shuffled_gene_ids[node] for node in shuffled_network.nodes()}
    nx.set_node_attributes(shuffled_network, shuffled_gene_ids, utils.gene_id_attribute_name())
    return shuffled_network


def generate_scale_free_network(ggi_network, seed):
    """Generate random scale free network with as many nodes and approximately as many edges as the original network.

    Parameters
    ----------
    ggi_network : nx.Graph
        Original GGI network.
    seed : convertible to np.uint32 or None
        Seed used by the generator.
    Returns
    -------
    scale_free_network : nx.Graph
        Random scale free network with as many nodes and approximately as many edges as the original network.
    """
    num_nodes = ggi_network.number_of_nodes()
    num_edges = ggi_network.number_of_edges()
    m = np.max(1, round(num_nodes / 2.0 - np.sqrt((num_nodes * num_nodes) / 4.0 - num_edges)))
    scale_free_network = nx.barabasi_albert_graph(num_nodes, m, seed=seed)
    gene_ids = nx.get_node_attributes(ggi_network, utils.gene_id_attribute_name())
    nx.set_node_attributes(scale_free_network, gene_ids, utils.gene_id_attribute_name())
    return scale_free_network


def generate_uniform_network(ggi_network, seed):
    """Generates random uniform network with as many nodes and edges as the original network.

    Parameters
    ----------
    ggi_network : nx.Graph
        Original GGI network.
    seed : convertible to np.uint32 or None
        Seed used by the generator.

    Returns
    -------
    uniform_network : nx.Graph
        Random uniform network with as many nodes and edges as the original network.
    """
    num_nodes = ggi_network.number_of_nodes()
    num_edges = ggi_network.number_of_edges()
    uniform_network = nx.gnm_random_graph(num_nodes, num_edges, seed=seed)
    gene_ids = nx.get_node_attributes(ggi_network, utils.gene_id_attribute_name())
    nx.set_node_attributes(uniform_network, gene_ids, utils.gene_id_attribute_name())
    return uniform_network


def generate_network(ggi_network, seed, network_generator_selector):
    """Generates random network based on the original network.

    Parameters
    ----------
    ggi_network : nx.Graph
        Original GGI network.
    seed : convertible to np.uint32 or None
        Seed used by the generator.
    network_generator_selector : utils.NetworkGeneratorSelector
        Specifies which generator should be used.

    Returns
    -------
    random_network : nx.Graph
        Random network that was generated based on the original GGI network.
    """
    if network_generator_selector == utils.NetworkGeneratorSelector.REWIRED:
        return generate_rewired_network(ggi_network, seed)
    elif network_generator_selector == utils.NetworkGeneratorSelector.SHUFFLED:
        return generate_shuffled_network(ggi_network, seed)
    elif network_generator_selector == utils.NetworkGeneratorSelector.SCALE_FREE:
        return generate_scale_free_network(ggi_network, seed)
    elif network_generator_selector == utils.NetworkGeneratorSelector.UNIFORM:
        return generate_uniform_network(ggi_network, seed)