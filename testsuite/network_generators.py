import networkx as nx
import numpy as np
import testsuite.utils as utils
import graph_tool.all as gt

def generate_RDPN(ggi_network, seed):
    """Generates random GGI network where nodes keep the exactly same degrees.

    Parameters
    ----------
    ggi_network : nx.Graph
        Original GGI network.
    seed : convertible to np.uint32 or None
        Seed used by the generator.

    Returns
    -------
    RDPN_network : nx.Graph
        Random GGI network where nodes keep their expected degrees.
    """
    d = nx.to_dict_of_lists(ggi_network)
    edges = [(i, j) for i in d for j in d[i]]
    GT = gt.Graph(directed=False)
    GT.add_vertex(sorted(ggi_network.nodes())[-1])
    GT.add_edge_list(edges)

    gt.random_rewire(GT,model = "constrained-configuration", n_iter = 100, edge_sweep = True)

    edges_new = list(GT.get_edges())
    edges_new = [tuple(x) for x in edges_new]
    rewired_network = nx.Graph()
    rewired_network.add_nodes_from(ggi_network.nodes())
    rewired_network.add_edges_from(edges_new)
    gene_ids = nx.get_node_attributes(ggi_network, 'GeneID')
    nx.set_node_attributes(rewired_network, gene_ids, 'GeneID')
    return(rewired_network)

def generate_rewired_network(ggi_network, seed):
    """Generates random GGI network where nodes keep their expected degrees.

    Parameters
    ----------
    ggi_network : nx.Graph
        Original GGI network.
    seed : convertible to np.uint32 or None
        Seed used by the generator.

    Returns
    -------
    rewired_network : nx.Graph
        Random GGI network where nodes keep their expected degrees.
    """
    degree_view = ggi_network.degree()
    degree_sequence = [degree_view[node] for node in ggi_network.nodes()]
    rewired_network = nx.expected_degree_graph(degree_sequence, seed=seed, selfloops=False)
    gene_ids = nx.get_node_attributes(ggi_network, 'GeneID')
    nx.set_node_attributes(rewired_network, gene_ids, 'GeneID')
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
    shuffled_gene_ids = list(nx.get_node_attributes(shuffled_network, 'GeneID').values())
    np.random.seed(seed=seed)
    np.random.shuffle(shuffled_gene_ids)
    shuffled_gene_ids = {node: shuffled_gene_ids[node] for node in shuffled_network.nodes()}
    nx.set_node_attributes(shuffled_network, shuffled_gene_ids, 'GeneID')
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
    m = np.max([1, round(num_nodes / 2.0 - np.sqrt((num_nodes * num_nodes) / 4.0 - num_edges))])
    scale_free_network = nx.barabasi_albert_graph(num_nodes, m, seed=seed)
    gene_ids = nx.get_node_attributes(ggi_network, 'GeneID')
    nx.set_node_attributes(scale_free_network, gene_ids, 'GeneID')
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
    gene_ids = nx.get_node_attributes(ggi_network, 'GeneID')
    nx.set_node_attributes(uniform_network, gene_ids, 'GeneID')
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
    if network_generator_selector == utils.NetworkGeneratorSelector.EXPECTED_DEGREE:
        return generate_rewired_network(ggi_network, seed)
    elif network_generator_selector == utils.NetworkGeneratorSelector.SHUFFLED:
        return generate_shuffled_network(ggi_network, seed)
    elif network_generator_selector == utils.NetworkGeneratorSelector.SCALE_FREE:
        return generate_scale_free_network(ggi_network, seed)
    elif network_generator_selector == utils.NetworkGeneratorSelector.UNIFORM:
        return generate_uniform_network(ggi_network, seed)
    elif network_generator_selector == utils.NetworkGeneratorSelector.REWIRED:
        return generate_RDPN(ggi_network, seed)