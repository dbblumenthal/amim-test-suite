import networkx as nx
import numpy as np
import testsuite.utils as utils

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
    edge_list = list(ggi_network.edges())
    A = nx.to_numpy_matrix(ggi_network)
    l = len(edge_list)
    seeds = np.arange(l*100)
    if seed == None:
        seed = 1
    np.random.seed(seed=seed)
    np.random.shuffle(seeds)
    for i in range(100 * l):
        np.random.seed(seed=seeds[i])
        idx = np.random.choice(l-1, 2)
        e1,e2 = [edge_list[i] for i in idx]
        stop = False
        for node1 in e1:
            for node2 in e2:
                if A[node1,node2]:
                    stop = True
        if not stop:
            A[e1[0],e2[1]] = 1
            A[e1[1],e2[0]] = 1
            A[e1[0],e1[1]] = 0
            A[e2[0],e2[1]] = 0
            del edge_list[idx[0]]
            del edge_list[idx[1]]
            edge_list.append((e1[0],e2[1]))
            edge_list.append((e1[1],e2[0]))
    rewired_network = nx.from_numpy_matrix(A)
    gene_ids = nx.get_node_attributes(ggi_network, 'GeneID')
    nx.set_node_attributes(rewired_network, gene_ids, 'GeneID')
    return(rewired_network)

# def generate_complete_network(ggi_network, seed = None):
#     """Generates a fully connected GGI network.
#
#     Parameters
#     ----------
#     ggi_network : nx.Graph
#         Original GGI network.
#     seed : convertible to np.uint32 or None
#         Seed used by the generator.
#
#     Returns
#     -------
#     rewired_network : nx.Graph
#         Random GGI network where nodes keep their expected degrees.
#     """
#     n = len(ggi_network.nodes)
#     A = np.ones((n,n))
#     complete_network = nx.from_numpy_matrix(A)
#     gene_ids = nx.get_node_attributes(ggi_network, 'GeneID')
#     nx.set_node_attributes(complete_network, gene_ids, 'GeneID')
#
#     return complete_network

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
    if network_generator_selector == utils.NetworkGeneratorSelector.REWIRED:
        return generate_rewired_network(ggi_network, seed)
    elif network_generator_selector == utils.NetworkGeneratorSelector.SHUFFLED:
        return generate_shuffled_network(ggi_network, seed)
    elif network_generator_selector == utils.NetworkGeneratorSelector.SCALE_FREE:
        return generate_scale_free_network(ggi_network, seed)
    elif network_generator_selector == utils.NetworkGeneratorSelector.UNIFORM:
        return generate_uniform_network(ggi_network, seed)
    elif network_generator_selector == utils.NetworkGeneratorSelector.RDPN:
        return generate_RDPN(ggi_network, seed)