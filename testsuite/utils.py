from enum import Enum, auto
import networkx as nx
import numpy as np
import itertools as itt


# todo: add one member for each condition
class ConditionSelector(Enum):
    """Enum specifying for which condition the tests should be run."""
    LUNG_CANCER = 'LUNG_CANCER'

    def __str__(self):
        return self.value


class GGINetworkSelector(Enum):
    """Enum specifying on which GGI network the tests should be run."""
    BioGRID = 'BioGRID'
    HPRD = 'HPRD'
    STRING = 'STRING'
    APID = 'APID'
    IID = 'IID'

    def __str__(self):
        return self.value


class NetworkGeneratorSelector(Enum):
    """Enum specifying which random generator should be used."""
    REWIRED = 'REWIRED'
    SHUFFLED = 'SHUFFLED'
    SCALE_FREE = 'SCALE_FREE'
    UNIFORM = 'UNIFORM'

    def __str__(self):
        return self.value


# todo: add one member for each algorithm
class AlgorithmSelector(Enum):
    """Enum specifying which network enrichment algorithm should be used."""
    DIAMOND = 'DIAMOND'

    def __str__(self):
        return self.value


def gene_id_attribute_name():
    """Returns the name of the gene ID attribute in the networkx graphs.

    Returns
    -------
    gene_id_attribute_name : str
        Name of the gene ID attribute in the networkx graphs.
    """
    return 'GeneID'


# todo: implement this method
def load_ggi_network(ggi_network_selector):
    """Loads the selected GGI network.

    Parameters
    ----------
    ggi_network_selector : GGINetworkSelector
        Specifies which GGI network should be loaded.

    Returns
    -------
    ggi_network : nx.Graph
        The selected GGI network as a networkx graph.
    """
    pass


# todo: implement this method
def load_phenotypes(condition_selector):
    """Loads the phenotypes for the selected condition.

    Parameters
    ----------
    condition_selector : ConditionSelector
        Specifies for which condition the phenotypes should be loaded.

    Returns
    -------
    phenotypes : phenotypes : np.array, shape (n_samples,)
        Phenotype data (indices are sample IDs).
    """
    pass


# todo: implement this method
def load_expression_data(condition_selector):
    """Loads the expression data for the selected condition.

    Parameters
    ----------
    condition_selector : ConditionSelector
        Specifies for which condition the phenotypes should be loaded.

    Returns
    -------
    expression_data : pd.DataFrame
        Expression data (indices are sample IDs, column names are gene IDs).
    """
    pass


# todo: implement this method
def load_pathways(condition_selector):
    """Loads the pathways associated to the selected condition.

    Parameters
    ----------
    condition_selector : ConditionSelector
        Specifies for which condition the associated pathways should be loaded.

    Returns
    -------
    pathways : list of str
        Names of phenotype-related pathways.
    """
    pass


# todo: implement this method
def get_algorithm_wrapper(algorithm_selector):
    """Returns the appropriate algorithm based on the selection.

    Parameters
    ----------
    algorithm_selector : AlgorithmSelector
        Specifies which algorithm should be used.
    """
    pass


# todo: implement this method
def extract_seed_genes(expression_data, phenotypes):
    """Extracts the seed genes from the expression data and the phenotypes.

    Parameters
    ----------
    expression_data : pd.DataFrame
        Expression data (indices are sample IDs, column names are gene IDs).
    phenotypes : phenotypes : np.array, shape (n_samples,)
        Phenotype data (indices are sample IDs).

    Returns
    -------
    seed_genes : list of str
            Seed genes (entries are gene IDs).
    """
    return []


def compute_seed_statistics(ggi_network, seed_genes):
    """Computes the seed genes' LCC ratio and the mean shortest distance between the seed genes.

    Parameters
    ----------
    ggi_network : nx.Graph
        Original GGI network.
    seed_genes : list of str
        Seed genes (entries are gene IDs).

    Returns
    -------
    lcc_ratio : float
        The ratio of nodes in subgraph induced by seed genes which are contained in largest connected component.
    mean_shortest_distance : float
        The mean shortest distance between the seed genes in the GGI network.
    """
    gene_ids = nx.get_node_attributes(ggi_network, gene_id_attribute_name())
    seed_nodes = [node for node in ggi_network.nodes() if gene_ids[node] in set(seed_genes)]
    subgraph = ggi_network.subgraph(seed_nodes)
    lcc_ratio = np.max([len(comp) for comp in nx.connected_components(subgraph)]) / subgraph.number_of_nodes()
    sum_shortest_distances = 0
    num_combinations = 0
    for source, target in itt.combinations(seed_genes, 2):
        sum_shortest_distances += nx.shortest_path_length(ggi_network, source=source, target=target)
        num_combinations += 1
    mean_shortest_distance = sum_shortest_distances / num_combinations
    return lcc_ratio, mean_shortest_distance
