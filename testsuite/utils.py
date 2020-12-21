from enum import Enum
import networkx as nx
import numpy as np
import pandas as pd
import scipy.stats as sps
import itertools as itt
from testsuite.clustex2_wrapper import ClustEx2Wrapper
from testsuite.diamond_wrapper import DIAMOnDWrapper
from testsuite.gxna_wrapper import GXNAWrapper
from testsuite.pinnaclez_wrapper import PinnacleZWrapper
from testsuite.giga_wrapper import GiGAWrapper
from testsuite.kpm_wrapper import KPMWrapper
from testsuite.grandforest_wrapper import GrandForestWrapper
from testsuite.cosine_wrapper import CosineWrapper
from testsuite.hotnet_wrapper import HotNetWrapper
from testsuite.domino_wrapper import DominoWrapper

from testsuite.netcore_wrapper import NetCoreWrapper
from testsuite.custom_wrapper import CustomWrapper
flatten = lambda l: [item for sublist in l for item in sublist]


# todo: add one member for each condition
class ConditionSelector(Enum):
    """Enum specifying for which condition the tests should be run."""
    ALS = 'GSE112680'
    LC = 'GSE30219'
    UC = 'GSE75214'
    CD = 'GSE75214_cd'
    HD = 'GSE3790'

    def __str__(self):
        return self.value


class GGINetworkSelector(Enum):
    """Enum specifying on which GGI network the tests should be run."""
    BIOGRID = 'BIOGRID'
    HPRD = 'HPRD'
    STRING = 'STRING'
    APID = 'APID'
    IID = 'IID'
    IID_BRAIN = "IID_BRAIN"
    IID_LUNG = "IID_LUNG"
    CUSTOM = "CUSTOM"


    def __str__(self):
        return self.value


class NetworkGeneratorSelector(Enum):
    """Enum specifying which random generator should be used."""
    ORIGINAL = 'ORIGINAL'
    EXPECTED_DEGREE = 'EXPECTED_DEGREE'
    SHUFFLED = 'SHUFFLED'
    SCALE_FREE = 'SCALE_FREE'
    UNIFORM = 'UNIFORM'
    REWIRED = "REWIRED"

    def __str__(self):
        return self.value


# todo: add one member for each algorithm
class AlgorithmSelector(Enum):
    """Enum specifying which network enrichment algorithm should be used."""
    KPM = 'KPM'
    DIAMOND = 'DIAMOND'
    GXNA = 'GXNA'
    CLUSTEX2 = 'CLUSTEX2'
    PINNACLEZ = 'PINNACLEZ'
    GIGA = 'GIGA'
    GF = 'GF'
    COSINE = 'COSINE'
    HOTNET = 'HOTNET'
    NETCORE = 'NETCORE'
    CUSTOM = "CUSTOM"
    DOMINO = "DOMINO"

    def __str__(self):
        return self.value


def load_ggi_network(ggi_network_selector, expression_data):
    """Loads the selected GGI network and removes all genes not contained in the expression data.

    Parameters
    ----------
    ggi_network_selector : GGINetworkSelector
        Specifies which GGI network should be loaded.
    expression_data : pd.DataFrame
        Expression data (indices are sample IDs, column names are gene IDs).

    Returns
    -------
    ggi_network : nx.Graph
        The selected GGI network as a networkx graph without genes not contained in the expression data.
    """
    if str(ggi_network_selector) == "CUSTOM":
        edges = pd.read_csv("../data/networks/CUSTOM.tsv", sep = "\t", header = None)
        edges = edges.values
        G = nx.Graph()
        G.add_edges_from(edges)
        all_genes = list(set(flatten(edges)))
        gene_ids = {x: x for x in all_genes}
        nx.set_node_attributes(G, gene_ids, "GeneID")
        nx.write_graphml(G, "../data/networks/CUSTOM.graphml")

    ggi_network = nx.read_graphml(f'../data/networks/{str(ggi_network_selector)}.graphml', node_type=int)
    gene_ids = nx.get_node_attributes(ggi_network, 'GeneID')
    selected_genes = set([str(x) for x in expression_data.columns])
    selected_nodes = [node for node in ggi_network.nodes() if gene_ids[node] in selected_genes]
    ggi_network = ggi_network.subgraph(selected_nodes).copy()
    return nx.convert_node_labels_to_integers(ggi_network)


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
    return np.load(f'../data/conditions/{str(condition_selector)}/phenotype.npy')


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
    expression_data = pd.read_csv(f'../data/expression/{str(condition_selector)}/expr.csv.zip', index_col=0)
    new_index = list(range(len(expression_data.index)))
    expression_data.index = new_index
    return expression_data


def get_pathways(condition_selector):
    """Returns the names of the KEGG pathways associated to the selected condition.

    Parameters
    ----------
    condition_selector : ConditionSelector
        Specifies for which condition the associated pathways should be loaded.

    Returns
    -------
    pathways : list of str
        Names of phenotype-related KEGG pathways.
    """
    if condition_selector == ConditionSelector.ALS:
        return ['hsa05014']
    elif condition_selector == ConditionSelector.LC:
        return ['hsa05223']
    elif condition_selector == ConditionSelector.UC:
        return ['hsa04060', 'hsa04630', 'hsa05321']
    elif condition_selector == ConditionSelector.HD:
        return ['hsa05016']
    elif condition_selector == ConditionSelector.CD:
        return ['hsa04621', 'hsa04060', 'hsa04630', 'hsa05321', 'hsa04140']


# todo: add cases for missing wrappers
def get_algorithm_wrapper(algorithm_selector):
    """Returns the appropriate algorithm based on the selection.

    Parameters
    ----------
    algorithm_selector : AlgorithmSelector
        Specifies which algorithm should be used.
    """
    if algorithm_selector == AlgorithmSelector.GXNA:
        return GXNAWrapper()
    elif algorithm_selector == AlgorithmSelector.CLUSTEX2:
        return ClustEx2Wrapper()
    elif algorithm_selector == AlgorithmSelector.DIAMOND:
        return DIAMOnDWrapper()
    elif algorithm_selector == AlgorithmSelector.PINNACLEZ:
        return PinnacleZWrapper()
    elif algorithm_selector == AlgorithmSelector.GIGA:
        return GiGAWrapper()
    elif algorithm_selector == AlgorithmSelector.KPM:
        return KPMWrapper()
    elif algorithm_selector == AlgorithmSelector.GF:
        return GrandForestWrapper()
    elif algorithm_selector == AlgorithmSelector.COSINE:
        return CosineWrapper()
    elif algorithm_selector == AlgorithmSelector.HOTNET:
        return HotNetWrapper()
    elif algorithm_selector == AlgorithmSelector.NETCORE:
        return NetCoreWrapper()
    elif algorithm_selector == AlgorithmSelector.DOMINO:
        return DominoWrapper()
    elif algorithm_selector == AlgorithmSelector.CUSTOM:
        return CustomWrapper()


def compute_gene_p_values(expression_data, phenotypes):
    """Computes p-values from the expression data and the phenotypes.

    Parameters
    ----------
    expression_data : pd.DataFrame
        Expression data (indices are sample IDs, column names are gene IDs).
    phenotypes : phenotypes : np.array, shape (n_samples,)
        Phenotype data (indices are sample IDs).

    Returns
    -------
    gene_p_values : dict of str: float
            The p-value of two-sided Mann-Whitney U test (keys are gene IDs).
    """
    cases = expression_data.loc[phenotypes == 1, ]
    controls = expression_data.loc[phenotypes == 0, ]
    genes = expression_data.columns
    return {gene: sps.mannwhitneyu(cases[gene], controls[gene], alternative='two-sided')[1] for gene in genes}


def p_value_to_score(p_value):
    """Transforms p-values to scores via -log-transformation.

    Parameters
    ----------
    p_value : float
        The p-value.

    Returns
    -------
    neg_log_p_value : float
        Negative log-transformed p-value.
    """
    return -np.log10(p_value)


# todo: implement this method
def extract_seed_genes(gene_p_values):
    """Extracts the seed genes from previously computed p-values.

    Parameters
    ----------
    gene_p_values : dict of str: float
            The p-value of two-sided Mann-Whitney U test (keys are gene IDs).

    Returns
    -------
    seed_genes : list of str
            List of genes with p-value < 0.001.
    """
    threshold = 0.001 / len(gene_p_values)
    return [gene for gene, p_value in gene_p_values.items() if p_value < threshold]


def compute_indicator_matrix(expression_data, phenotypes):
    """Transforms the expression data to an indicator matrix.

        Parameters
        ----------
        expression_data : pd.DataFrame
            Expression data (indices are sample IDs, column names are gene IDs).
        phenotypes : phenotypes : np.array, shape (n_samples,)
            Phenotype data (indices are sample IDs).

        Returns
        -------
        indicator_matrix : pd.DataFrame
            Indicator matrix obtained from expression data.
        """
    cases = expression_data.loc[phenotypes == 1, ]
    controls = expression_data.loc[phenotypes == 0, ]
    means = np.mean(controls)
    stds = np.std(controls)
    return (np.fabs(cases - means) > 1.5 * stds) * 1


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
    gene_ids = nx.get_node_attributes(ggi_network, 'GeneID')
    seed_nodes = [node for node in ggi_network.nodes() if gene_ids[node] in set(seed_genes)]
    subgraph = ggi_network.subgraph(seed_nodes)
    lcc_ratio = np.max([len(comp) for comp in nx.connected_components(subgraph)]) / subgraph.number_of_nodes()
    sum_shortest_distances = 0
    num_combinations = 0
    for source, target in itt.combinations(seed_nodes, 2):
        try:
            sum_shortest_distances += nx.shortest_path_length(ggi_network, source=source, target=target)
            num_combinations += 1
        except nx.exception.NetworkXNoPath:
            continue
    mean_shortest_distance = sum_shortest_distances / num_combinations
    return lcc_ratio, mean_shortest_distance
