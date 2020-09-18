import sklearn.feature_selection as skf
import numpy as np


def compute_mean_mutual_information(expression_data, phenotypes, result_genes):
    """Computes mean mutual information between expression of selected genes and the phenotypes.

    Parameters
    ----------
    expression_data : pd.DataFrame
        Expression data (indices are sample IDs, column names are gene IDs).
    phenotypes : phenotypes : np.array, shape (n_samples,)
        Phenotype data (indices are sample IDs).
    result_genes : list of str
        Set of genes computed by network enrichment algorithm.

    Returns
    -------
    mutual_information_score : float
        Mean mutual information between expression of selected genes and the phenotypes.
    """
    mutual_information = skf.mutual_info_classif(expression_data.loc[:, result_genes], phenotypes)
    return np.mean(mutual_information)


# todo: implement this method
def compute_neg_log_gsea_p_value(pathways, result_genes):
    """Computes gene set enrichment score for result genes w.r.t. phenotype-related pathways.

    Parameters
    ----------
    pathways : list of str
        Names of phenotype-related pathways.
    result_genes : list of str
        Set of genes computed by network enrichment algorithm.

    Returns
    -------
    neg_log_gsea_p_value : float
        Negative log-transformed p-value of gene set enrichment analysis.
    """
    return 0.0
