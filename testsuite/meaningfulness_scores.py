import sklearn.feature_selection as skf
import numpy as np
import gseapy
import mygene
import subprocess



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
    mutual_information = skf.mutual_info_classif(expression_data.loc[:, result_genes], phenotypes, discrete_features=False)
    return np.mean(mutual_information)




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
    mean neg_log_gsea_p_value : float
        Negative log-transformed p-value of gene set enrichment analysis.
    """
    mg = mygene.MyGeneInfo()

    out = mg.querymany(result_genes, scopes= 'entrezgene', fields='symbol', species='human', verbose=False)
    gene_names = []
    for line in out:
        try:
            gene_names.append(line["symbol"])
        except KeyError:
            pass


    res = gseapy.enrichr(gene_list=gene_names, description='pathway', gene_sets='KEGG_2016', cutoff=0.05,
                         outdir="../temp/enrichment", organism='Human', no_plot=True)
    full_results = res.results
    p_values = []
    for tup in full_results.itertuples():
        if tup[10][-8:] in pathways:
            p_values.append(-np.log10(tup[1]))
    subprocess.call('rm ../temp/enrichment/*', shell=True)
    return np.mean(p_values)
