from testsuite.algorithm_wrapper import AlgorithmWrapper
import networkx as nx
import subprocess
import pandas as pd
import mygene
flatten = lambda l: [item for sublist in l for item in sublist]

#
# import os
# from testsuite.utils import *
#
# os.chdir('/home/olga/Dropbox/testing-onfah/testsuite')
# expression_data = load_expression_data("GSE3790")
# phenotypes = load_phenotypes("GSE3790")
# ggi_network = load_ggi_network("HPRD", expression_data)
# prefix = 'HPRD_ORIGINAL'
# pvs = compute_gene_p_values(expression_data, phenotypes)
# seed_genes = extract_seed_genes(pvs)

class DominoWrapper(AlgorithmWrapper):

    def run_algorithm(self, ggi_network, expression_data, phenotypes, seed_genes, p_values, indicator_matrix, prefix):
        """Runs the algorithm.

        Parameters
        ----------
        ggi_network : nx.Graph
            Possibly randomized GGI network.
        expression_data : pd.DataFrame
            Expression data (indices are sample IDs, column names are gene IDs).
        phenotypes : np.array, shape (n_samples,)
            Phenotype data (indices are sample IDs).
        seed_genes : list of str
            Seed genes (entries are gene IDs).
        p_values : dict of str: float
            P-values for all genes (keys are gene IDs).
        indicator_matrix : pd.DataFrame
            Indicator matrix obtained from expression data (indices are sample IDs, column names are gene IDs).
        prefix : str
            Prefix to be used for temporary files and directories.

        Returns
        -------
        result_genes : list of str
            Set of genes computed by the algorithm.
        mean_degree : float
            Mean degree of the result genes.
        """


        genes = list(expression_data.columns)
        mg = mygene.MyGeneInfo()
        out = mg.querymany(genes, scopes="entrezgene", fields='ensembl.gene', species='human', verbose=False)
        mapping = dict()
        for line in out:
            try:
                res = line["ensembl"]
                if len(res) == 1:

                    mapping[line["query"]] = res["gene"]
                else:
                    mapping[line["query"]] = res[0]["gene"]

            except KeyError:
                mapping[line["query"]] = "pass"
        # 1. Write GGI network in format required by your methodmapping[line["query"]] = res["gene"]
        path_to_network = f'../temp/{prefix}_domino_ggi.cif'
        gene_ids = nx.get_node_attributes(ggi_network, 'GeneID')

        with open(path_to_network, 'w') as edge_list_file:
            edge_list_file.write(f'node_1\t combined_score \tnode_2\n')
            for u, v in ggi_network.edges():
                edge_list_file.write(f'{mapping[gene_ids[u]]}\t ppi \t{mapping[gene_ids[v]]}\n')
        outfile = f'../../temp/{prefix}_domino_slices.txt'
        path_to_network = f'../../temp/{prefix}_domino_ggi.cif'

        command = f'cd ../algorithms/DOMINO; slicer --network_file {path_to_network} --output_file {outfile}'
        subprocess.call(command, shell = True)

        # 2. Write expression data and phenotype in the format required by your method

        seeds_ens = [mapping[x] for x in seed_genes if mapping[x] != "pass"]
        path_to_seeds = f'../temp/{prefix}_domino_seeds.txt'
        with open(path_to_seeds, 'w') as seeds_file:
            for g in seeds_ens:
                seeds_file.write(f'{g}\n')
        path_to_seeds = f'../../temp/{prefix}_domino_seeds.txt'

        # 3. Insert the command to run your method, direct the output to path_to_output
        path_to_output = f'../../temp/'
        command = f'cd ../algorithms/DOMINO; domino --active_genes_files {path_to_seeds} --network_file {path_to_network} --slices_file {outfile} --output_folder {path_to_output} --visualization false'
        subprocess.call(command, shell = True)

        # 4. Process results such that they are formatted as a list of strings (entez IDs)
        path_to_seeds = f'../temp/{prefix}_domino_seeds.txt'

        result_genes = []
        with open(path_to_seeds[:-4]+"/modules.out", 'r') as results:
            for line in results:
                result_genes.append(line.strip())
        res =  flatten([x[1:-1].split(", ") for x in result_genes])
        # Delete temporary data.
        res = list(set(res))
        mapping_rev = {v: k for k,v in mapping.items()}
        result_genes = [mapping_rev[x] for x in res]
        subprocess.call(f'rm -r ../temp/{path_to_seeds[:-4]}', shell=True)

        subprocess.call(f'rm ../temp/{prefix}_domino_*', shell=True)


        return result_genes, AlgorithmWrapper.mean_degree(ggi_network, result_genes)
