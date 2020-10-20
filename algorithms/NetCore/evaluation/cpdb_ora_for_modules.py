### Gal Barel
### 8.11.2019
### CPDB ORA

### Note : this script needs to run using Python2 (beacuse of some libraries of CPDB_WSDL)

import sys , os
import pandas as pd
from argparse import ArgumentParser


sys.path.append('/project/Toxicogenomics/Toxicogenomics/CPDB_WSDL/')
from cpdb_services import *

# Parse arguments.
def get_parser():
	description = 'Run CPDB ORA'
	parser = ArgumentParser(description=description)
	parser.add_argument('-op', '--output_path', type=str, required=True, help='Path to save the results')
	parser.add_argument('-mf', '--modules_file', type=str, required=True, help='File that includes all the modules')
	parser.add_argument('-bg', '--background_file', type=str, required=True, help='File that includes a list of nodes that will be the background for the enrichment analysis')
	parser.add_argument('-sf', '--seed_file', type=str, default=None, required=False, help='File that includes a list of seed nodes')
	parser.add_argument('-ms', '--min_size', type=int, default=10, required=False, help='Minimum size of module for doing the enrichment analysis')
	return parser

def run_ora_modules(ora_path,modules_file,background_file,seed_file,min_size):
	"""
		function to open the modules file and run CPDB ORA for each module above a given size
	"""

	print("Getting CPDB ORA for modules...")

	# open the file of the all the genes and use as background genes
	background_df = pd.read_csv(background_file, sep="\t", header=None)
	back_genes = background_df[1].tolist()
	### NOTE: assuming the file is in the same format as the vertex index file of the network ###

	# open file and read line by line
	with open(modules_file) as f:
		line_count=0
		for line in f:
			if line.startswith("#"):
				# this is a comment line - > skip !
				continue
			if ',' in line:
				# this is a csv file , need to get the modules from the first part
				module_genes_part = line.split(',')[0]
				module_genes = module_genes_part.split("\t")
				if len(module_genes)==1:
					# this is a header line
					continue
			else:
				# all of the line are the genes
				module_genes = line.split("\t")
			
			if len(module_genes)>min_size:
				print("Module " + str(line_count))
				ora_paths = run_ora(module_genes,back_genes)
				ora_paths.drop('details', axis=1, inplace=True) # remove the details col
				c_file = os.path.join(ora_path,"module" + str(line_count) + "_pathways.txt")
				ora_paths.to_csv(c_file, sep='\t',index=False)

				# if seed nodes are given -> check results without seed nodes and of seed nodes only
				if seed_file is not None:

					print("Getting ORA for seed and non-seed nodes...")

					# get the seed nodes
					seed_genes_df = pd.read_csv(seed_file, sep="\t", header=None)
					seed_genes = seed_genes_df[0].tolist()

					# get the module genes that are not seed nodes
					module_genes_seed = [g for g in module_genes if g in seed_genes]
					module_genes_non_seed = [g for g in module_genes if g not in seed_genes]

					# run ORA using genes from seed list only 
					if len(module_genes_seed)>min_size:
						ora_path_seed = run_ora(module_genes_seed,back_genes)
						ora_path_seed.drop('details', axis=1, inplace=True)
						ora_path_seed_file = os.path.join(ora_path,"module" + str(line_count) + "_seed_only_pathways.txt")
						ora_path_seed.to_csv(ora_path_seed_file, sep='\t',index=False)

					# run ORA using the genes that are NOT in the seed list
					if len(module_genes_non_seed)>min_size:
						ora_path_non_seed = run_ora(module_genes_non_seed,back_genes)
						ora_path_non_seed.drop('details', axis=1, inplace=True)
						ora_path_non_seed_file = os.path.join(ora_path,"module" + str(line_count) + "_non_seed_pathways.txt")
						ora_path_non_seed.to_csv(ora_path_non_seed_file, sep='\t',index=False)

			# update module count
			line_count+=1

	print("###DONE###")


"""
Function to run ORA via CPDB
input : - gene_names : a list of gene names to run the analysis for  (gene_names need to be hgnc-symbols)
		- background_genes : a list of background genes to use as a background for the ORA
output : data fame with over represented pathways
"""
def run_ora(gene_names,background_genes):

	loc = cpdbLocator()
	proxy = loc.getcpdb_portType(tracefile=sys.stdout)

	# the very first connection fails so do it twice
	req = getCpdbVersionRequest()
	try:
		response = proxy.getCpdbVersion(req)
	except:
		print("First connection failed... trying again")
	try:
		response = proxy.getCpdbVersion(req)
	except:
		print("Cannot Conncet")
		return

	print("connection to " + str(response._cpdbVersion) + " is successful!")

	## first map a list of gene_names and background_genes to cpdbIds
	gene_names_cpdbIds = map_gene_names_to_cpdb_ids(gene_names,proxy) #TODO : move this outside the function so that it's only done once
	background_genes_cpdbIds = map_gene_names_to_cpdb_ids(background_genes,proxy)

	req = overRepresentationAnalysisRequest()
	req._entityType = 'genes'
	req._fsetType = 'P' # request for pathways
	req._cpdbIdsFg = gene_names_cpdbIds
	req._cpdbIdsBg = background_genes_cpdbIds # need to set up the background IDs and then it will return all pathways
	req._pThreshold = 0.01
	response = proxy.overRepresentationAnalysis(req)
	result = zip(response._name, response._details, response._overlappingEntitiesNum, response._allEntitiesNum, response._pValue, response._qValue)
	results_df = pd.DataFrame(result,columns=['path_name','details','overlapp','size','pvalue','qvalue']) #save the results as a dataframe
	return results_df

"""
Function to map gene names (HGNC symbols) to CPDB IDs
"""
def map_gene_names_to_cpdb_ids(gene_names,proxy):

	req = mapAccessionNumbersRequest()
	req._accType = 'hgnc-symbol'
	req._accNumbers = gene_names
	response = proxy.mapAccessionNumbers(req)
	mapping = zip(response._accNumber, response._cpdbId)
	cpdbIds = []
	for r in mapping:
		if r[1]:
			cpdbIds.append(r[1].split(',')[0])  ## here we take only one of the mappings for simplicity

	return cpdbIds


def run(args):

	run_ora_modules(args.output_path,args.modules_file,args.background_file,args.seed_file,args.min_size)


if __name__ == "__main__":
	run(get_parser().parse_args(sys.argv[1:]))