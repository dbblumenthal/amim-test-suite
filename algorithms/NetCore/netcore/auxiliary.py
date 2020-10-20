### Gal Barel
### 4.4.2019
### auxiliary functions

import pandas as pd
import numpy as np
import sys, os

def check_input(restart_prob,normalisaton_method,output_dir,weights_file,seed_file):
	assert restart_prob>0.0 and restart_prob<=1.0 , "restart probability must be between 0 and 1"
	assert normalisaton_method in ['deg_norm','core_norm','deg_core_diff','deg_core_ratio'] , 'normalisation method must be one of the follwing: deg_norm/core_norm/deg_core_diff/deg_core_ratio'
	assert weights_file is not None or seed_file is not None, "Missing input: either weights file or seed file must be provided"
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)


def create_seeds_from_weights(weights_sorted,nodes_order):
	"""
		function to create seed nodes based on input weights
	"""
	
	print("Creating seed nodes based on input weights...")

	# make a df 
	weights_df = pd.DataFrame({"node":nodes_order,"weight":weights_sorted})

	# check if the weights are binary
	weights_count = weights_df['weight'].value_counts()

	if len(weights_count.index)==2:
		# there are only 2 values for the weights (0 and 1)
		if (weights_count.index[0]==0 and weights_count.index[1]==1) or (weights_count.index[0]==1 and weights_count.index[1]==0):
			# the given weights are binary -> seed nodes will be the ones with a weight of 1
			seed_nodes = weights_df.loc[weights_df['weight']==1]['node'].tolist()

	else:
		# sort the df
		sorted_weights_df = weights_df.sort_values('weight',ascending=False)

		# take the top 100 nodes if all of their weights are above 0
		if sorted_weights_df.iloc[100]['weight']>0:
			seed_nodes = sorted_weights_df.iloc[0:100]['node'].tolist()
			### NOTE : even if the dataframe contains less than 100 rows, pandas still returns all the rows avilable from 0 until end ###
		else:
			# take all the nodes with a weight != 0 
			seed_nodes = sorted_weights_df.loc[sorted_weights_df['weight']>0]['node'].tolist()
		

	return seed_nodes


def get_seed_nodes(seed_file,vertex_index):

	seed_list = pd.read_csv(seed_file,sep=r"\s+", header=None , error_bad_lines=False)

	if seed_list.isnull().values.any():
		sys.exit('Seed nodes file contains NA values')

	if seed_list.shape[1]!=1 :
		sys.exit('Seed nodes file should contain only one column with node names')

	seed_list.columns = ["node"]
	
	# get the index of the seed nodes
	seed_list_index = vertex_index.loc[vertex_index['node'].isin(seed_list['node'])]

	# print if there are seed nodes that are not in the network
	if seed_list.shape[0]!=seed_list_index.shape[0]:
		print(str(seed_list.shape[0]-seed_list_index.shape[0]) + " seed nodes are not the network -> They will not be included in the final modules!")

	if seed_list_index.shape[0]==0:
		# none of the seed nodes are in the network! Can't run!
		sys.exit('None of the seed nodes are in the given network -> Cannot run without seed nodes')

	return seed_list_index['index'].tolist()


def get_edges(edge_file,vertex_index=None):
	# open the edge list file
	edges = pd.read_csv(edge_file,sep=r"\s+", header=None)
	is_edge_weights=False

	if edges.isnull().values.any():
		sys.exit('Edges file contains NA values')

	if edges.shape[1]==2:
		print("Network without edge weights")
		edges.columns = ['node1','node2']

	elif edges.shape[1]==3:
		print("Network with edge weights (assumed third column of input file as edge weights)")
		edges.columns = ['node1','node2','weight']
		assert edges.dtypes[2]==int or edges.dtypes[2]==float, "Weights must be numeric"
		is_edge_weights=True

	else:
		sys.exit('Edges file can contain only 2 or 3 columns: node1 , node2 ,(weight)')

	edges = edges.astype("object")
	assert edges.dtypes[0]==object and edges.dtypes[1]==object, "Edge list must contain node names as strings"

	if vertex_index is None:
		# create a vertex index
		node1_unique = edges['node1'].unique()
		node2_unique = edges['node2'].unique()
		all_nodes_unique = list(set(node1_unique).union(set(node2_unique)))
		# list must be sorted because union will return a different order of nodes every time
		all_nodes_unique.sort()
		all_nodes_unique_df = pd.DataFrame({"node":all_nodes_unique,"index":range(1,len(all_nodes_unique)+1)})
	else:
		# vertex index is given
		all_nodes_unique_df = vertex_index

	# create an edge list using the vertex index
	edges["node1"] = edges["node1"].astype(str)
	edges["node2"] = edges["node2"].astype(str)
	all_nodes_unique_df["node"] = all_nodes_unique_df["node"].astype(str)

	node1_in_vertex_index = pd.merge(edges,all_nodes_unique_df,left_on=['node1'], right_on = ['node'], how = 'left')
	node1_vertex_index = node1_in_vertex_index['index'].tolist()

	node2_in_vertex_index = pd.merge(edges,all_nodes_unique_df,left_on=['node2'], right_on = ['node'], how = 'left')
	node2_vertex_index = node2_in_vertex_index['index'].tolist()

	edge_list_vertex_index = pd.DataFrame({"node1":node1_vertex_index,"node2":node2_vertex_index})
	if is_edge_weights:
		# add the weights
		edge_list_vertex_index['weight'] = edges['weight']

	return edge_list_vertex_index , is_edge_weights , all_nodes_unique_df


def get_vertex(vertex_file):
	# open the index vertex file
	vertex_index = pd.read_csv(vertex_file,sep=r"\s+", header=None, error_bad_lines=False)

	if vertex_index.isnull().values.any():
		sys.exit('Vertex file contains NA values')

	if vertex_index.shape[1]!=2 :
		sys.exit('Vertex file shpuld contain only 2 columns: index , node')

	vertex_index.columns = ['index','node']
	vertex_index.loc[:,'index'] = vertex_index['index'].apply(lambda x: int(x))
	vertex_index.loc[:,'node'] = vertex_index['node'].apply(lambda x: str(x))

	return vertex_index


def update_vertex_index(vertex_index,removed_nodes):
	"""
		update the vertex index list after removing some nodes that are not in the biggest connected component in the network
	"""
	vertex_index_removed = vertex_index.loc[~vertex_index['index'].isin(removed_nodes)]
	return vertex_index_removed


def get_weights(weights_file):
	
	# open the weights file 
	node_weights = pd.read_csv(weights_file, sep='\t', header=None)
	
	if node_weights.isnull().values.any():
		sys.exit('Weight file contains NA values')

	if node_weights.shape[1]==1:
		# the file is only a list of genes, give all the genes a weight of 1
		print("No input weights given -> all nodes in the weights file will have a weight of 1")
		node_weights['weight'] = 1

	elif node_weights.shape[1]==0 or node_weights.shape[1]>2:
		sys.exit('Weight file should contain only 2 columns: node , weight (TAB separated only)')
	
	node_weights.columns = ['node','weight']		
	node_weights['node'] = node_weights['node'].astype(str)
	try:
		node_weights['weight'] = node_weights['weight'].astype(float)
	except ValueError:
		print('Weight file -> second colum must contain ints')
		raise

	return node_weights


def process_weights(node_weights,vertex_index,G,nodes_order):
	# process the input weights and give weight 0 to any node in the network that doesn't have an input weight

	print("Processing weights...")

	# merge the weights df with the index df
	vertex_index["node"] = vertex_index["node"].astype(str)
	node_weights_index = pd.merge(node_weights, vertex_index, on=['node'], how='outer')

	# remove nodes that have a weight but not an index in the network
	if node_weights_index['index'].isnull().values.any():

		nodes_with_weight_but_no_index = node_weights_index.loc[node_weights_index['index'].isnull(),]
		print(str(nodes_with_weight_but_no_index.shape[0]) + " nodes with an input weight are not in the network -> They will be removed from further analysis!")
		node_weights_index.drop(nodes_with_weight_but_no_index.index.tolist(),inplace=True)

	# give weight 0 to nodes in the network that don't have an input weight
	if node_weights_index['weight'].isnull().values.any():

		nodes_with_index_but_no_weight = node_weights_index.loc[node_weights_index['weight'].isnull(),]
		print(str(nodes_with_index_but_no_weight.shape[0]) + " nodes in the network don't have an input weight -> They will be assigned a weight of 0!")
		node_weights_index.loc[node_weights_index['weight'].isnull(),'weight'] = 0.0

	# convert index to int
	node_weights_index.loc[:,'index'] = node_weights_index['index'].apply(lambda x : int(x))

	# create a node and weight list
	nodes = node_weights_index.iloc[:]['index'].tolist()
	weights = node_weights_index.iloc[:]['weight'].tolist()

	# create a dict for the weights
	node_weight_dict = dict(zip(nodes,weights))

	# order the weights of the nodes according to the nodes order
	weights_sorted = [node_weight_dict[node] for node in nodes_order]

	# create a diagonal matrix of the weights
	weights_mat = np.diag(weights_sorted)

	return weights_sorted , weights_mat


def write_edgelist_with_node_names(G,vertex_index_dict,is_edge_weights,edgelist_file_name):
	"""
		function to save the edgelist such that the names of the nodes are used instead of 
	"""

	with open(edgelist_file_name, 'w') as f:

		if is_edge_weights:
			# write edges with weight data
			for u, v, d in G.edges(data=True):
				f.write("{0}\t{1}\t{2}\n".format(vertex_index_dict[u],vertex_index_dict[v],d['weight']))

		else:
			for u, v in G.edges(data=False):
				f.write("{0}\t{1}\n".format(vertex_index_dict[u],vertex_index_dict[v]))

	f.close()


def write_adjlist_with_node_names(G,vertex_index_dict,adjlist_file_name):
	"""
		This function will write a network in the format of an adjecency list 
		The nodes will be saved using their names (instead of index)
		In order to open the file using networkx:
		nx.read_adjlist(file_name, create_using = nx.Graph(), nodetype = str)
	"""

	delimiter="\t"
	with open(adjlist_file_name, 'w') as f:
		for s, nbrs in G.adjacency():
		    line = vertex_index_dict[s] + delimiter
		    for t, data in nbrs.items():
		        line += vertex_index_dict[t] + delimiter
		    f.write(line[:-len(delimiter)]+"\n")
	f.close()
