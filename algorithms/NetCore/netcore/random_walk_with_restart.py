### Gal Barel
### 4.4.2019
### This file contains the functions for running a random walk with restart

import networkx as nx
import numpy as np
import pandas as pd
import scipy

from auxiliary import get_edges , update_vertex_index

def edges_to_net(edge_list,vertex_index=None):
	"""read the edges file and convert into a nx network"""

	# check the edge list file
	edge_list_df , is_edge_weights, vertex_index = get_edges(edge_list,vertex_index)

	# create network from edge list file
	if is_edge_weights:
		G=nx.from_pandas_edgelist(edge_list_df,'node1','node2',['weight'])
		# when the weights are floats then the node names are converted to floats as well
		# need to convert back from float to int 
		mapping = dict(zip(G,[int(n) for n in G.nodes()]))
		G = nx.relabel_nodes(G, mapping)
	else:
		G=nx.from_pandas_edgelist(edge_list_df,'node1','node2')
		# in this case the nodes are read as int 

	# remove self edges (loops)
	G.remove_edges_from(nx.selfloop_edges(G))

	# remove nodes that are not connected to the biggest SCC
	G_connected = max(nx.connected_component_subgraphs(G), key=len)

	# get node that were removed from the network
	removed_nodes = [n for n in G.nodes() if n not in G_connected.nodes()]

	# update the vertex_index list
	vertex_index = update_vertex_index(vertex_index,removed_nodes)
	
	if len(removed_nodes)>0:
		print(str(len(removed_nodes)) + " nodes were removed from the network since they are not connected to the biggest component")

	return G_connected, removed_nodes, is_edge_weights , vertex_index


def make_network(edge_list):

	print("Getting network information...")

	G_connected, removed_nodes, is_edge_weights , vertex_index = edges_to_net(edge_list)

	# order the nodes
	nodes_order = sorted(G_connected.nodes())

	# get the adjacency matrix of the graph, ordered accodring to the nodes
	adj_mat = nx.adjacency_matrix(G_connected,nodelist=nodes_order)
	adj_mat_dense = adj_mat.todense()

	# get the node degree and k-shell
	node_deg = G_connected.degree()
	node_deg_order = np.array([node_deg[n] for n in nodes_order])
	node_core = nx.core_number(G_connected)
	node_core_order = np.array([node_core[n] for n in nodes_order])

	return vertex_index, G_connected, adj_mat_dense , nodes_order , node_deg_order , node_core_order, removed_nodes, is_edge_weights


def get_random_walk_matrix(normalisaton_method,is_edge_weights,adj_mat_dense,node_deg_order,node_core_order,restart_prob,output_dir):

	print("Creating random walk matrix...")

	# normalize the adj matrix by the degree
	if normalisaton_method=='deg_norm':
		norm_adj_mat = norm_by_degree(adj_mat_dense,node_deg_order,is_edge_weights)
		subnets_file = output_dir+"deg_norm_subnetworks.txt"

	if normalisaton_method=='core_norm':
		norm_adj_mat = norm_by_core(adj_mat_dense,node_core_order,is_edge_weights)
		subnets_file = output_dir+"core_norm_subnetworks.txt"

	# normalize by degree-core difference
	if normalisaton_method=='deg_core_diff':
		norm_adj_mat = norm_by_deg_core_diff(adj_mat_dense,node_deg_order,node_core_order)
		subnets_file = output_dir+"deg_core_diff_subnetworks.txt"

	if normalisaton_method=="deg_core_ratio":
		norm_adj_mat = norm_by_deg_core_ratio(adj_mat_dense,node_deg_order,node_core_order)
		subnets_file = output_dir+"deg_core_ratio_subnetworks.txt"
	
	# get the random walk matrix
	I = np.matrix(np.identity(norm_adj_mat.shape[0]))
	random_walk_mat = restart_prob*scipy.linalg.fractional_matrix_power(I-((1-restart_prob)*norm_adj_mat),-1.0)

	return norm_adj_mat , random_walk_mat, subnets_file


def norm_by_degree(adj_mat,node_deg_order,is_edge_weights):
	""" given the adj matrix and degree matrix normalise the adj matrix by the degree:
		W = AD^-1
	"""
	
	# get the diagonal matrix of the degrees
	node_deg_order_rev = 1.0 / node_deg_order # 1/D = D^-1
	# if degree is 0, rev value is inf, so change back to zero
	node_deg_order_rev[node_deg_order_rev == np.inf] = 0.0
	deg_mat = np.diag(node_deg_order_rev)	
	
	# multiple the adj matrix by the reverse degree diagonal matrix
	deg_norm_adj_mat = adj_mat.dot(deg_mat)

	if is_edge_weights:
		# need to normalize the columns such that they sum up to 1
		col_sum = deg_norm_adj_mat.sum(axis=0)
		col_sum[col_sum==0.0] = -1.0 #nodes with degree 0 will be set to sum of -1 to avoid division by 0

		deg_norm_adj_mat = deg_norm_adj_mat.astype(float) /  col_sum
		deg_norm_adj_mat[deg_norm_adj_mat<0.0] = 0.0 #replace the negative values back to 0

	return deg_norm_adj_mat


def norm_by_core(adj_mat,k_shell,is_edge_weights):

	""" given an adjacncy matrix 
		construct a k_shell matrix
		where col i corresponds to node i
		M(i,j) = kj if i and j are neighbours
		where kj is the k-shell of node j
		@ adj_mat is of type matrix
		@ k_shell is a list of the same order as the rows and cols in adj_mat
	"""
	
	k_shell_mat = adj_mat.copy()
	for col in range(adj_mat.shape[1]):
		# get the index of the neighbours
		node_col = adj_mat[:,col]
		neighbours = np.where(node_col>0)[0]
		neighbours_k_shell = k_shell[neighbours]
		if is_edge_weights:
			# need to multiply by the edge weights
			k_shell_mat[neighbours,col] = list(k_shell_mat[neighbours,col].flat)*neighbours_k_shell
		else:
			k_shell_mat[neighbours,col] = neighbours_k_shell

	# normalise each col to make a stochiometric matrix
	col_sum = k_shell_mat.sum(axis=0)
	col_sum[col_sum==0.0] = -1.0 #nodes with degree 0 will be set to sum of -1 to avoid division by 0

	k_shell_mat = k_shell_mat.astype(float) /  col_sum
	k_shell_mat[k_shell_mat<0.0] = 0.0 #replace the negative values back to 0

	return k_shell_mat


def norm_by_deg_core_diff(adj_mat,node_deg_order,k_shell):

	""" given an adjacncy matrix, the nodes degree and k-shell values
		normalize the adj matrix by both the degree of the node and by the k-shell of its neighbours

	"""

	diff_norm_mat = adj_mat.copy().astype(float)

	for col in range(diff_norm_mat.shape[1]):
		
		# get the index of the neighbours
		node_col = diff_norm_mat[:,col]
		neighbours = np.where(node_col>0)[0]

		# get the degree and k-shell of the neighbours 
		neighbours_deg = node_deg_order[neighbours]
		neighbours_k_shell = k_shell[neighbours]

		# normalize by the difference of the degree and the kshell
		deg_kshell_diff = (neighbours_deg - neighbours_k_shell) + 1 #add 1 to avoid division by 0
		diff_norm_mat[neighbours,col] = np.multiply(np.transpose(diff_norm_mat[neighbours,col]),1.0/deg_kshell_diff)
		
	# normalise each col to make a stochiometric matrix
	col_sum = diff_norm_mat.sum(axis=0)
	col_sum[col_sum==0.0] = -1.0 #nodes with degree 0 will be set to sum of -1 to avoid division by 0

	diff_norm_mat = diff_norm_mat.astype(float) /  col_sum
	diff_norm_mat[diff_norm_mat<0.0] = 0.0 #replace the negative values back to 0

	return diff_norm_mat


def norm_by_deg_core_ratio(adj_mat,node_deg_order,k_shell):

	""" given an adjacncy matrix, the nodes degree and k-shell values
		normalize the adj matrix by both the degree of the node and by the k-shell of its neighbours

	"""

	ratio_norm_mat = adj_mat.copy().astype(float)

	for col in range(ratio_norm_mat.shape[1]):
		
		# get the index of the neighbours
		node_col = ratio_norm_mat[:,col]
		neighbours = np.where(node_col>0)[0]

		# get the degree and k-shell of the neighbours 
		neighbours_deg = node_deg_order[neighbours]
		neighbours_k_shell = k_shell[neighbours]

		# normalize by the ratio of the degree and the kshell
		ratio = np.array(neighbours_k_shell,dtype=float)/np.array(neighbours_deg,dtype=float)
		ratio_norm_mat[neighbours,col] = np.multiply(np.transpose(ratio_norm_mat[neighbours,col]),ratio)

		
	# normalise each col to make a stochiometric matrix
	col_sum = ratio_norm_mat.sum(axis=0)
	col_sum[col_sum==0.0] = -1.0 #nodes with degree 0 will be set to sum of -1 to avoid division by 0

	ratio_norm_mat = ratio_norm_mat.astype(float) /  col_sum
	ratio_norm_mat[ratio_norm_mat<0.0] = 0.0 #replace the negative values back to 0

	return ratio_norm_mat