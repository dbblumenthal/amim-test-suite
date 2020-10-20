### Gal Barel
### 29.05.2019
### permutation test script

# import for swapping implementation
from __future__ import division
import math
import random

# impot 
from random_walk_with_restart import *
from auxiliary import write_edgelist_with_node_names
import os,copy
import multiprocessing as mp


def make_network_permutations(net_file,net_name,output_path,num_perm=100,swap_factor=100,num_cores=1):
	"""
		generate degree preserving network permutations using edge swap alg.
		can be done in a parallel fashion 
	"""

	# get network
	print("Getting network")
	G, removed_nodes, is_edge_weights, vertex_index = edges_to_net(net_file)

	# create output directorty
	output_dir = output_path + net_name + "_edge_permutations/"
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)

	# set number of swaps
	num_swap = G.number_of_edges()*swap_factor

	if num_cores>1:
		# run in parallel
		run_edge_swap_parallel(net_name,G,vertex_index,is_edge_weights,num_swap,num_perm,num_cores,output_dir)
	else:
		# run not parallel
		for i in range(num_perm):
			swap_edges([G,vertex_index,is_edge_weights,num_swap,i,output_dir,net_name])


def get_permutation_weights(perm_dir,perm_name,net_name,num_perm,norm_adj_mat,vertex_index,nodes_order,node_deg_order,node_core_order,weights_sorted,normalisaton_method,is_edge_weights,alpha):
	"""
	get the permutated weights after the random walk
	return a data frame where each colum are the weights using a permuted network
	"""

	print("Running permutations test...")

	# create random walk matrix for each permutation

	I = np.matrix(np.identity(norm_adj_mat.shape[0]))

	# create a df that will contain all of the propagation weights
	perm_weights = pd.DataFrame({'node':nodes_order})

	for i in range(0,num_perm+1):
		print("Permutation " + str(i))
		
		# get the permuted network

		if net_name in perm_name:
			# don't need to add the network name to the permutation file name
			perm_file = perm_dir + perm_name + str(i)
		else:
			perm_file = perm_dir+ net_name + "_" + perm_name + str(i)

		if os.path.isfile(perm_file):
			perm_net , removed_nodes, edge_weights, perm_vertex_index = edges_to_net(perm_file,vertex_index)
		else:
			# file not found, try with a txt extention
			perm_file = perm_file + ".txt"
			if os.path.isfile(perm_file):
				perm_net , removed_nodes, edge_weights, perm_vertex_index = edges_to_net(perm_file,vertex_index)
			else:
				print("No permutation file found: " + perm_file)
				continue

		# get the adjacency matrix of the network, ordered accodring to the nodes
		perm_adj_mat = nx.adjacency_matrix(perm_net,nodelist=nodes_order)
		perm_adj_mat_dense = perm_adj_mat.todense()

		perm_norm_adj_mat, perm_norm_random_walk_mat, subnets_file= get_random_walk_matrix(normalisaton_method,is_edge_weights,perm_adj_mat_dense,node_deg_order,node_core_order,alpha,perm_dir)
		
		# get the network propagation weights
		perm_norm_random_walk_weights = perm_norm_random_walk_mat.dot(weights_sorted)
		
		# add weights to the big df
		perm_weights[str(i)] = perm_norm_random_walk_weights

	return perm_weights


def get_permutation_pvals(perm_weights,random_walk_weights,nodes_order):
	"""
	get the p-values based on the permutation weights
	"""
	# create return vector 

	print("Getting permutation pvalues...")

	p_vals = [0]*len(nodes_order)

	# the first colum of the dataframe is the node , remove it
	perm_weights_only = perm_weights.drop("node", axis=1, inplace=False)

	for i in range(len(nodes_order)):
		# count how many times the random weight was larger than the "real" one
		num_larger = np.where(perm_weights_only.loc[i]>=random_walk_weights[i])
		# the pval is the number of times the random weights are larger, out of all permutations
		pval = float(len(num_larger[0])+1) / float(perm_weights_only.shape[1]+1)
		p_vals[i] = pval

	return p_vals



def swap_edges(args):
	"""
		use the networkx alg to swap edges and save the permuted network
	"""
	
	# get arguments
	G,vertex_index,is_edge_weights,num_swap,index,swap_dir,net_name = args
	
	print("Running permutation " + str(index))

	# copy network before applying swapp alg
	G_swap = copy.deepcopy(G)
	connected_double_edge_swap_with_weights(G_swap, is_edge_weights, nswap=num_swap, seed=random)

	# get the nodes names to save in the edge list file (instead of the index)
	vertex_index_dict = dict(zip(vertex_index['index'], vertex_index['node']))
	
	# save the network into an edgelist file using the node names
	edgelist_file_name = swap_dir + net_name + "_edgelist"+str(index)+".txt"
	write_edgelist_with_node_names(G_swap,vertex_index_dict,is_edge_weights,edgelist_file_name)


def run_edge_swap_parallel(net_name,net,vertex_index,is_edge_weights,num_swap,num_perm,num_cores,output_dir):
	"""
		run in parallel edge swapping
	"""

	print("Making Edge permutations (in parallel)")
	# run parallel
	with mp.Pool(processes = num_cores) as p:
		results = p.map(swap_edges, [(net,vertex_index,is_edge_weights,num_swap,i,output_dir,net_name) for i in range(num_perm)])


def connected_double_edge_swap_with_weights(G, is_edge_weights, nswap=1, _window_threshold=3, seed=None):
	"""
		Added edge weights in the original implementation of the alg. from networkx.
		The weight is randomly assigned after the swapp to the new pairs of edges
	"""
	if not nx.is_connected(G):
		raise nx.NetworkXError("Graph not connected")
	if len(G) < 4:
		raise nx.NetworkXError("Graph has less than four nodes.")
	n = 0
	swapcount = 0
	deg = G.degree()
	# Label key for nodes
	dk = list(n for n, d in G.degree())
	cdf = nx.utils.cumulative_distribution(list(d for n, d in G.degree()))
	discrete_sequence = nx.utils.discrete_sequence
	window = 1
	while n < nswap:
		wcount = 0
		swapped = []
		# If the window is small, we just check each time whether the graph is
		# connected by checking if the nodes that were just separated are still
		# connected.
		if window < _window_threshold:
			# This Boolean keeps track of whether there was a failure or not.
			fail = False
			while wcount < window and n < nswap:
				# Pick two random edges without creating the edge list. Choose
				# source nodes from the discrete degree distribution.
				(ui, xi) = discrete_sequence(2, cdistribution=cdf)
				# If the source nodes are the same, skip this pair.
				if ui == xi:
					continue
				# Convert an index to a node label.
				u = dk[ui]
				x = dk[xi]
				# Choose targets uniformly from neighbors.
				v = seed.choice(list(G.neighbors(u)))
				y = seed.choice(list(G.neighbors(x)))
				# If the target nodes are the same, skip this pair.
				if v == y:
					continue
				if x not in G[u] and y not in G[v]:
					if is_edge_weights:
						G.add_edge(u, x, weight=G[u][v]['weight'])
						G.add_edge(v, y, weight=G[x][y]['weight'])
					else:
						G.add_edge(u, x)
						G.add_edge(v, y)
					G.remove_edge(u, v)
					G.remove_edge(x, y)
					swapped.append((u, v, x, y))
					swapcount += 1
				n += 1
				# If G remains connected...
				if nx.has_path(G, u, v):
					wcount += 1
				# Otherwise, undo the changes.
				else:
					if is_edge_weights:
						G.add_edge(u, v, weight=G[u][x]['weight'])
						G.add_edge(x, y, weight=G[v][y]['weight'])
					else:
						G.add_edge(u, v)
						G.add_edge(x, y)
					G.remove_edge(u, x)
					G.remove_edge(v, y)
					swapcount -= 1
					fail = True
			# If one of the swaps failed, reduce the window size.
			if fail:
				window = int(math.ceil(window / 2))
			else:
				window += 1
		# If the window is large, then there is a good chance that a bunch of
		# swaps will work. It's quicker to do all those swaps first and then
		# check if the graph remains connected.
		else:
			while wcount < window and n < nswap:
				# Pick two random edges without creating the edge list. Choose
				# source nodes from the discrete degree distribution.
				(ui, xi) = nx.utils.discrete_sequence(2, cdistribution=cdf)
				# If the source nodes are the same, skip this pair.
				if ui == xi:
					continue
				# Convert an index to a node label.
				u = dk[ui]
				x = dk[xi]
				# Choose targets uniformly from neighbors.
				v = seed.choice(list(G.neighbors(u)))
				y = seed.choice(list(G.neighbors(x)))
				# If the target nodes are the same, skip this pair.
				if v == y:
					continue
				if x not in G[u] and y not in G[v]:
					if is_edge_weights:
						G.add_edge(u, x, weight=G[u][v]['weight'])
						G.add_edge(v, y, weight=G[x][y]['weight'])
					else:
						G.add_edge(u, x)
						G.add_edge(v, y)
					G.remove_edge(u, v)
					G.remove_edge(x, y)
					swapped.append((u, v, x, y))
					swapcount += 1
				n += 1
				wcount += 1
			# If the graph remains connected, increase the window size.
			if nx.is_connected(G):
				window += 1
			# Otherwise, undo the changes from the previous window and decrease
			# the window size.
			else:
				while swapped:
					(u, v, x, y) = swapped.pop()
					if is_edge_weights:
						G.add_edge(u, v, weight=G[u][x]['weight'])
						G.add_edge(x, y, weight=G[v][y]['weight'])
					else:
						G.add_edge(u, v)
						G.add_edge(x, y)
					G.remove_edge(u, x)
					G.remove_edge(v, y)
					swapcount -= 1
				window = int(math.ceil(window / 2))
	return swapcount