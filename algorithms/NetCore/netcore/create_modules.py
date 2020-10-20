### Gal Barel
### 4.4.2019
### This file contains the functions to find network modules

import os 
import numpy as np
import networkx as nx
import itertools

from auxiliary import write_adjlist_with_node_names

# set graphics
from matplotlib import pyplot
import seaborn as sns
sns.set(style="whitegrid")


def get_seed_modules(network,seed_nodes,random_walk_weights_pvals,pval_thresh,weight_thresh,max_subnet_size,vertex_index,output_file,output_dir):
	"""
		function to get the modules based on the seed nodes
			- start from the subnetwork that is induced by the seed nodes
			- add gradually nodes with a significant pval and weight
			- return modules in the subnetwork (connected components)
	"""

	print("Getting seed modules...")

	# get the subnetwork that is induced by the seed nodes
	seed_subnet = network.subgraph(seed_nodes)

	# sort the propagation results by pvalue and then weight
	sorted_random_walk_weights_pvals = random_walk_weights_pvals.sort_values(['pvalue', 'prop_weight'], ascending=[True, False])

	# get the nodes to add to the subnet
	nodes_to_add , weight_thresh = add_nodes_to_seed_subnet(network,seed_subnet,sorted_random_walk_weights_pvals,pval_thresh,weight_thresh,max_subnet_size)

	if len(nodes_to_add)==0:
		print("No nodes were found that could be added to the seed subnetwork!")

	# create a subnet with both seed nodes and added nodes
	seed_subnet_with_added_nodes = network.subgraph(seed_nodes+nodes_to_add)

	# get the connected components ( at least 2 nodes)
	seed_subnet_comp = [cc for cc in nx.connected_component_subgraphs(seed_subnet_with_added_nodes) if cc.number_of_nodes()>1]

	# plot the subnetwork
	subnet_file=output_dir+"extended_seed_subnet"
#	plot_subet(seed_subnet_with_added_nodes,seed_nodes,nodes_to_add,random_walk_weights_pvals,pval_thresh,weight_thresh,subnet_file+".pdf")

	# save the subnetwork
	# write the edge list with node names
	# get the nodes names to save in the edge list file (instead of the index)
	vertex_index_dict = dict(zip(vertex_index['index'], vertex_index['node']))
	write_adjlist_with_node_names(seed_subnet_with_added_nodes,vertex_index_dict,subnet_file+".adjlist")
	# write_edgelist_with_node_names(seed_subnet_with_added_nodes,vertex_index_dict,False,subnet_file+".edgelist")

	# save the modules (will also plot each module)
	save_modules(seed_subnet_comp,seed_nodes,random_walk_weights_pvals,vertex_index,output_file,output_dir)


def add_nodes_to_seed_subnet(G,seed_subnet,sorted_results,pval_thresh,weight_thresh,max_subnet_size):

	'''
		gradually add nodes to seed subnetwork based on pvalue and propagation weight
	'''

	# get the size of the seed subnet
	seed_subnet_size=seed_subnet.number_of_nodes()

	if max_subnet_size is not None:
		if max_subnet_size<seed_subnet_size:
		# the requested maximum subnet size is smaller than just the seed subnet size
			print("Maximum subnet size requested is smaller than the size of the seed induced subnetwork alone -> calculating the maximum size based on the seed nodes...")
			max_subnet_size=None

	# set the maximum subnet size (unless given)
	if max_subnet_size is None:
		# allow to incease the size of the subnet by 50%
		max_subnet_size = seed_subnet_size + round(seed_subnet_size/float(2))
		print("The maximum subnet size is " + str(max_subnet_size))

	nodes_to_add = []

	# only add nodes that have a pvalue lower than threshold
	sig_sorted_results = sorted_results.loc[sorted_results['pvalue']<=pval_thresh]

	if sig_sorted_results.shape[0]==0:
		# There are no nodes with a significan pvalue
		print("There are no nodes with a pvalue <=" + str(pval_thresh))
		return nodes_to_add , 0.0

	# based on the distribution of the significant nodes, get the minimum weight threshold
	if weight_thresh is None:
		# the minimum weight is calculated based on the propagation weights of the significant nodes that are not from the seed nodes
		sig_sorted_results_non_seed = sig_sorted_results.loc[~sig_sorted_results['node_index'].isin(list(seed_subnet.nodes()))]
		weight_thresh = np.percentile(sig_sorted_results_non_seed.prop_weight,75)


	for index, row in sig_sorted_results.iterrows():
		
		if len(nodes_to_add)>=(max_subnet_size-seed_subnet_size):
			# stop when the number of nodes to add reached the limit
			print("Number of seed induced subnetwork with propagation nodes has reached the maximum size of " +str(max_subnet_size) + "-> no more nodes will be added to the subnetwork!")
			break

		if row['prop_weight']>weight_thresh:
			# consider only nodes with a sig pvalue and high propagation weight
			# check if the node has neighbours in the subnetwork of the seed nodes
			node = row['node_index']
			if node not in seed_subnet.nodes():
				node_neighbours = G.neighbors(node)
				neighbours_in_seed_subnet = set(node_neighbours) & set(seed_subnet.nodes())
				if len(neighbours_in_seed_subnet)>0:
					# add node to seed subnet only if it has at least 1 neighbour there
					nodes_to_add.append(node)

	return nodes_to_add , weight_thresh


def change_range(weights,NewMin=1,NewMax=100):
	"""
		change the range of the weights to be between 1 and 100 (for plotting size)
	"""
	OldMax = max(weights)
	OldMin = min(weights)
	OldRange = (OldMax - OldMin)  
	NewRange = (NewMax - NewMin)
	if OldRange==0:
		# there is only one value for the seed weights
		OldRange=NewRange
	sizes = [0]*len(weights)
	for i,w in enumerate(weights):   
		NewValue = (((w - OldMin) * NewRange) / OldRange) + NewMin
		sizes[i] = NewValue
	return sizes


def plot_subet(subnet,seed_nodes,nodes_to_add,random_walk_weights_pvals,pval,weight,subnet_file):

	# plot the network with the significant genes
	pos = nx.spring_layout(subnet,k=5/np.sqrt(subnet.number_of_nodes()),iterations=50)

	if len(nodes_to_add)>0:
		# get the weights of the nodes that should be added
		nodes_to_add_weights_df = random_walk_weights_pvals.loc[random_walk_weights_pvals['node_index'].isin(nodes_to_add)]
		#sort the weights accoding to the order they will by plotted in
		nodes_to_add_weights_df = nodes_to_add_weights_df.set_index('node_index')
		nodes_to_add_df_sorted = nodes_to_add_weights_df.loc[nodes_to_add]
		nodes_to_add_weights = nodes_to_add_df_sorted['prop_weight'].tolist()

	# get the weights of the seed nodes
	seed_weights_df = random_walk_weights_pvals.loc[random_walk_weights_pvals['node_index'].isin(seed_nodes)]
	#sort the weights accoding to the order they will by plotted in
	seed_weights_df = seed_weights_df.set_index('node_index')
	seed_weights_df_sorted = seed_weights_df.loc[seed_nodes]
	seed_weights = seed_weights_df_sorted['prop_weight'].tolist()

	# calculate the nodes sizes based on their weights
	all_nodes_weights = change_range(seed_weights + nodes_to_add_weights)
	# seed nodes are first and then nodes_to_add
	seed_sizes = all_nodes_weights[:len(seed_weights)]
	nodes_to_add_sizes = all_nodes_weights[len(seed_weights)+1:]


	# add the original seed nodes that are not significant
	nc = nx.draw_networkx_nodes(subnet,pos,nodelist=seed_nodes,node_color="darkorange",node_shape = "o",node_size=seed_sizes)

	# add the rest of the nodes
	if len(nodes_to_add)>0:
		nc = nx.draw_networkx_nodes(subnet,pos,nodelist=nodes_to_add,node_color="darkgray", node_shape = "o",node_size=nodes_to_add_sizes)

	# add edges
	ec = nx.draw_networkx_edges(subnet,pos,restart_prob=0.5,edge_color="lightgray")

	pyplot.title("Extended seed subnetwork with p<" + str(pval) + " and w>" + str(round(weight,3)) ,fontsize=16)
	pyplot.axis('off')
	pyplot.savefig(subnet_file)
	pyplot.close()


def plot_module(module,module_file,seed_nodes,vertex_index):
	"""
	function to plot a network module , with a different shape for seed nodes 

	"""
	# get module and seed nodes
	nodes = module.nodes()
	non_seed_nodes = [n for n in nodes if n not in seed_nodes]

	# set plotting params 
	nodes_labels = vertex_index.loc[vertex_index['index'].isin(nodes)]
	nodes_labels_dict = dict(zip(nodes_labels['index'],nodes_labels['node']))
	
	# set layout
	pos = nx.spring_layout(module,k=10/np.sqrt(module.number_of_nodes()),iterations=50)

	node_size=200
	if module.number_of_nodes()>300:
		node_size = (module.number_of_nodes()/2)

	# add the original seed noes
	nc = nx.draw_networkx_nodes(module,pos,nodelist=seed_nodes,node_color="darkorange",node_size=node_size)
	
	# add the rest of the nodes
	nc = nx.draw_networkx_nodes(module,pos,nodelist=non_seed_nodes,node_color="darkgray",node_size=node_size)
	
	# add node label
	nx.draw_networkx_labels(module, pos, nodes_labels_dict,font_size=6,font_color='black')
	
	# add edges
	ec = nx.draw_networkx_edges(module, pos, restart_prob=0.5, edge_color="lightgray")
	
	pyplot.axis('off')
	pyplot.savefig(module_file)
	pyplot.close()


def save_modules(modules,seed_nodes,random_walk_weights_pvals,vertex_index,output_file,output_dir):
	"""
		save a file with all the nodes in each module
		plot the modules, mark seed nodes if avilable
		@ modules - a list of subnetworks   
	"""

	print("Saving modules...")

	# create directory for plots
	modules_plot_dir = output_dir+"/modules"
	if not os.path.exists(modules_plot_dir):
		os.makedirs(modules_plot_dir)

	with open(output_file,'w') as of:
		# write header for output file 
		of.write('nodes'+','+'size'+','+'sum_weights'+'\n')
		for i, m in enumerate(modules):
			# get the nodes names 
			subnet_nodes = vertex_index.loc[vertex_index['index'].isin(m.nodes())]['node'].tolist()
			# get the nodes propagation weights
			subnet_nodes_prop_weights = random_walk_weights_pvals.loc[random_walk_weights_pvals['node_index'].isin(m.nodes())]['prop_weight'].tolist()
			# write all nodes to file , separated by tab
			of.write('\t'.join(subnet_nodes))
			# write the size of the network and the sum of the propagation weights
			of.write(','+str(len(subnet_nodes))+','+str(round(sum(subnet_nodes_prop_weights),3)))
			of.write('\n')

			# plot module
			module_file = modules_plot_dir + "/module" + str(i) + ".pdf"
			
			# get the seed nodes for this module
			module_seed_nodes = [n for n in m.nodes() if n in seed_nodes]

#			plot_module(m,module_file,module_seed_nodes,vertex_index)

	of.close() 