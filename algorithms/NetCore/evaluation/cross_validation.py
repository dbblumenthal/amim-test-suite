### Gal Barel
### 6.6.2019
### compare normalization methods using cancer consensus genes

# impot 
import sys , os
from argparse import ArgumentParser
from random_walk_with_restart import *
from permutations_test import *
from auxiliary import *


# Parse arguments.
def get_parser():
	description = 'Run k-fold cross validation on network propagation starting from seed nodes (classification only)'
	parser = ArgumentParser(description=description)
	parser.add_argument('-e', '--edge_list', type=str, required=True, help='Edge list of the input network')
	parser.add_argument('-s', '--seed_list', type=str, required=True, help="List of seed nodes for evaluating the cross validations. If no weights are given, then all of these nodes will get a weights of 1.")
	parser.add_argument('-w', '--weights_file', type=str, required=False, default=None, help='Vertex weight file')
	parser.add_argument('-p', '--restart_prob', type=float, required=False, default=0.5, help='Random walk restart probability')
	parser.add_argument('-n', '--normalisaton_method', type=str,required=False, default='deg_norm',help='Which method to use for the normalisation of the random walk matrix. Options: deg_norm - degree only (same as Hotnet2), core_norm - normalize by the kshell, deg_core_diff - normalize by the difference of the degree and the kshell, deg_core_ratio - normalize by the ratio of the degree and the kshell')
	parser.add_argument('-fn', '--fold_num', type=int, required=False, default=0,help='Cross validation fold number (when running only one)')
	parser.add_argument('-nc', '--num_cross_valid_sets', type=int, required=False, default=5,help='Number of cross validation sets to be used')
	parser.add_argument('-np', '--num_permutations', type=int, required=False, default=100,help='Number of permutations')
	parser.add_argument('-nn', '--net_name', type=str, required=True, help='network name (for getting the permutations files')
	parser.add_argument('-pd', '--permut_dir', type=str, required=True, help='permutations directory that includes the edge lists for all network permutations')
	parser.add_argument('-pn', '--perm_name', type=str, required=False, default="edgelist", help='name of file of the edge list of the permuted networks')
	parser.add_argument('-od', '--output_dir', type=str, required=True, help='Output directory')
	return parser

def run_cross_valid_set(vertex_index,seed_nodes,weights_sorted,random_walk_mat,norm_adj_mat,nodes_order,node_deg_order,node_core_order,normalisaton_method,is_edge_weights,restart_prob,permut_dir,perm_name,net_name,num_perm):
	"""
		run random walk with restart for one cross validation set
	"""

	# create the random walk weights
	random_walk_weights = random_walk_mat.dot(weights_sorted)
	# run permutation test
	set_permutations = get_permutation_weights(permut_dir,perm_name,net_name,num_perm,norm_adj_mat,vertex_index,nodes_order,node_deg_order,node_core_order,weights_sorted,normalisaton_method,is_edge_weights,restart_prob)
	# get the pvalues
	set_p_vals = get_permutation_pvals(set_permutations,random_walk_weights.tolist(),nodes_order)

	return set_p_vals , set_permutations


# Run script.
def run(args):

	# create network and get properties
	vertex_index, G, adj_mat_dense , nodes_order , node_deg_order , node_core_order, removed_nodes, is_edge_weights = make_network(args.edge_list)

	# get the weights 
	if args.weights_file is None:
		# create binary weights based on seed list
		seed_nodes = get_seed_nodes(args.seed_list,vertex_index)
		# create binary weights - seed nodes will get a weight of 1, others 0
		binary_weights = [1 if node in seed_nodes else 0 for node in vertex_index['index']]
		weights = pd.DataFrame({'node':vertex_index['node'],'weight':binary_weights})
	else:
		# weights are given
		weights = get_weights(args.weights_file)

	# get the weights for all the nodes in the network
	weights_sorted , weights_mat = process_weights(weights,vertex_index,G,nodes_order)

	# get seed nodes
	if args.seed_list is None:
		# create seed nodes based on input weights
		seed_nodes = create_seeds_from_weights(weights_sorted,nodes_order)
	else:
		if args.weights_file is not None:
			seed_nodes = get_seed_nodes(args.seed_list,vertex_index)


	# get the random walk matrix 
	norm_adj_mat, random_walk_mat, subnets_file= get_random_walk_matrix(args.normalisaton_method,is_edge_weights,adj_mat_dense,node_deg_order,node_core_order,args.restart_prob,args.output_dir)
	
	print("Running cross validation fold " + str(args.fold_num))

	# create a cross validation set
	size = round(len(seed_nodes)/5)
	# randomly choose the training and testing sets
	test_nodes = np.random.choice(seed_nodes,size=size,replace=False)
	train_nodes = [node for node in seed_nodes if node not in test_nodes]
	group_vec = ["test" if n in test_nodes else "train" if n in train_nodes else "none" for n in nodes_order]

	# run random walk on the set and get the pvalues
	set_pvals , set_permutations = run_cross_valid_set(vertex_index,train_nodes,weights_sorted,random_walk_mat,norm_adj_mat,nodes_order,node_deg_order,node_core_order,args.normalisaton_method,is_edge_weights,args.restart_prob,args.permut_dir,args.perm_name,args.net_name,args.num_permutations)
	
	# save the pvalues for this set
	set_pvals_df = pd.DataFrame({"gene_index":nodes_order, "group":group_vec, "pvalue":set_pvals})
	set_pvals_df.to_csv(args.output_dir+args.normalisaton_method+"_fold"+str(args.fold_num)+".txt",sep="\t",header=True,index=False)

	print("Finished cross validation fold " + str(args.fold_num))

if __name__ == "__main__":
	run(get_parser().parse_args(sys.argv[1:]))