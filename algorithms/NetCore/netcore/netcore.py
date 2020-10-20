### Gal Barel
### 4.4.2019
### This file contains the main function fo running the algorithm

# import packages
import sys
from argparse import ArgumentParser
import pickle
import warnings

# import scripts
from auxiliary import *
from random_walk_with_restart import *
from create_modules import *
from permutations_test import *

# Parse arguments.
def get_parser():
	description = 'Network propagation with degree/k-shell correction'
	parser = ArgumentParser(description=description)
	parser.add_argument('-e', '--edge_list', type=str, required=True, help='Edge list of the input network')
	parser.add_argument('-w', '--weights_file', type=str, required=False, default=None, help='Vertex weight file')
	parser.add_argument('-s', '--seed_file', type=str, required=False, default=None, help="List of seed nodes for finding network modules")
	parser.add_argument('-p', '--restart_prob', type=float, required=False, default=0.5, help='Random walk restart probability')
	parser.add_argument('-n', '--normalisaton_method', type=str,required=False, default='core_norm',help='Which method to use for the normalisation of the random walk matrix. Options: deg_norm - degree only (same as Hotnet2), core_norm - normalize by the kshell, deg_core_diff - normalize by deg - core difference, deg_core_ratio - normalize by the deg-core ratio')
	parser.add_argument('-pd', '--permut_dir', type=str, required=True, help='permutations directory that includes the edge lists for all network permutations')
	parser.add_argument('-pn', '--perm_name', type=str, required=False, default="edgelist", help='name of file of the edge list of the permuted networks')
	parser.add_argument('-np', '--num_permutations', type=int, required=False, default=100,help='Number of permutations')
	parser.add_argument('-nn', '--net_name', type=str, required=False, default=None, help='network name (for getting the permutations files')
	parser.add_argument('-pt', '--pval_thresh', type=float, required=False, default=0.01, help='Pvalue threshold for adding nodes from the propagation results to the seed modules')
	parser.add_argument('-wt', '--weight_thresh', type=float, required=False, default=None, help='Propagaion weight threshold for adding nodes from the propagation results to the seed modules')
	parser.add_argument('-ms', '--max_subnet_size', type=int, required=False, default=None, help='Maximum number of nodes in the network that is induced by the seed nodes and the added nodes from the propagation')
	parser.add_argument('-o', '--output_dir', type=str, required=True, help='Output directory')
	return parser


# Run script.
def run(args):

	# check input is correct
	print("Processing input...")

	# check input
	check_input(args.restart_prob,args.normalisaton_method,args.output_dir,args.weights_file,args.seed_file)

	# get network and node properties

	vertex_index, G, adj_mat_dense , nodes_order , node_deg_order , node_core_order, removed_nodes, is_edge_weights = make_network(args.edge_list)
	vertex_index = vertex_index.astype({'node': 'object'})
	# get the weights 
	if args.weights_file is None:
		# create binary weights based on seed list
		seed_nodes = get_seed_nodes(args.seed_file,vertex_index)
		# create binary weights - seed nodes will get a weight of 1, others 0
		binary_weights = [1 if node in seed_nodes else 0 for node in vertex_index['index']]
		weights = pd.DataFrame({'node':vertex_index['node'],'weight':binary_weights})
	else:
		# weights are given
		weights = get_weights(args.weights_file)

	# get the weights for all the nodes in the network
	weights_sorted , weights_mat = process_weights(weights,vertex_index,G,nodes_order)

	# get seed nodes
	if args.seed_file is None:
		# create seed nodes based on input weights
		seed_nodes = create_seeds_from_weights(weights_sorted,nodes_order)
	else:
		if args.weights_file is not None:
			seed_nodes = get_seed_nodes(args.seed_file,vertex_index)

	# get the random walk matrix 
	norm_adj_mat, random_walk_mat, subnets_file= get_random_walk_matrix(args.normalisaton_method,is_edge_weights,adj_mat_dense,node_deg_order,node_core_order,args.restart_prob,args.output_dir)

	# get the weights and weights matrix
	random_walk_weights = random_walk_mat.dot(weights_sorted)

	# assume net name is the same as edge list name (if not given)
	if args.net_name is None:
		net_name = os.path.basename(args.edge_list).split('.')[0]
	else:
		net_name = args.net_name

	# get the pvalues based on a permutation test
	permutations = get_permutation_weights(args.permut_dir,args.perm_name,net_name,args.num_permutations,norm_adj_mat,vertex_index,nodes_order,node_deg_order,node_core_order,weights_sorted,args.normalisaton_method,is_edge_weights,args.restart_prob)

	# get the pvalues
	p_vals = get_permutation_pvals(permutations,random_walk_weights.tolist(),nodes_order)

	# create a DF with the propagation weights and pvalues
	nodes_order_name = vertex_index.loc[vertex_index['index'].isin(nodes_order)]['node'].tolist()
	weights_df = pd.DataFrame({"node_index":nodes_order,"node":nodes_order_name,"prop_weight":random_walk_weights,"pvalue":p_vals})
	weights_df_sorted = weights_df.sort_values(['pvalue', 'prop_weight'], ascending=[True, False])

	# save the random walk weights and pvalues
	weights_df_sorted.to_csv(args.output_dir+"random_walk_weights.txt",sep="\t",header=True,index=False)

	get_seed_modules(G,seed_nodes,weights_df,args.pval_thresh,args.weight_thresh,args.max_subnet_size,vertex_index,subnets_file,args.output_dir)

	print("DONE")

if __name__ == "__main__":
	warnings.filterwarnings("ignore", category=UserWarning)
	run(get_parser().parse_args(sys.argv[1:]))