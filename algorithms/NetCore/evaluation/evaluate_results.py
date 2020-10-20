### Gal Barel	
### 15.04.2019
### evaluate network propagation results

# import 
import sys
import os
from argparse import ArgumentParser
import pandas as pd
import networkx as nx
import numpy as np
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib import pyplot as plt
from sklearn import metrics 

# Parse arguments
def get_parser():
	description = 'Comapre the results of random walk with degree correction to the results of hotnet2 and hierarchical hotnet'
	parser = ArgumentParser(description=description)
	parser.add_argument('-co', '--consensus_file', type=str, required=True, help='File with gene symbols of the genes that are in the consensus list')
	parser.add_argument('-cn', '--candidate_file', type=str, required=True, help='File with gene symbols of the genes that are in the candidate list.')
	parser.add_argument('-e', '--edge_list', type=str, required=True, help='Edge list of the input network')
	parser.add_argument('-v', '--vertex_list', type=str, required=True, help='Vertex index file')
	parser.add_argument('-w', '--weights_file', type=str, required=True, help='Vertex weight file')
	parser.add_argument('-r', '--results_dir', type=str, required=True, help='The directory of the results files')
	parser.add_argument('-rn', '--results_names', type=str, nargs='+', required=True, help='List of the results names')
	parser.add_argument('-s', '--min_size', type=int, required=False, default=1, help='Minimum subnetwork size for evaluation')
	parser.add_argument('-o', '--output_dir', type=str, required=False, default=None, help='Output directory. If not mentioned is the same as results_dir')
	return parser


def create_consensus_list(consensus_file):
	# create consensus lists for comparison
	# file needs to contain the gene symbols in the first column
	consensus_gene = pd.read_csv(consensus_file,sep='\t',header=None)
	consensus_gene_list = consensus_gene[0].tolist()
	return consensus_gene_list


def create_candidate_list(candidate_file,G,vertex_index):
	# create a candidate list 
	# file needs to contain the gene symbols in the first column
	candidate_gene = pd.read_csv(candidate_file,sep='\t',header=None)
	candidate_gene_list = candidate_gene[0].tolist()
	return candidate_gene_list


def get_results(res_file,min_size):
	# get the results of the method , create a list of all the nodes that were identified
	
	# open the subnetworks file
	all_lines=[]
	for line in open(res_file):
		l=line.strip()
		if not l.startswith("#"):
			# skip commented lines
			all_lines.append(l)
	lines = [x.strip().split('\t') for x in all_lines]

	# create a list of all nodes that are in subnetworks of size > min_size
	all_nodes = []
	for l in lines:
		if len(l)>min_size:
			for n in l:
				all_nodes.append(n)
	
	# return only the unique nodes in the subnetworks 
	return list(set(all_nodes))


def eval_method(network_nodes,method_nodes,gs_list):
	# compare the nodes of the method to the nodes in the gold standard (gs) list

    # create a binary vector for the nodes, implying if they are in the gs list or not
    nodes_order = sorted(network_nodes)
    nodes_binary = [False]*len(nodes_order)
    for n in range(len(nodes_order)):
        if nodes_order[n] in gs_list:
            nodes_binary[n]=True
	
	# create a binary vector for the method nodes
    node_score = [0.0]*len(nodes_order)
    for n in range(len(nodes_order)):
        if nodes_order[n] in method_nodes:
            node_score[n] = 1.0
    
	# get the precision-recall matrix
	# the confusion matrix returns 4 values:
	# 		        Actual
	# 		  |  0  |  1  |
	# 	    ___________________
	# pred  0 |  TN | FP |
	# 	    1 |  FN | TP |

    predcit_group_mat = metrics.confusion_matrix(np.array(nodes_binary), np.array(node_score))

    tn, fp, fn, tp =  predcit_group_mat.ravel()

    if (tp==0 and fp==0):
    	# all of the predictions are negatives, there are no positive predictions
    	# print("All predictions are negatives!")
    	return 0,0

	
    precision = float(tp) / float(tp+fp)
    recall = float(tp) / float(tp+fn)
	
    return precision , recall


def get_all_results_eval(results_dir,results_names,min_size,network_nodes,consensus_list,candidate_list):

	# get the results from the results directory, according to the results names
	# the results directory should contain subdirectories that are have the given results names

	# create a dataframe to save all the results
	res_df = pd.DataFrame(columns=['method','consensus_pre','candidate_pre','consensus_rec','candidate_rec'])

	for n in results_names:
		res_dir = results_dir + "/" + n + "/"

		# get all of the results in for this subdirectory
		res_dir_all = os.listdir(res_dir)
		res_file = n +"_subnetworks.txt"

		for r in res_dir_all:
			r_full_path_file = res_dir + r + '/' + res_file
			r_nodes = get_results(r_full_path_file,min_size)
			
			# evaluate results - with consensus
			r_consensus_pre , r_consensus_rec = eval_method(network_nodes,r_nodes,consensus_list)

			# evaluate results - with candidate
			r_candidate_pre , r_candidate_rec = eval_method(network_nodes,r_nodes,candidate_list)	

			# add the results to the dataframe
			res_df = res_df.append({'method':n,'cond':float(r) ,'consensus_pre':r_consensus_pre,'candidate_pre':r_candidate_pre,'consensus_rec':r_consensus_rec,'candidate_rec':r_candidate_rec}, ignore_index=True)

	return res_df

def plot_precision(df,num_methods,plot_file):
    # plot the precision of all methods under 2 conditions
    # input : a df of the x and y vectors, grouped by the methods

    # figure settings
	fig, ax = plt.subplots()
	markers=['o', 'd','x','+','s','^','h','*','.','2'] 
	marker_count = 0
	ax.margins(0.05) #  adds 5% padding to the autoscaling
	df_groups = df.groupby('method')

	for name, group in df_groups:

		# color according to condition
		conds = group['cond'].tolist()
		cmap = cm.Blues(np.linspace(0,1,20))
		cmap = colors.ListedColormap(cmap[10:,:-1])
		norm = colors.Normalize(vmin=min(conds), vmax=max(conds))

		ax.scatter(group.iloc[:,2], group.iloc[:,3], marker=markers[marker_count], label=name, c = conds, cmap=cmap, norm=norm)
		marker_count += 1

	# add legend
	ax.legend(loc='upper right',numpoints = 1,fontsize=8)

	# set the color in the legend
	leg = ax.get_legend()
	for l in leg.legendHandles:
		l.set_color('navy')

	# plot
	plt.xlabel(df.columns[2], fontsize=10)
	plt.ylabel(df.columns[3], fontsize=10)
	fig.savefig(plot_file)
	

# Run script.
def run(args):

	# set output dir if needed
	if args.output_dir==None:
		args.output_dir=args.results_dir
	
	# create network
	print("Getting network information...")
	G=nx.read_edgelist(args.edge_list)
	
	# get the list that mapps the vertex index to their gene name
	vertex_index = pd.read_csv(args.vertex_list,sep=r"\s+", header=None)
	vertex_index.columns = ['index','node'] 
	vertex_index['index'] = vertex_index['index'].apply(lambda x : str(int(x)))

	# get all the nodes from the network for the evaluation 
	vertex_index_in_graph = vertex_index.loc[vertex_index['index'].isin(G.nodes())]
	vertex_index_in_graph = vertex_index_in_graph.dropna()
	network_nodes = vertex_index_in_graph['node'].tolist()
	
	# get the consensus list
	print("Getting consensus and candidate gene lists...")
	consensus_list = create_consensus_list(args.consensus_file)
	
	# create a candidate list
	candidate_list = create_candidate_list(args.candidate_file,G,vertex_index)

	# remove genes that are in consensus list from the candidate list
	common_in_lists = list(set(consensus_list) & set(candidate_list))
	candidate_list = [n for n in candidate_list if n not in common_in_lists]

	# remove genes from both lists that are NOT in the network !!!
	consensus_list = [n for n in consensus_list if n in network_nodes]
	candidate_list = [n for n in candidate_list if n in network_nodes]

	print("removed genes from consensus and candidate lists that are not in the network")
	
	# get the nodes from all methods
	print("Evaluating results...")
	res_pre_rec_df = get_all_results_eval(args.results_dir,args.results_names,args.min_size,network_nodes,consensus_list,candidate_list)

	# plot pre-rec for consensus
	pre_rec_consensus = res_pre_rec_df.loc[:,['method','cond','consensus_rec','consensus_pre']]
	plot_precision(pre_rec_consensus,len(args.results_names),args.output_dir+"/precision_recall_consensus.pdf")

	# plot pre-rec for candidate
	pre_rec_candidate = res_pre_rec_df.loc[:,['method','cond','candidate_rec','candidate_pre']]
	plot_precision(pre_rec_candidate,len(args.results_names),args.output_dir+"/precision_recall_candidate.pdf")

	# plot pre-pre of consensus vs candidate
	pre_consensus_candidate = res_pre_rec_df.loc[:,['method','cond','consensus_pre','candidate_pre']]
	plot_precision(pre_consensus_candidate,len(args.results_names),args.output_dir+"/precision_comparison.pdf")

	# plot rec-rec of consensus vs candidate
	rec_consensus_candidate = res_pre_rec_df.loc[:,['method','cond','consensus_rec','candidate_rec']]
	plot_precision(rec_consensus_candidate,len(args.results_names),args.output_dir+"/recall_comparison.pdf")

	print("DONE")
	
if __name__ == "__main__":
	run(get_parser().parse_args(sys.argv[1:]))