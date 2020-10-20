### Gal Barel
### 17.6.2019
### evaluate cross validation results 

# import 
from auxiliary import get_seed_nodes
from random_walk_with_restart import *

import sys
from argparse import ArgumentParser
import seaborn as sns
from sklearn import metrics
from matplotlib import pyplot


# Parse arguments.
def get_parser():
	description = 'Run k-fold cross validation on network propagation starting from seed nodes (classification only)'
	parser = ArgumentParser(description=description)
	parser.add_argument('-e', '--edge_list', type=str, required=True, help='Edge list of the input network')
	parser.add_argument('-s', '--seed_list', type=str, required=True, help="List of seed nodes for starting the propagation from")
	parser.add_argument('-m', '--methods', type=str, nargs='+', required=False, default=["deg_norm","core_norm","deg_core_diff","deg_core_ratio"], help='Types of normalisation methods to be evaluated and compared')
	parser.add_argument('-p', '--pval_thresh', type=list, required=False, default=np.linspace(0.01,1.0,100), help='Pvalues thresholds for the evaluation')
	parser.add_argument('-f', '--fold_num', type=int, required=True, help='Number of cross validation folds')
	parser.add_argument('-d', '--cross_valid_dir', type=str, required=True, help='permutations directory that includes the edge lists for all network permutations')
	return parser


def get_test_negatives(all_negatives,test_nodes,G):
	"""
		choose a set of radnom nodes from all the negatives that will match the degree of the test node
	"""	
	test_negatives = []
	for n in test_nodes:
		n_deg = G.degree(n)
		# get a range of +-100 from the degree
		deg_range = [n_deg-100,n_deg+100]
		if deg_range[0]<0:
			deg_range[0]=0
		# get negative nodes with a degree in this range
		negatives_in_range = [n for n in all_negatives if G.degree(n)>deg_range[0] and G.degree(n)<deg_range[1]]
		if len(negatives_in_range)<2:
			# there are not enough negatives in the range
			# skip this node
			# FIXME : think of a better solution
			continue
		match_neg = np.random.choice(negatives_in_range,size=1)
		while(match_neg[0] in test_negatives):
			# don't add a negative that was already added
			# remove this one from the negatives in range 
			negatives_in_range.remove(match_neg[0])
			if len(negatives_in_range)==0:
				# cannot find more negatives that are not already in the list
				# skip this node
				break
			match_neg = np.random.choice(negatives_in_range,size=1)
		if len(negatives_in_range)==0:
			# cannot find more negatives that are not already in the list
			# skip this node
			continue
		test_negatives.append(match_neg[0])
	return test_negatives


def get_class_vals(pval,fold_pvals,test_nodes,test_negatives):
	"""
		get the classification evaluation values: FPR, recall and precision
		returns : 3 values for this specfic fold and pvalue
	"""
	 # get the nodes with a p-value under the threshold
	under_thresh = fold_pvals[fold_pvals['pvalue']<pval]['gene_index'].tolist()

	if len(under_thresh)<1:
		# there are no nodes with a p-value smaller than threshold
		# this pvalue cannot be evaluated 
		print("There are no nodes with a p-value below " + str(pval))
		return 0,0,0
				
	# count number of positives and negatives
	TP = len([n for n in test_nodes if n in under_thresh])
	FP = len([n for n in test_negatives if n in under_thresh])

	# get precision recall and TPR
	FPR = FP / len(test_negatives)
	recall = TP / len(test_nodes)
	precision = TP / len(under_thresh)

	return FPR, recall, precision


def eval_fold(fold_num,cross_valid_dir,method,pval_thresh,all_negatives,G):
	"""
		evaluate the results of this corss validation fold
		returns: 3 lists of the values for all different p-values thresholds
	"""

	# get fold pvalues
	fold_file=cross_valid_dir+method+"_fold"+str(fold_num)+".txt"
	fold_pvals = pd.read_csv(fold_file,sep="\t")

	# convert gene index from int to string
	fold_pvals.gene_index = fold_pvals.gene_index.apply(int)

	# get the training and test sets
	test_nodes = fold_pvals[fold_pvals['group']=="test"]['gene_index'].tolist()

	# get the negative set for this test set
	test_negatives = get_test_negatives(all_negatives,test_nodes,G)

	FPR = [0]*len(pval_thresh)
	recall = [0]*len(pval_thresh)
	precision = [0]*len(pval_thresh)

	for i,p in enumerate(pval_thresh):
		FPR[i], recall[i], precision[i] = get_class_vals(p,fold_pvals,test_nodes,test_negatives)

	return FPR, recall, precision


def eval_methods(methods,num_fold,cross_valid_dir,pval_thresh,all_negatives,G):
	"""
		evaluate all the cross validation folds for all methods
	"""

	# create df for all the results
	methods_fold_eval = pd.DataFrame()

	for method in methods:
		print("Evaluating " + method)
		# fold numbers start from 1 not 0
		for f in range(1,num_fold+1):
			FPR, recall, precision = eval_fold(f,cross_valid_dir,method,pval_thresh,all_negatives,G)
			f_df = pd.DataFrame({"method":[method]*len(pval_thresh), "fold":[f]*len(pval_thresh), "pval":pval_thresh,"FPR":FPR,"recall":recall,"precision":precision})
			methods_fold_eval = methods_fold_eval.append(f_df,ignore_index = True)

	methods_fold_eval = get_cross_valid_mean(methods,methods_fold_eval)

	# change the fold to be a categorical variable
	methods_fold_eval.fold = methods_fold_eval.fold.astype('category')

	return methods_fold_eval


def get_cross_valid_mean(methods,methods_fold_eval):
	"""
		get the means of all cross validation folds for each method
	"""

	for m in methods:
		method_pvals = methods_fold_eval[methods_fold_eval['method']==m]
		method_pvals_fold_mean = method_pvals.groupby(['pval']).mean()
		# add the mising columns
		method_pvals_fold_mean['pval'] = method_pvals_fold_mean.index
		method_pvals_fold_mean['fold']="mean"
		method_pvals_fold_mean['method']=m
		# add the means to the DF
		methods_fold_eval = methods_fold_eval.append(method_pvals_fold_mean,
													 ignore_index = True)

	return methods_fold_eval


def plot_method_fold(method,num_fold,methods_fold_eval,cross_valid_dir):
	"""
		plot the ROC and precision-recall curves of all cross validation folds (and mean) for one method
	"""

	# colors for plot
	# the mean will be black
	fold_colors = sns.cubehelix_palette(num_fold, start=.5, rot=-.75)
	fold_colors.append([0.7340253748558246, 0.16608996539792387, 0.20261437908496732])

	# precision recall curves
	pyplot.figure(figsize=(16, 12))
	ax = sns.lineplot(data=methods_fold_eval[methods_fold_eval['method']==method],
					  x="recall", y="precision",hue="fold",legend="full",
					 palette=fold_colors,style="fold",dashes=False,
					  markers=["o","o","o","o","o","D"])

	fig_file=cross_valid_dir+method+"_pre_rec.pdf"
	pyplot.savefig(fig_file)
	pyplot.close()

	# ROC curve
	pyplot.figure(figsize=(16, 12))
	ax = sns.lineplot(data=methods_fold_eval[methods_fold_eval['method']==method],
				  x="FPR", y="recall",hue="fold",legend="full",
				 palette=fold_colors,style="fold",dashes=False,
				  markers=["o","o","o","o","o","D"])

	fig_file=cross_valid_dir+method+"_roc.pdf"
	pyplot.savefig(fig_file)
	pyplot.close()


def plot_all_methods(methods,methods_fold_eval,cross_valid_dir):
	"""
		plot the ROC and precision-recall curves of all methods (mean only)
	"""

	# plot the means of all methods
	methods_fold_eval_means = methods_fold_eval[methods_fold_eval['fold']=="mean"]
	
	# get the area under the precision recall for each curve
	methods_auc = pd.DataFrame()
	for m in methods:
		m_mean = methods_fold_eval_means.loc[methods_fold_eval_means['method']==m,]
		m_auprc = metrics.auc(m_mean['recall'],m_mean['precision'])
		m_auroc = metrics.auc(m_mean['FPR'],m_mean['recall'])
		m_auc = pd.DataFrame({"method":m,"auprc":m_auprc,"auroc":m_auroc},index=[0])
		methods_auc = methods_auc.append(m_auc,ignore_index=True)

	# precision recall curve
	pyplot.figure(figsize=(16, 12))
	ax = sns.lineplot(data=methods_fold_eval_means,
					  x="recall", y="precision",hue="method",legend="full",
					 palette=sns.color_palette("Set2",len(methods)),style="method",dashes=False,
					 markers=True)
	pyplot.legend(title='Method - AUPRC', labels=methods_auc["method"] + " " + round(methods_auc['auprc'],3).map(str))

	fig_file=cross_valid_dir+"method_compare_pre_rec.pdf"
	pyplot.savefig(fig_file)
	pyplot.close()

	pyplot.figure(figsize=(16, 12))
	ax = sns.lineplot(data=methods_fold_eval_means,
				  x="FPR", y="recall",hue="method",legend="full",
				 palette=sns.color_palette("Set2",len(methods)),style="method",dashes=False,
				 markers=True)
	pyplot.legend(title='Method - AUROC', labels=methods_auc["method"] + " " + round(methods_auc['auroc'],3).map(str))

	fig_file=cross_valid_dir+"method_compare_roc.pdf"
	pyplot.savefig(fig_file)
	pyplot.close()



def run(args):
	# create network and get properties
	vertex_index, G, adj_mat_dense , nodes_order , node_deg_order , node_core_order, removed_nodes, is_edge_weights = make_network(args.edge_list)

	# get the seed nodes
	seed_nodes = get_seed_nodes(args.seed_list,vertex_index)

	# create a negative set - all nodes that are not seed nodes
	all_negatives = [n for n in nodes_order if n not in seed_nodes]

	# TODO : generate p-values thresholds based on the min and max of the pvalues 

	# evaluate the methods 
	methods_fold_eval = eval_methods(args.methods,args.fold_num,args.cross_valid_dir,args.pval_thresh,all_negatives,G)

	# save all values to file
	methods_fold_eval.to_csv(args.cross_valid_dir+"method_compare.txt",sep="\t",header=True,index=False)

	# plot for each method
	for m in args.methods:
		plot_method_fold(m,args.fold_num,methods_fold_eval,args.cross_valid_dir)

	# plot means of all methods together
	plot_all_methods(args.methods,methods_fold_eval,args.cross_valid_dir)

if __name__ == "__main__":
	sns.set(style="whitegrid")
	run(get_parser().parse_args(sys.argv[1:]))