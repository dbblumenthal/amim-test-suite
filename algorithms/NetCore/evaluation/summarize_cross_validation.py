### Gal Barel
### 2.7.2019
### summary of the cross validation results

# impot 
import sys, os
from argparse import ArgumentParser
import numpy as np
import pandas as pd
from sklearn import metrics
from matplotlib import pyplot as plt
import seaborn as sns
sns.set()
sns.set(style="whitegrid")


# Parse arguments.
def get_parser():
	description = 'Summarize the results of the cross validation over all gene-disease sets'
	parser = ArgumentParser(description=description)
	parser.add_argument('-dp', '--disease_path', type=str, required=True, help='Path for the directory with the results of the diseases')
	return parser


def get_auc(disease_path,disese_name):
	"""
		get the AUC values for this disease
	"""
	
	# get file with cross valid values
	try:
		eval_values = pd.read_csv(disease_path+"method_compare.txt",sep="\t")
		methods = eval_values.method.unique().tolist()
	except FileNotFoundError:
		print("Disease " + disese_name + " has no method_compare file:")
		print("Disease will be skipped!")
		return pd.DataFrame()
	
	# use only the mean values of the corss validations
	eval_values_mean = eval_values.loc[eval_values['fold']=="mean"]
	
	# create a dataframe to save the auc values
	methods_auc = pd.DataFrame()
	
	for m in methods:
		m_values_mean = eval_values_mean.loc[eval_values_mean['method']==m]
		# AUC values for PR and ROC curves
		auc_roc = metrics.auc(m_values_mean['FPR'],m_values_mean['recall'])
		auc_pr = metrics.auc(m_values_mean['recall'],m_values_mean['precision'])
		m_auc = pd.DataFrame({"method":m,"roc":auc_roc,"pr":auc_pr},index=[0])
		methods_auc = methods_auc.append(m_auc)
	
	# add disease name 
	methods_auc['disease'] = disese_name
	
	return methods_auc


def plot_auc(all_auc,curve_type,disease_path):
	"""
		plot the AUC values for all diseaes for all methods
	"""

	if curve_type=="roc":
		plot_title="AUROC"
	else:
		plot_title="AUPRC"

	# create plot
	ax = sns.stripplot(x="method", y=curve_type, data=all_auc,jitter=True,size=5,palette="Set2")
	ax = sns.boxplot(x="method", y=curve_type, data=all_auc,whis=np.inf,width=0.5,palette="Set2")
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_ticks_position('bottom')
	plt.title(plot_title)
	ax.set(xlabel='', ylabel='')
	ax.set_ylim(1e-06)

	# change the color of the boxplot whiskers
	for i,artist in enumerate(ax.artists):
		# Set the linecolor on the artist to the facecolor, and set the facecolor to None
		col = artist.get_facecolor()
		artist.set_edgecolor(col)
		artist.set_facecolor('None')

		# Each box has 6 associated Line2D objects (to make the whiskers, fliers, etc.)
		# Loop over them here, and use the same colour as above
		for j in range(i*6,i*6+6):
			line = ax.lines[j]
			line.set_color(col)
			line.set_mfc(col)
			line.set_mec(col)

	plt.savefig(disease_path+"methods_"+curve_type+".pdf")
	plt.close()


def plot_mean_rank(ranking,disease_path):
	"""
		plot the mean of the ranking for each method in both categories
	"""

	# get the mean of the ranking
	roc_rank_summary = ranking[["method","roc_rank"]].groupby("method").agg(['mean', 'sum'])
	pr_rank_summary = ranking[["method","pr_rank"]].groupby("method").agg(['mean', 'sum'])

	# join roc and pr 
	rank_summary = roc_rank_summary.join(pr_rank_summary.reindex(roc_rank_summary.index, level=0))
	rank_summary.columns = [' '.join(col).strip() for col in rank_summary.columns.values]
	rank_summary['method'] = rank_summary.index

	# arrange df
	rank_summary_flat = rank_summary.melt(id_vars="method")
	rank_summary_flat['type'] = rank_summary_flat.variable.str.split(' ', expand = True)[0]
	rank_summary_flat['stat'] = rank_summary_flat.variable.str.split(' ', expand = True)[1]

	# plot the mean of the ranks for both 
	ax = sns.barplot(x="method",y="value",data=rank_summary_flat.loc[rank_summary_flat['stat']=='mean'],
					hue ="variable")
	plt.savefig(disease_path+"methods_mean_rank.pdf")
	plt.close()


def plot_rank_freq(ranking,disease_path):
	"""
		get the frequncey of the ranking for each method and plot
	"""

	# get the rank frequency for each method 
	method_rank_count = pd.DataFrame()
	for m in ranking.method.unique():
		m_rank = ranking.loc[ranking['method']==m]
		roc_rank_count = m_rank['roc_rank'].value_counts()
		pr_rank_count = m_rank['pr_rank'].value_counts()
		all_ranks = [1,2,3,4]

		if roc_rank_count.shape[0]!=4:
			# not all ranks are avilable for this method
			missing_ranks = set(all_ranks)-set(roc_rank_count.index)
			# add the missing ranks
			for i in missing_ranks:
				roc_rank_count.loc[i] = 0

		if pr_rank_count.shape[0]!=4:
			# not all ranks are avilable for this method
			missing_ranks = set(all_ranks)-set(pr_rank_count.index)
			# add the missing ranks
			for i in missing_ranks:
				pr_rank_count.loc[i] = 0

		m_rank_count = pd.DataFrame({"method":m, 'rank':all_ranks,"roc_count":roc_rank_count.sort_index(),'pr_count':pr_rank_count.sort_index()})
		method_rank_count = method_rank_count.append(m_rank_count)

	# plot ROC
	ax = sns.barplot(x="method",y="roc_count",data=method_rank_count,hue ="rank",palette="Blues_d")
	plt.title("AUROC")
	ax.set(xlabel='', ylabel='Frequency')
	plt.savefig(disease_path+"methods_rank_roc.pdf")
	plt.close()

	# plot ROC
	ax = sns.barplot(x="method",y="pr_count",data=method_rank_count,hue ="rank",palette="Greens_d")
	plt.title("AUPRC")
	ax.set(xlabel='', ylabel='Frequency')
	plt.savefig(disease_path+"methods_rank_pr.pdf")
	plt.close()


def run(args):

	diseases = [dI for dI in os.listdir(args.disease_path) if os.path.isdir(os.path.join(args.disease_path,dI))]

	print("Getting cross validation values")

	# get auc values from all diseases
	all_auc = pd.DataFrame()
	for d in diseases:
		auc_disease = get_auc(args.disease_path+d+"/",d)
		all_auc = all_auc.append(auc_disease)

	# plot AUROC and AUPRC
	plot_auc(all_auc,"roc",args.disease_path)
	plot_auc(all_auc,"pr",args.disease_path)

	# save the AUC values into dataframe
	all_auc.to_csv(args.disease_path+"methods_auc.txt",sep='\t',header=True, index=False)

	print("Getting method ranking")

	# get a ranking for each method in each disease
	ranking = pd.DataFrame()
	for d in diseases:
		# get the randking for this disease
		d_auc = all_auc.loc[all_auc['disease']==d]
		# rank tie breaker is min -> both values will get the lower option
		d_auc['roc_rank'] = d_auc['roc'].rank(ascending=False,method="min")
		d_auc['pr_rank'] = d_auc['pr'].rank(ascending=False,method="min")
		ranking = ranking.append(d_auc)

	plot_mean_rank(ranking,args.disease_path)
	plot_rank_freq(ranking,args.disease_path)

if __name__ == "__main__":
	sns.set(style="whitegrid")
	run(get_parser().parse_args(sys.argv[1:]))