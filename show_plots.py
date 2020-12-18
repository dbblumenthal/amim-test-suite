import pandas as pd
import numpy as np
from os import listdir
import seaborn as sns
import matplotlib.pyplot as plt
import mygene
import gseapy
import time
flatten = lambda l: [item for sublist in l for item in sublist]

mypath = "results/"
all_results = listdir(mypath)
def get_pathways(condition):
    if condition == "ALS":
        return ['hsa05014']
    elif condition == "LC":
        return ['hsa05223']
    elif condition== "UC":
        return ['hsa04060', 'hsa04630', 'hsa05321']
    elif condition== "HD":
        return ['hsa05016']
    elif condition== "CD":
        return ['hsa04621', 'hsa04060', 'hsa04630', 'hsa05321', 'hsa04140']
new_files = [x for x in all_results if "CUSTOM" in x]

for i in range(len(new_files)):
    if new_files[i].endswith(".csv"):
    if i == 0:
        results = pd.read_csv(mypath + new_files[i])
    else:
        temp = pd.read_csv(mypath + new_files[i])
        temp = temp.replace({"REWIRED": "EXPECTED_DEGREE"})
        temp = temp.replace({"RDPN": "REWIRED"})
        results = pd.concat([results,temp], axis=0, ignore_index=True)
results = results.replace({"GSE112680":"ALS","GSE30219": "LC","GSE75214": "UC", "GSE75214_cd": "CD", "GSE3790": "HD"})


all_genes = list(results.result_genes)
genes = []
for i in range(len(all_genes)):
    try:
        genes.append(all_genes[i].split(","))
    except:
        pass


genes = list(set(flatten(genes)))
mg = mygene.MyGeneInfo()
out = mg.querymany(genes, scopes="entrezgene", fields='symbol', species='human', verbose=False)
mapping = dict()
for line in out:
    try:
        mapping[line["query"]] = line["symbol"]
    except KeyError:
        print("{0} was not mapped to any gene name".format(line["query"]))
        mapping[line["query"]] = line["query"]

kegg_pval = []

for i in results.index :
    if results.neg_log_gsea_p_value[i] == -1:
        gs = str(results.result_genes[i])
        gs = gs.split(",")
        gs = [mapping[x] for x in gs]

        res = False
        while res == False:
            try:
                enr = gseapy.enrichr(gene_list=gs, description='pathway', gene_sets='KEGG_2016', outdir="out")
                res = True
            except:
                print("sleeping")
                time.sleep(1)
        full_results = enr.results
        terms = list(full_results.Term)
        terms = [x.split(' ')[-1] for x in terms]
        p_values = []
        pathways = get_pathways(results.condition_name[i])
        for j in range(len(terms)):
            if terms[j] in pathways:
                p_values.append(-np.log10(full_results['Adjusted P-value'][j]))

        kegg_pval.append(np.mean(p_values))

        print(i)

next = 0
for i in results.index:


    if results.neg_log_gsea_p_value[i] == -1:
        num = kegg_pval[next]
        if np.isnan(num):
            num = 0
        results.neg_log_gsea_p_value[i] = 0


disgenet = pd.read_csv("data/networks/disgenet.csv")
dis_ids = {"ALS":"C0002736", "LC": "C1737250", "UC": "C0009324","HD":"C0020179","CD":"C0021390"}
disgenet = disgenet[disgenet["disease_id"].isin(list(dis_ids.values()))]
overlaps = []
for i in results.index:
    try:
        gs = results.result_genes[i]
        condition = results["condition_name"][i]
        dis = disgenet[disgenet.disease_id == dis_ids[condition]]
        dis_genes = list(dis.gene)
        dis_genes = [str(x) for x in dis_genes]
        gs = gs.split(",")
        oc = len(set(gs).intersection(set(dis_genes)))/min(len(set(gs)), len(set(dis_genes)))
        overlaps.append(oc)
    except:
        overlaps.append(0)

results["disgenet_overlap"] = overlaps


fig, axes = plt.subplots(3, 1)
fig.subplots_adjust(hspace=0.8, wspace=0.4)
kv = {"x":"network_generator_name", "data": results,"hue": "condition_name", "palette": "Accent"}

sns.boxplot(ax = axes[0], y = "neg_log_gsea_p_value",**kv)
axes[0].set_title(r"$\bf{" + "a)" + "}$" + ' KEGG gene set enrichment', loc = "left")
axes[0].set(xlabel="")

sns.boxplot(ax = axes[1], y = "mean_mutual_information",**kv)
axes[1].set_title(r"$\bf{" + "b)" + "}$" +' Mean mutual information with the phenotype', loc = "left")
axes[1].legend_.remove()


sns.boxplot(ax = axes[2], y = "disgenet_overlap",**kv)
axes[2].set_title(r"$\bf{" + "c)" + "}$" +' Overlap with DisGeNET disease genes', loc = "left")
axes[2].legend_.remove()
axes[2].set(xlabel="")



plt.tight_layout()
# plt.savefig("../img/basic/all.png")
# plt.savefig("../img/basic/all.pdf")


plt.show()

