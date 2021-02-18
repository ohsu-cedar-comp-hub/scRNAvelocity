import matplotlib.pyplot as plt
import numpy as np
from sklearn.manifold import TSNE
import os
import datetime
import scvelo as scv
import scanpy
import sys
import pandas as pd
from scipy import stats
from statannot import add_stat_annotation
from functools import reduce
import seaborn as sns
import statsmodels.api as sm

curr_dir = os.getcwd()
begin_time = datetime.datetime.now().timestamp()

sys.stderr.write("beginning analyze!")

in_object         = snakemake.input.in_object
contrast		  = snakemake.params.contrast_name
out_object        = snakemake.output.out_object
out_dir           = os.path.dirname(out_object)

genes_of_interest = snakemake.params.genes
condition         = snakemake.params.seurat_status
cluster           = snakemake.params.seurat_cluster
Color_hex         = snakemake.params.color_hex
Order_plot        = snakemake.params.order_plot

color_dict        = dict(zip(Order_plot,Color_hex)) 


###################################################################

adata             = scv.read(in_object)
os.chdir(out_dir)
cond_levels       = list(pd.unique(adata.obs[condition]))

top_n_genes = 30

stats_order = []
for i in range(len(cond_levels)):
	for j in range(i+1,len(cond_levels)):
		stats_order.append((cond_levels[i],cond_levels[j]))

if not os.path.exists("figures"):
	os.makedirs("figures")

#loop through the clusters
#analyze the condition
#adata.obs['cluster_Condition'] = adata.obs.apply(lambda x: "{}_{}".format(x['cluster'],x[condition]), axis=1)

scv.tl.rank_velocity_genes(adata, groupby="cluster", n_genes = adata.n_vars)

scanpy.tl.rank_genes_groups(adata, groupby = "cluster",layer = 'velocity', method = "wilcoxon", use_raw = False)
results = adata_subset.uns['rank_genes_groups']
groups = results['names'].dtype.names
keys = list(results.keys())
keys = keys[1:]
all_dicts = {group : pd.DataFrame({"{}".format(key) : results[key][group] for key in keys}) for group in groups}
for df in all_dicts.keys():
	all_dicts[df]["cluster"] = df
top_n_genes_cluster = {key : list(value.nlargest(top_n_genes,"scores").names) for key,value in all_dicts.items()}

##########################################


def get_velocity_table(adata):
	"""
	input: a scvelo or scanpy object
	
	Returns a merged pandas data frame of 'rank_velocity_genes'
	
	Must run 'scv.tl.rank_velocity_genes(adata)' before running
	"""
	df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
	df_scores = scv.DataFrame(adata.uns['rank_velocity_genes']['scores'])
	all_dicts = dict()
	for col in df.columns:
		all_dicts[col] = pd.DataFrame()
		all_dicts[col]["genes"] = df[col]
		all_dicts[col]["scores"] = df_scores[col]
	all_list = [value for key,value in all_dicts.items()] 
	df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['genes'],how='inner'), all_list)
	df_merged.columns = ["genes"] + list(df.columns)
	return(df_merged)


df = get_velocity_table(adata)
df.to_csv("top_cluster_velocity_genes.tsv",sep="\t")

unique_clusters = pd.unique(adata.obs["cluster"])
n_num = 4
genes_cluster = dict()
for clust in unique_clusters:
	df_sort = df.sort_values(by = [clust], ascending = False)
	genes_list = df_sort['genes'].head(n=n_num).to_list() + df_sort['genes'].tail(n=n_num).to_list()
	genes_cluster[clust] = genes_list

all_dicts = {}
for cont in pd.unique(adata.obs[condition]):
	adata_subset = adata[adata.obs[condition].isin([str(cont)])]
	scv.tl.rank_velocity_genes(adata_subset, groupby="cluster", n_genes = adata_subset.n_vars)
	df = get_velocity_table(adata_subset)
	all_list[cont] = df
	#df.to_csv("top_{}_cluster_velocity_genes.tsv".format(cont),sep="\t")
all_list = [value for key,value in all_list.items()]
all_df = reduce(lambda  left,right: pd.merge(left,right,on=['genes'],how='inner'), all_list)
all_df.columns = [x.replace("_x","_HD").replace("_y","_FPD") for x in list(all_df.columns)]
all_df = pd.melt(all_df, id_vars = ["genes"], value_vars = list(all_df.columns)[1:])
all_df["cell_type"] = [x.split("_")[0] for x in all_df["variable"]]
all_df["condition"] = [x.split("_")[1] for x in all_df["variable"]]
pplot = sns.violinplot(x = 'cell_type', y = "value", data = all_df, inner = "box", rotation = 90, hue = "condition",split=True, palette = color_dict)
pplot.set_xticklabels(pplot.get_xticklabels(),rotation = 90)
pplot.set(ylabel = "velocity score")
pplot.get_figure().savefig("figures/test_cell_type_violin_velocity_scores.png",bbox_inches='tight')

keys = 'velocity_length', 'velocity_confidence'


def return_root(adata, root_cluster = "HSC",method="topleft", obs_column = "cluster"):
	"""
	adata: scvelo or scanpy object
	root_cluster: [str] name of cluster in the 'obs_column' of the metadata to find a ce
	method: the method of returning a cell within the cluster
	
	Retuns a cell name
	"""
	
	potential_roots = [i for i,x in enumerate(adata.obs['cluster']) if x == root_cluster]
	potential_umap = np.vstack([x for i,x in enumerate(adata.obsm["umap"]) if i in potential_roots])
	ysd_value = np.std(potential_umap[:,1])
	xsd_value = np.std(potential_umap[:,0])
	max_std = np.argmax(np.array(xsd_value,ysd_value))
	ymax_value = np.max(potential_umap[:,1])
	xmin_value = np.min(potential_umap[:,0])
	xmean_value = np.mean(potential_umap[:,0])
	ymean_value = np.mean(potential_umap[:,1])
	xmed_value = np.median(potential_umap[:,0])
	ymed_value = np.median(potential_umap[:,1])
	
	def center_return(potential_umap,max_std,ymed_value,xmed_value,xmean_value,ymean_value,i,j):
		
		putative_root = potential_umap[np.where( ((xmed_value-(xsd_value*i)) <= potential_umap[:,0]) & (potential_umap[:,0] <= (xmed_value+(i*xsd_value))) & (potential_umap[:,1] <= (ymed_value*(j*ysd_value))) & (potential_umap[:,1] > (ymed_value-(j*ysd_value))) & (potential_umap[:,0] >= (xmean_value-(xsd_value*i))) & (potential_umap[:,1] <= (ymean_value+(ysd_value*j))) )]
		if putative_root.shape[0] < min(10,int((len(potential_umap)/100)+1)):
			return center_return(potential_umap,max_std,ymed_value,xmed_value,xmean_value,ymean_value,i+0.1,j+0.1)
		elif putative_root.shape[0] > int(potential_umap.shape[0]/2):
			print(i)
			return np.argwhere(adata.obsm["umap"] == putative_root[int(putative_root.shape[0]/2)])[0][0]
		else:
			print(i,j)
			return np.argwhere(adata.obsm["umap"] == putative_root[int(putative_root.shape[0]/2)])[0][0]
	def max_expression(potential_roots):
		rowsums_cells = adata.layers['matrix'][potential_roots,:].sum(axis = 1)
		rowsums_array = np.squeeze(np.asarray(rowsums_cells))
		rowmax = np.argmax(rowsums_array)
		return potential_roots[rowmax]
		
	if type(method) == int:
		if method > len(potential_roots):
			method = len(potential_roots)-1
		if method < 0:
			method = 0
		return potential_roots[method]
	elif method == "topleft":
		return top_left_return(potential_umap,max_std,ymax_value,xmin_value,xmean_value,ymean_value,1,1)
	elif method == "center":
		return center_return(potential_umap,max_std,ymed_value,xmed_value,xmean_value,ymean_value,1,1)
	elif method == "max_x":
		return max_expression(potential_roots)
	elif method == "any":
		return(potential_roots[0])
	elif method == "random":
		return potential_roots[np.random.randint(1,len(potential_roots))-1]
	else:
		return np.argwhere(adata.obsm["umap"] == potential_roots[0])[0][0]

key_root = return_root(adata, root_cluster = "HSC", method = "max_x")
key_name = adata.obs.index[key_root]
scv.tl.velocity_pseudotime(adata, root_key = key_root)
scv.pl.scatter(adata, color='velocity_pseudotime', color_map='gnuplot', save = "pseudotime.png")
scv.tl.score_genes_cell_cycle(adata)
adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']
scv.tl.paga(adata, groups='cluster',root_key = key_name)
df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
df.to_csv("PAGA_cluster_transitions.tsv",sep="\t")

scv.pl.paga(adata, basis='umap', size=50, alpha=.1,min_edge_width=2, node_size_scale=1.5, save = "PAGA.png")

for clust in unique_clusters:
	adata_subset = adata[adata.obs["cluster"].isin([str(clust)])]
	scv.tl.rank_velocity_genes(adata_subset, groupby=condition, n_genes = adata_subset.n_vars)
	df = scv.DataFrame(adata_subset.uns['rank_velocity_genes']['names'])
	df_scores = scv.DataFrame(adata_subset.uns['rank_velocity_genes']['scores'])
	
	all_dicts = dict()
	for col in df.columns:
		all_dicts[col] = pd.DataFrame()
		all_dicts[col]["genes"] = df[col]
		all_dicts[col]["scores"] = df_scores[col]
	#df_all = pd.concat([df,df_scores], axis = 1, join = 'inner', ignore_index = False, keys = ["genes","scores"], sort = False)	
	all_list = [value for key,value in all_dicts.items()] 
	df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['genes'],how='inner'), all_list)
	df_merged.columns = ["genes"] + list(df.columns)
	#df_merged = df_merged.sort_values(by = df.columns[0])
	df_merged.to_csv("top_{}_{}_velocity_genes.tsv".format(clust.replace("/","."),condition) ,sep="\t")
	df.to_csv("top_{}_velocity_genes.tsv".format(clust.replace("/",".")), sep="\t")
	more_genes = set(pd.unique(df.head(4).values.ravel('K'))) 
	less_genes = set(pd.unique(df.tail(4).values.ravel('K'))) 
	#more_genes = set(list(df_merged['genes']))
	genes_of_interest = list(set(genes_of_interest) | more_genes | less_genes | set(genes_cluster[clust]) )
	phase_colors = {'G1':'red', 'G2M':'blue', 'S':'green'}
	
	for gene in genes_of_interest:
		try:
			scv.pl.velocity(adata_subset,str(gene), dpi = 120, figsize = (7,5), color = "phase",legend_loc = 'best',save = "scatter_phase_{}_{}.png".format(str(clust),gene), palette = phase_colors)
			scv.pl.velocity(adata_subset,str(gene), dpi = 120, figsize = (7,5), color = condition,legend_loc = 'best',save = "scatter_gene_cluster_{}_{}.png".format(str(clust),gene), palette = color_dict)
		except:
			sys.stderr.write("{} not included in {}, ".format(gene,str(clust)))
	sys.stderr.write("\n")
	try:
		scv.pl.scatter(adata_subset, c=keys, cmap='coolwarm', perc=[5, 95], save = "scatter_confidence_{}.png".format(str(clust)))
		scv.pl.velocity_embedding_stream(adata_subset, basis='umap', color = condition, palette = color_dict,save = "scvelo_stream_{}.png".format(str(clust)))
		#scv.pl.velocity_embedding(adata_subset, basis='umap', color = 'Condition', palette = color_dict,arrow_length=0, arrow_size=0, dpi=120, save = "scvelo_embedding_{}.png".format(str(clust)))
		scanpy.pl.violin(adata_subset, keys = "velocity_length", groupby = condition, save = "velocity_length_{}.png".format(str(clust)), palette = color_dict)
		scanpy.pl.violin(adata_subset, keys = "velocity_length", groupby = "orig_ident", save = "orig_velocity_length_{}.png".format(str(clust)), rotation = 90)
	except:
		sys.stderr.write("Error in one of the plots\n")
	
	Count = adata_subset.obs.groupby(condition)['velocity_length'].count()
	Max = adata_subset.obs["velocity_length"].max()
	Means = adata_subset.obs.groupby(condition)["velocity_length"].mean()
	Median = adata_subset.obs.groupby(condition)["velocity_length"].median()
	p_plot = scanpy.pl.violin(adata_subset, keys = "velocity_length", groupby = condition, show = False, inner = "box", size = 0, order = Order_plot, palette = color_dict)
	
	fig = add_stat_annotation(p_plot, data = adata_subset.obs, x=condition, y = "velocity_length", box_pairs = [tuple(Order_plot)], test = 't-test_welch', text_format = "full", loc = 'outside')
	
	add_fig = fig[0]
	for i,x in enumerate(Count):
		add_fig.text(Order_plot.index(Count.index.values[i]),Max+1,"{}".format(x))
	add_fig.get_figure().savefig("figures/stats_velocity_length_{}.png".format(str(clust)))
	
	### additional testing specific to two conditions (ie. wild-type vs mutant)
	#adata_subset.layers['velocity'] is the velocity residual for each gene per cell. axis 0 = gene, axis 1 = cell
	valid_var = [np.isnan(adata_subset.layers['velocity'][:,x]).any() for x,i in enumerate(adata_subset.var_names)] #any velocities in which genes have NAN ?  False means there is no NAN and has a velocities for the cells in that feature
	valid_gene_names = [adata_subset.var_names[i] for i,x in enumerate(valid_var) if x == False]
	valid_index = [i for i,x in enumerate(valid_var) if x == False]
	valid_gene_index = dict(zip(valid_gene_names,valid_index)) 
	condish = list(adata_subset.obs[condition])
	df2 = {"sum_values" : [], "mean_values" : [], "condition" : []}
	valid_velocities = adata_subset.layers['velocity'][:,[i for i,x in enumerate(valid_var) if x == False]] # -> [array of cells : array of velocities of valid genes]
	top_genes = top_n_genes_cluster[clust]
	condition_dict = {}
	for j,y in enumerate(condish): #for y in in conditions
		#append the mean/sum of the velocities of all valid genes. valid genes being ones that have velocities / spliced/unspliced
		df2["sum_values"].append(sum(valid_velocities[j,:]))
		df2["mean_values"].append(np.mean(valid_velocities[j,:]))
		df2["condition"].append(y)
		if y not in condition_dict:
			condition_dict[y] = [j]
		else:
			condition_dict[y].append(j)
	cond_mean_genes = {}
	for key,value in condition_dict.items():
		cond_mean_genes[key] = np.mean(valid_velocities[value,:], axis = 0)  #mean of each genes velocity across the FPD cells of this cluster
	cond_mean_genes["gene"] = valid_gene_names
	
	df = pd.DataFrame(df2)
	pplot = sns.violinplot(x = 'condition', y = 'sum_values', data = df, inner = "box")
	fig2 = add_stat_annotation(pplot, data = df, x='condition', y = "sum_values", box_pairs = [tuple(Order_plot)], test = 't-test_welch', text_format = "full", loc = 'outside')
	fig2[0].get_figure().savefig("figures/violin_sum_{}_plot.png".format(str(clust)))
	fig2[0].get_figure().clf()
	pplot.get_figure().clf()
	
	pplot = sns.violinplot(x = 'condition', y = 'mean_values', data = df, inner = "box", palette = color_dict,order = Order_plot)
	fig2 = add_stat_annotation(pplot, data = df, x='condition', y = "sum_values", box_pairs = [tuple(Order_plot)], test = 't-test_welch', text_format = "full", loc = 'outside')
	fig2[0].get_figure().savefig("figures/violin_mean_{}_plot.png".format(str(clust)))
	fig2[0].get_figure().clf()
	pplot.get_figure().clf()
	
	df3 = pd.DataFrame(cond_mean_genes)

	try:
		stat,p = stats.f_oneway(valid_velocities[condition_dict[list(condition_dict.keys())[0]],:],valid_velocities[condition_dict[list(condition_dict.keys())[1]],:])
		df3['anova'] = stat
		df3['p_val'] = p
		df3['p_val_Adjusted'] = sm.stats.multipletests(p,method = "fdr_bh")[1]
		
	except:
		sys.stderr.write("Error in Anova\n")
	
	#wilcoxon rank-sum test on velocity values
	scanpy.tl.rank_genes_groups(adata_subset, groupby = condition,layer = 'velocity', method = "wilcoxon", use_raw = False)
	results = adata_subset.uns['rank_genes_groups']
	groups = results['names'].dtype.names
	keys = list(results.keys())
	keys = keys[1:]
	results = pd.DataFrame({"{}".format(key) : results[key][group] for group in groups for key in keys if group == contrast})
	results.to_csv("wilcoxon_{}_velocity.tsv".format(clust.replace("/",".")), sep="\t")
	
	#wilcoxon rank-sum test on total RNA
	scanpy.tl.rank_genes_groups(adata_subset, groupby = condition, method = "wilcoxon")
	RNA_results = adata_subset.uns['rank_genes_groups']
	RNA_results = pd.DataFrame({"{}".format(key) : RNA_results[key][group] for group in groups for key in keys if group == contrast})
	
	# merge velocity and total RNA tests
	RNA_velocity = pd.merge(results,RNA_results,on = ['names'])
	RNA_velocity.columns = [k.replace("x","Vel").replace("y","RNA") for k in list(RNA_velocity.columns)]
	RNA_velocity.to_csv("wilcoxon_{}_RNA_velocity_merged.tsv".format(clust.replace("/",".")), sep="\t")
	# all genes 
	df3.to_csv("anova_{}_mean_velocity_genes.tsv".format(clust.replace("/",".")), sep="\t")
	adata.uns["DE_velo_RNA_{}".format(clust.replace("/","."))] = RNA_velocity
	
adata.write_h5ad(in_object)
os.chdir(curr_dir)
adata.obs.to_csv(out_object,sep="\t")


# get comparision for genes between clusters
"""

Color_hex = ["#67CAD4","#B0207E"] #purple for FPD, teal for HD
Order_plot = ["HD","FPD"]

color_dict        = dict(zip(Order_plot,Color_hex))



heatmap_data = pd.pivot_table(all_df, values='value', 
                              index='genes', 
                              columns='cell_type')
#heat_plot = sns.clustermap(heatmap_data, annot=False,center=0,xticklabels=True, yticklabels=True,figsize=(12,20), col_cluster = False)
		#heat_plot.ax_heatmap.set_ylabel("")
#heat_plot.cax.set_visible(False)
#heat_plot.savefig("figures/test_heatmap_velocity_scores.png",bbox_inches='tight')


	
all_dfs = [value for key,value in all_dicts.items()] 
all_keys = [key for key,value in all_dicts.items()]
all_merged = reduce(lambda  left,right: pd.merge(left,right,on=['genes'],how='inner', suffixes=('_x', '_y')), all_dfs)
reduce(lambda  left,right: pd.merge(left,right,on=['genes'],how='inner', suffixes=('_x', '_y')), zip(all_dfs,all_keys))
 #reduce(lambda left,right : pd.merge(all_dicts[left],all_dicts[right],on=['names'],how='inner', suffixes=('_{}'.format(left), '_{}'.format(right))), list(all_dicts))
results = pd.DataFrame({"{}".format(key) : results[key][group] for group in groups for key in keys if group == contrast})
results_df = pd.DataFrame({"{}".format(key) : results[key][group] for group in groups for key in keys})
top_30_genes = {x : list(FPD_df.nlargest(30,x).genes) for x in FPD_df.columns[1:]}
top_30_genes2 = {x : list(HD_df.nlargest(30,x).genes) for x in HD_df.columns[1:]}
top_30_genes_cluster = {key : list(value.nlargest(30,"scores").names) for key,value in all_dicts.items()}


	

def compare_top_n_genes(adata,cluster1,cluster2,top_dict,condition = "Condition"):
	"""
	#compare two clusters velocity scores 
	#returns a dataframe of ranksum scores of top genes
	#example usage: df_test = compare_top_n_genes(adata,"HSC","GMP",top_30_genes_cluster)
	"""
	top1 = top_dict[cluster1]
	top2 = top_dict[cluster2]
	adata_subset1 = adata.obs["cluster"].isin([str(cluster1)])
	index_subset1 = [i for i,x in enumerate(adata_subset1) if x == True]
	adata_subset2 = adata.obs["cluster"].isin([str(cluster2)])
	index_subset2 = [i for i,x in enumerate(adata_subset2) if x == True]
	genes_in_both_clusters = [x for x in top1 if x in top2]
	all_genes = top1 + top2
	df = {}
	df["gene"] = []
	df["stat"] = []
	df["p"] = []
	for gene in genes_in_both_clusters:
		df["gene"].append(gene)
		stat,p = stats.ranksums(adata.layers['velocity'][index_subset1,adata.var_names.get_loc(gene)],adata.layers['velocity'][index_subset2,adata.var_names.get_loc(gene)])
		df["stat"].append(stat)
		df["p"].append(p)
	df["adj_p"] = sm.stats.multipletests(df["p"],method = "fdr_bh")[1]
	
	
	condish = list(adata.obs[condition])
	condition_dict = {}
	adata.obs[condition]
	for cond in set(condish):
		condition_dict[cond] = [i for i,x in enumerate(condish) if x == cond]
	list(set(condition_dict[cond]) & set(index_subset2))
		subset1 = adata.layers['velocity'][index_subset1,adata.var_names.get_loc(gene)]
		subset2 = adata.layers['velocity'][index_subset2,adata.var_names.get_loc(gene)]
		data = {"condition" : list(adata.obs[condition][index_subset1]) + list(adata.obs[condition][index_subset2]), 
		"cluster" : list(adata.obs["cluster"][index_subset1]) + list(adata.obs["cluster"][index_subset2]),
		"velocity" : list(subset1) + list(subset2)}
		maxsize = max([len(a) for a in data.values()])
		data_pad = {k:np.pad(v, pad_width=(0,maxsize-len(v),), mode='constant', constant_values=np.nan) for k,v in data.items()}
		df = pd.DataFrame(data_pad)
		
		data_f = pd.DataFrame(data)
		pplot = sns.violinplot(x = "cluster", y = "velocity", data = data,inner = "box", hue = "condition",split=True, palette = color_dict)
		pplot.set_xticklabels(pplot.get_xticklabels(),rotation = 90)
		pplot.set_title(gene)
		pplot.set(ylabel = "velocity score")
	
	return df



for top_genes
plot_gene_velocity(adata,top_genes["HSC"][0], condition = "Condition", hue_order = Order_plot, palette = color_dict)


pplot = sns.violinplot(x = "cluster", y = "velocity", data = df,inner = "box", hue = "condition",split=True, hue_order = Order_plot, palette = color_dict)	
	
	
data = {{"condition" : adata.obs[condition][indexes] for indexes in dict_of_clusters.values()}  , "cluster" : [adata.obs["cluster"][indexes] for indexes in dict_of_clusters.values() ], "velocity" : [adata.layers['velocity'][cluster_subset,adata.var_names.get_loc(gene)] for cluster_subset in dict_of_clusters.values() ]}

pplot = sns.violinplot(x = "cluster", y = "velocity", data = pd.DataFrame(data),inner = "box", hue = "condition",split=True, hue_order = Order_plot, palette = color_dict)

scv.pl.scatter(adata, color='velocity_pseudotime', color_map='gnuplot', save = "test_pseudotime.png")




			
df = pd.DataFrame(data)
for key,value in top_genes.items():
	current_key = df[df[
	pplot = sns.violinplot(x = "cluster", y = "velocity", data = df,inner = "box", hue = "condition",split=True, hue_order = Order_plot, palette = color_dict)
	pplot.set_xticklabels(pplot.get_xticklabels(),rotation = 90)
	pplot.set_title(gene)
	pplot.set(ylabel = "velocity score")
	pplot.get_figure().savefig("figures/violin_gene_{}_cluster_{}_velocity_scores.png".format(gene,key),bbox_inches='tight')
	pplot.get_figure().clf()

key = "HSC"
pplot.get_figure().savefig("figures/violin_cluster_{}_velocity_scores.png".format(key),bbox_inches='tight')
sns.violinplot(x = "cluster", y = "velocity", data = pd.DataFrame(data),inner = "box", hue = "condition",split=True, hue_order = Order_plot, palette = color_dict)
g = sns.FacetGrid(HSC, col = "gene", hue = "condition", palette = color_dict,hue_order = Order_plot)
pplot = g.map(sns.violinplot, "cluster", "velocity", split = True, inner = "box")
g.savefig("figures/violin_cluster_{}_velocity_scores.png".format(key), format="png")


for top_key,top_value in top_genes.items():
	dict_of_clusters = {each_cluster : [i for i,x in enumerate(adata.obs[cluster]) if x == each_cluster] for each_cluster in pd.unique(adata.obs[cluster])}
	for gene in top_value:
		data = {"condition" : [], "velocity" : [], "cluster" : []}
		for key,value in dict_of_clusters.items():
			n_len = len(list(adata.obs[groupby][value]))
			data["condition"] = data["condition"] + list(adata.obs[groupby][value])
			data["cluster"] = data["cluster"] + list(adata.obs[cluster][value])
			data["velocity"] = data["velocity"] + list(adata.layers['velocity'][value,adata.var_names.get_loc(gene)])

for top_key,top_value in top_genes.items():
	for gene in top_value:
		plot_gene_velocity(adata,gene,groupby = "Condition",cluster = "cluster", clust = top_key,hue_order = Order_plot, palette = color_dict)
		#pplot.get_figure().savefig("figures/violin_gene_{}_cluster_{}_velocity_scores.png".format(gene,top_key.replace("/",".")),bbox_inches='tight')
		#pplot.get_figure().clf()
		
for gene in GOI:
	plot_gene_velocity(adata,gene = gene, groupby = "Condition", cluster = "cluster",hue_order = Order_plot, palette = color_dict)
	
gene = "CD74"	
data = {"condition" : [], "velocity" : [], "cluster" : [], "sample" : []}
for key,value in dict_of_clusters.items():
	data["condition"] = data["condition"] + list(adata.obs[groupby][value])
	data["cluster"] = data["cluster"] + list(adata.obs[cluster][value])
	data["sample"] = data["sample"] + list(adata.obs['orig_ident'][value])
	data["velocity"] = data["velocity"] + list(adata.layers['velocity'][value,adata.var_names.get_loc(gene)])
	
HSC = df[df["cluster"] == "HSC"]
pplot = sns.violinplot(x = "sample", y = "velocity", data = df,inner = "box", hue = "condition",split=True, hue_order = Order_plot, palette = color_dict)
pplot.set_xticklabels(pplot.get_xticklabels(),rotation = 90)
pplot.set_title(gene)
pplot.get_figure().savefig("figures/test_violin_samples_CD74.png".format(gene),bbox_inches='tight')
pplot.get_figure().clf()


CD14_Monocyte = ["VCAN","FCN1","S100A8","S100A9","CD14","HLA-DRA", "HLA-DRB1"]
HSC = ["CD34", "KIT", "CD59", "THY1", "SOX4"]
Erythroid_cells = ["GYPA","TFRC", "ITGA4","HBB", "HBA1","ANK1", "ICAM4", "BCAM", "SLC4A1", "ACKR1"]
Neutrophil = ["FUT4", "MPO", "CSF3R","FCGR3A", "FCGR3B", "CEACAM8"]
Plasma = ["CD38", "XBP1", "CD27", "SLAMF7","TNFRSF17", "TNFRSF13B","IGHA1", "IGHG1"]
Cytotoxic_T_cell = ["CD3D", "CD3E", "CD3G","TRAC","CD8A", "CD8B"]
Memory_B = ["CD19","MS4A1","CD79A", "CD79B","MZB1","HLA-DRA", "HLA-DRB1","CD27"]
Stromal = ["THY1","PDGFRA","PDGFRB","MCAM","NT5E","ENG","VCAM1","CXCL12","SDF1","FGF2","CXC3CL1","CD90","FL","GAS6","GAS1","RARRES2"]
MAST_cells = ["TPSAB1", "TPSB2", "CPA3", "KIT", "FCER1G"]
reduced_stromal = ["THY1","PDGFRA","PDGFRB","VCAM1"]

heatmap_genes = ["MEIS1","EGR1","MSI2","CD34","EGFL7","CD38","MPO","ELANE","CTSG","AZU1","LYZ","CEBPD","MNDA","FCER1G","FCN1","CD14","C5AR1","CLEC4A","CLEC10A","FCER1A","CLEC4C","PTPRS","IRF8","TCF4","CSF1","KIT","HBB","HBD","GYPA","CD24","EBF1","MME","CD19","CD79A","MS4A1","BANK1","MZB1","IGLL5","JCHAIN","CD3D","CD3G","IL32","IL7R","CCL5","GZMK","CD8A","KLRB1","KLRD1","GZMB","NCAM1"]
hm_genes = [gene for gene in heatmap_genes if gene in adata.var_names]

order = ["HSC","Progenitor","MKP","GMP","Pro-Mono","Mono","pDC","cDC","Early Eryth","Late Eryth","CLP","Pro-B","B","Plasma","CTL","NK","T","Stroma"]

top_genes = get_top_n_genes_per_cluster(adata, n_genes = 5)
genes_of_interest = ["RUNX1","CD74","MIF","FOS","CCL2", "PU.1", "TLR4", "TLR2","CD44", "SRC", "PI3K","AKT", "TEN", "JAB1", "CXCL8", "MAPK", "ERK", "SPI1", "STAT3", "STAT5", "NFKb", "CXCR1/2","CXCR1","CXCR2",  "JUN", "GATA1"]
GOI = list(set([gene for genes in top_genes.values() for gene in genes] + heatmap_genes + genes_of_interest))

#heat_plot = sns.clustermap(heatmap_data, annot=False,center=0,xticklabels=True, yticklabels=True,figsize=(12,20), col_cluster = False)
		#heat_plot.ax_heatmap.set_ylabel("")
#heat_plot.cax.set_visible(False)
#heat_plot.savefig("figures/test_heatmap_velocity_scores.png",bbox_inches='tight')


groupby = "Condition"
cluster = "cluster"
data = {"condition" : [], "velocity" : [], "cluster" : [], "gene" : [], "sample" : [], "cell" : [], "Ms" : []}
for gene in GOI:
	if gene not in adata.var_names:
		continue
	else:
		dict_of_clusters = {each_cluster : [i for i,x in enumerate(adata.obs[cluster]) if x == each_cluster] for each_cluster in pd.unique(adata.obs[cluster])}
		for key,value in dict_of_clusters.items():
			n_len = len(list(adata.obs[groupby][value]))
			data["condition"] = data["condition"] + list(adata.obs[groupby][value])
			data["cluster"] = data["cluster"] + list(adata.obs[cluster][value])
			data["velocity"] = data["velocity"] + list(adata.layers['velocity'][value,adata.var_names.get_loc(gene)])	
			data["Ms"] = data["Ms"] + list(adata.layers['Ms'][value,adata.var_names.get_loc(gene)])	
			data["gene"] = data["gene"] + [gene for i in range(0,n_len)]
			#data["top_cluster"] = data["top_cluster"] + [top_key for i in range(0,n_len)]
			data["sample"] = data["sample"] + list(adata.obs["orig_ident"][value])
			data["cell"] = data["cell"] + list(adata.obs.index[value])


	
df = pd.DataFrame(data)
weights = np.ones(100)/100
df["Ms_scaled"] = np.convolve(df["Ms"], weights, mode = 'same')  
df["cluster"].cat.reorder_categories(list(order),ordered = True)
heatmap_data = pd.pivot_table(df, values='Ms_scaled', 
                              index='gene', 
                              columns=['cluster','cell'])
heatmap_data.columns = heatmap_data.columns.droplevel(1)		
df["cluster"].cat.reorder_categories(list(order),ordered = True)
pplot = sns.heatmap(heatmap_data.drop("HBB"), annot = False, yticklabels = True)	
pplot = sns.heatmap(heatmap_data, annot = False, yticklabels = True)	
pplot.get_figure().savefig("figures/test_heatmap.png",bbox_inches='tight')					  
pplot.get_figure().clf()

scv.tl.latent_time(adata,root_key = key_name)
scv.pl.heatmap(adata, var_names = heatmap_genes, save = "test_heatmap.png")
scv.pl.heatmap(adata, var_names = heatmap_genes, col_color = "cluster", col_cluster = True, n_convolve =30, sort = True, save = "test_heatmap.png", yticklabels = True, figsize = (14,10))


scv.pl.heatmap_deprecated(adata, var_names = hm_genes, groupby = "cluster", layers = "velocity", save = "test_heatmap.png")

adata.obs['clusters'] = adata.obs['cluster'].cat.reorder_categories(list(order), ordered=True)
scanpy.pl.heatmap(adata, var_names = hm_genes, swap_axes = True, var_group_positions = [(1,3)], standard_scale = "var",groupby = "clusters", layer = "Mu", save = "scanpy_Mu_heatmap.png")
scanpy.pl.heatmap(adata, var_names = hm_genes, swap_axes = True, var_group_positions = [(1,3)], standard_scale = "var",groupby = "clusters", layer = "Ms", save = "scanpy_Ms_heatmap.png")

FPD = adata[adata.obs["Condition"].isin(["FPD"])]
HD = adata[adata.obs["Condition"].isin(["HD"])]
scanpy.pl.heatmap(FPD, var_names = hm_genes, swap_axes = True, var_group_positions = [(1,3)], standard_scale = "var",groupby = "clusters", layer = "Ms", save = "scanpy_FPD_heatmap.png")
scanpy.pl.heatmap(HD, var_names = hm_genes, swap_axes = True, var_group_positions = [(1,3)], standard_scale = "var",groupby = "clusters", layer = "Ms", save = "scanpy_HD_heatmap.png")
scanpy.pl.heatmap(FPD, var_names = ordered_genes , swap_axes = True, var_group_positions = [(1,3)], standard_scale = "var",groupby = "clusters", layer = "Ms", save = "scanpy_Ms_top_genes_FPD.png", figsize = (14,10), show_gene_labels=True)
scanpy.pl.heatmap(HD, var_names = ordered_genes , swap_axes = True, var_group_positions = [(1,3)], standard_scale = "var",groupby = "clusters", layer = "Ms", save = "scanpy_Ms_top_genes_HD.png", figsize = (14,10), show_gene_labels=True)
ordered_genes = []
for name in order:
	for gene in top_genes[name]:
		ordered_genes.append(gene)
scanpy.pl.heatmap(adata, var_names = ordered_genes , swap_axes = True, var_group_positions = [(1,3)], standard_scale = "var",groupby = "clusters", layer = "Ms", save = "scanpy_Ms_top_genes_heatmap.png", figsize = (14,10), show_gene_labels=True)
scanpy.pl.heatmap(adata, var_names = ordered_genes , swap_axes = True, var_group_positions = [(1,3)], standard_scale = "var",groupby = "clusters", layer = "Mu", save = "scanpy_Mu_top_genes_heatmap.png", figsize = (14,10), show_gene_labels=True)
scv.tl.differential_kinetic_test(adata, var_names=hm_genes, groupby='cluster')


df = scv.get_df(adata[:, hm_genes], ['fit_diff_kinetics', 'fit_pval_kinetics'], precision=2)
kwargs = dict(linewidth=2, add_linfit=True, frameon=False)
scv.pl.scatter(adata, basis = df.index[:15], add_outline='fit_diff_kinetics', ncols=5, **kwargs,color = "cluster", save = "test_kinetics.png", legend_loc = "best")

diff_clusters=list(set(adata[:, df.index].var['fit_diff_kinetics']))
diff_clusters = [cluster for cluster in diff_clusters if cluster in pd.unique(adata.obs["cluster"])]
scv.pl.scatter(adata, basis = "umap",legend_loc='right', size=60, title='diff kinetics', color = "cluster" , add_outline=diff_clusters, outline_width=(.8, .2), save = "test_outline_clusters.png")

scv.tl.recover_dynamics(adata)

adata.write("scvelo_analysis_object.h5ad")
top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:100]
scv.tl.differential_kinetic_test(adata, var_names=top_genes, groupby='clusters')

scanpy.pl.clustermap(adata, obs_key = "cluster", save = "test_clustermap.png")


## check if limiting the number of clusters resolves the root key
key_clusters = ["HSC","Progenitor","MKP","GMP"]
subset_adata = adata[adata.obs["cluster"].isin(key_clusters)]

key_root = return_root(subset_adata, root_cluster = "HSC", method = "max_x")
key_name = subset_adata.obs.index[key_root]
scv.tl.velocity_pseudotime(subset_adata, root_key = key_root)
scv.pl.scatter(subset_adata, color='velocity_pseudotime', color_map='gnuplot', save = "pseudotime.png")
scv.tl.score_genes_cell_cycle(subset_adata)
subset_adata.uns['neighbors']['distances'] = subset_adata.obsp['distances']
subset_adata.uns['neighbors']['connectivities'] = subset_adata.obsp['connectivities']
scv.tl.paga(subset_adata, groups='cluster')
df = scv.get_df(subset_adata, 'paga/transitions_confidence', precision=2).T
df.to_csv("PAGA_cluster_transitions.tsv",sep="\t")

scv.pl.paga(subset_adata, basis='umap', size=50, alpha=.1,min_edge_width=2, node_size_scale=1.5, save = "test_PAGA.png")

#####################################################################################################################################
def plot_violin_compares_velocity(adata,gene, A,B , groupby = "Condition",cluster = "cluster", **kwargs):
	"""
	#args: adata  = scanpy or scvelo anndata object
	#	  gene = gene in the matrix
	#	  groupby = is the condition splitting the data
	#	  cluster = obs column of the clusters or other categorical data
	#	  A = cluster to compare with B
	#	  B = cluster to compare with B
	#	  **kwargs to be using in seaborn violinplot
	#returns nothing but saves an png image with stats comaring gene between clustsers A and B within groupby 
	#example usage:  plot_violin_compares_velocity(adata,gene = "MSI2", A = "HSC",B = "GMP",groupby = "Condition", cluster = "cluster",hue_order = Order_plot, palette = color_dict)
	"""
	groups = pd.unique(adata.obs[groupby])
	clusts = [A,B]
	samps = list(([(x,y) for x in clusts for y in groups])) #get base comparison list
	compares = [(x,y) for i,x in enumerate(samps) for j,y in enumerate(samps) if j > i and x != y] #get all comparisons to be made
	if gene not in adata.var_names:
		sys.stderr.write("{} not included in the object\n".format(gene))
		return None
	dict_of_clusters = {each_cluster : [i for i,x in enumerate(adata.obs[cluster]) if x == each_cluster] for each_cluster in [A,B]}
	data = {"condition" : [], "velocity" : [], "cluster" : []}
	for key,value in dict_of_clusters.items():
		data["condition"] = data["condition"] + list(adata.obs[groupby][value])
		data["cluster"] = data["cluster"] + list(adata.obs[cluster][value])
		data["velocity"] = data["velocity"] + list(adata.layers['velocity'][value,adata.var_names.get_loc(gene)])
	data_f = pd.DataFrame(data)
	sns.set(style="whitegrid")
	pplot = sns.violinplot(x = "cluster", y = "velocity", data = data_f,inner = "box", hue = "condition", **kwargs)
	pplot1 = add_stat_annotation(pplot, data = data_f, x='cluster', y = "velocity", hue="condition", box_pairs = compares, test = 't-test_welch', text_format = "full", loc = 'outside')
	#pplot1[0].set_xticklabels(pplot.get_xticklabels(),rotation = 90)
	pplot1[0].set_title(gene)
	path = "figures/violin/{}/{}".format(A.replace("/","."),B.replace("/","."))
	os.makedirs(os.path.join(os.getcwd(),path),exist_ok = True)
	pplot1[0].get_figure().savefig("figures/violin/{}/{}/violin_gene_{}_cluster_{}_vs_{}_plot.png".format(A.replace("/","."),B.replace("/","."),gene,A.replace("/","."),B.replace("/",".")),bbox_inches='tight')	
	pplot1[0].get_figure().clf()
	pplot.get_figure().clf()
	

for gene in GOI:
	for i,clust in enumerate(pd.unique(adata.obs["cluster"])):
		for j,clustB in enumerate(pd.unique(adata.obs["cluster"])):
			if clust == clustB or i < j:
				continue
			else:
				plot_violin_compares_velocity(adata,gene = gene, A = clust,B = clustB,groupby = "Condition", cluster = "cluster",hue_order = Order_plot, palette = color_dict)
				#print(gene,clust,clustB)

def plot_gene_velocity(adata,gene,groupby = "Condition",cluster = "cluster", **kwargs):
	"""
	#args: adata  = scanpy or scvelo anndata object
	#	  gene = gene in the matrix
	#	  groupby = is the condition splitting the data
	#	  cluster = obs column of the clusters or other categorical data
	#	  **kwargs to be using in seaborn violinplot
	#returns nothing but saves an png image
	#example usage:  plot_gene_velocity(adata,gene = "MSI2", groupby = "Condition", cluster = "cluster",hue_order = Order_plot, palette = color_dict)
	"""
	if gene not in adata.var_names:
		sys.stderr.write("{} not included in the object\n".format(gene))
		return None
	dict_of_clusters = {each_cluster : [i for i,x in enumerate(adata.obs[cluster]) if x == each_cluster] for each_cluster in pd.unique(adata.obs[cluster])}
	data = {"condition" : [], "velocity" : [], "cluster" : []}
	for key,value in dict_of_clusters.items():
		data["condition"] = data["condition"] + list(adata.obs[groupby][value])
		data["cluster"] = data["cluster"] + list(adata.obs[cluster][value])
		data["velocity"] = data["velocity"] + list(adata.layers['velocity'][value,adata.var_names.get_loc(gene)])
	data_f = pd.DataFrame(data)
	sns.set(style="whitegrid")
	pplot = sns.violinplot(x = "cluster", y = "velocity", data = data_f,inner = "box", hue = "condition", **kwargs)
	pplot.set_xticklabels(pplot.get_xticklabels(),rotation = 90)
	pplot.set_title(gene)
	#pplot.set(ylabel = "velocity score")
	pplot.get_figure().savefig("figures/violin_gene_{}_clusters_velocity_scores.png".format(gene),bbox_inches='tight')
	pplot.get_figure().clf()
			
			
def get_top_n_genes_per_cluster(adata,n_genes = 10, cluster = "cluster"):
	"""
	#args: adata  = scanpy or scvelo anndata object
	#	  n_genes = number of genes to select from
	#	  cluster = obs column of the clusters or other categorical data
	#returns a dictionary of the top n genes per cluster according to each genes velocity 
	#keys are the cluster names and values are the list of top genes
	#example usage:
	#	top_genes = get_top_n_genes_per_cluster(adata, n_genes = 10)
	"""
	if "top_cluster_velocity" not in adata.uns:
		scanpy.tl.rank_genes_groups(adata, groupby = cluster,layer = 'velocity', method = "wilcoxon", use_raw = False,key_added = "top_cluster_velocity")
	results = adata.uns['top_cluster_velocity']
	groups = results['names'].dtype.names
	keys = list(results.keys())
	keys = keys[1:]
	all_dicts = {group : pd.DataFrame({"{}".format(key) : results[key][group] for key in keys}) for group in groups}
	for df in all_dicts.keys():
		all_dicts[df][cluster] = df
	top_n_genes_cluster = {key : list(value.nlargest(n_genes,"scores").names) for key,value in all_dicts.items()}
	return top_n_genes_cluster	
	


"""
############################################## HEATMAPS


Color_hex = ["#67CAD4","#B0207E"] #purple for FPD, teal for HD
Order_plot = ["HD","FPD"]

color_dict        = dict(zip(Order_plot,Color_hex))


heatmap_genes = ["MEIS1","EGR1","MSI2","CD34","EGFL7","CD38","MPO","ELANE","CTSG","AZU1","LYZ","CEBPD","MNDA","FCER1G","FCN1","CD14","C5AR1","CLEC4A","CLEC10A","FCER1A","CLEC4C","PTPRS","IRF8","TCF4","CSF1","KIT","HBB","HBD","GYPA","CD24","EBF1","MME","CD19","CD79A","MS4A1","BANK1","MZB1","IGLL5","JCHAIN","CD3D","CD3G","IL32","IL7R","CCL5","GZMK","CD8A","KLRB1","KLRD1","GZMB","NCAM1"]
hm_genes = [gene for gene in heatmap_genes if gene in adata.var_names]

order = ["HSC","Progenitor","MKP","GMP","Pro-Mono","Mono","pDC","cDC","Early Eryth","Late Eryth","CLP","Pro-B","B","Plasma","CTL","NK","T","Stroma"]

top_genes = get_top_n_genes_per_cluster(adata, n_genes = 10)
genes_of_interest = ["RUNX1","CD74","MIF","FOS","CCL2", "PU.1", "TLR4", "TLR2","CD44", "SRC", "PI3K","AKT", "TEN", "JAB1", "CXCL8", "MAPK", "ERK", "SPI1", "STAT3", "STAT5", "NFKb", "CXCR1/2","CXCR1","CXCR2",  "JUN", "GATA1"]
GOI = list(set([gene for genes in top_genes.values() for gene in genes] + heatmap_genes + genes_of_interest))

adata.obs['clusters'] = adata.obs['cluster'].cat.reorder_categories(list(order), ordered=True)

################

##########
scanpy.pl.heatmap(adata, var_names = hm_genes, swap_axes = True, var_group_positions = [(1,3)], standard_scale = "var",groupby = "clusters", layer = "Mu", save = "scanpy_Mu_heatmap.png")
scanpy.pl.heatmap(adata, var_names = hm_genes, swap_axes = True, var_group_positions = [(1,3)], standard_scale = "var",groupby = "clusters", layer = "Ms", save = "scanpy_Ms_heatmap.png")

FPD = adata[adata.obs["Condition"].isin(["FPD"])]
HD = adata[adata.obs["Condition"].isin(["HD"])]
scanpy.pl.heatmap(FPD, var_names = hm_genes, swap_axes = True, var_group_positions = [(1,3)], standard_scale = "var",groupby = "clusters", layer = "Ms", save = "scanpy_FPD_heatmap.png")
scanpy.pl.heatmap(HD, var_names = hm_genes, swap_axes = True, var_group_positions = [(1,3)], standard_scale = "var",groupby = "clusters", layer = "Ms", save = "scanpy_HD_heatmap.png")
scanpy.pl.heatmap(FPD, var_names = ordered_genes , swap_axes = True, var_group_positions = [(1,3)], standard_scale = "var",groupby = "clusters", layer = "Ms", save = "scanpy_Ms_top_genes_FPD.png", figsize = (14,10), show_gene_labels=True)
scanpy.pl.heatmap(HD, var_names = ordered_genes , swap_axes = True, var_group_positions = [(1,3)], standard_scale = "var",groupby = "clusters", layer = "Ms", save = "scanpy_Ms_top_genes_HD.png", figsize = (14,10), show_gene_labels=True)


####################################################




############


def gene_correlation(adata,gene, layer = "velocity", n = 10):
	"""
	# args: adata = scanpy or scvelo anndata object
	#	layer = layer in data object
	#
	#
	"""
	if gene not in adata.var_names:
		return "gene not in var names"
	gene_exp = adata.layers[layer][:,adata.var_names.get_loc(gene)]
	if np.isnan(gene_exp).any():
		return "invalid gene"
	data = {"Correlation" : [],"Gene" : []}
	for i in range(0,adata.layers[layer].shape[1]):
		if np.isnan(adata.layers[layer][:,i]).any():
			continue
		elif i != adata.var_names.get_loc(gene):
			gene_cor = stats.pearsonr(gene_exp, adata.layers[layer][:,i])[0]
			data["Correlation"].append(gene_cor)
			data["Gene"].append(adata.var_names[i])
	df = pd.DataFrame(data)
	top_n = df.nlargest(n,"Correlation")
	btm_n = df.nsmallest(n,"Correlation")
	Heatmap_genes = list(top_n["Gene"])
	Heatmap_genes.append(gene)
	Heatmap_genes = Heatmap_genes + list(btm_n["Gene"])
	scanpy.pl.heatmap(adata, var_names = Heatmap_genes , swap_axes = True, var_group_positions = [(1,3)], standard_scale = "var",groupby = "clusters", layer = layer, save = "scanpy_correl_{}_{}_heatmap.png".format(gene,layer), figsize = (14,10), show_gene_labels=True)
#	results = adata.layers[layer][:,[adata.var_names.get_loc(hmg) for hmg in Heatmap_genes]]
	
for gene in GOI:
	if "heatmapscanpy_correl_{}_Ms_heatmap.png".format(gene) not in os.listdir("figures"):
		gene_correlation(adata, gene = gene, layer = "Ms", n = 10)
#scanpy.pl.heatmap(adata, var_names = Heatmap_genes , swap_axes = True, var_group_labels = ["cor_gene: {}".format(gene)],var_group_positions = [(n-1,n)], standard_scale = "var",groupby = "clusters", layer = layer, save = "scanpy_correl_heatmap.png", figsize = (14,10), show_gene_labels=True)	


###########################################################
#########           violin plots clusters     #############
###########################################################
unique_clusters = pd.unique(adata.obs["clusters"])
cluster = "clusters"
scanpy.tl.rank_genes_groups(adata, groupby = cluster, method = "wilcoxon")
RNA_results = adata.uns['rank_genes_groups']
groups = RNA_results['names'].dtype.names
keys = list(RNA_results.keys())
keys = keys[1:]
all_dicts = {group : pd.DataFrame({"{}".format(key) : RNA_results[key][group] for key in keys}) for group in groups}
for df in all_dicts.keys():
	all_dicts[df]["cluster"] = df
top_n_genes = 10	
top_n_genes_cluster = {key : list(value.nlargest(top_n_genes,"scores").names) for key,value in all_dicts.items()}

layer = "velocity"
for clust in top_n_genes_cluster.keys():
	dict_of_clusters = {each_cluster : [i for i,x in enumerate(adata.obs[cluster]) if x == each_cluster] for each_cluster in [clust]}
	data = {"condition" : [], layer : [], "cluster" : []}
	for key,value in dict_of_clusters.items():
		data["condition"] = data["condition"] + list(adata.obs[groupby][value])
		data["cluster"] = data["cluster"] + list(adata.obs[cluster][value])
		results = adata.layers[layer][value,:]
		results = np.nanmean(results[:,[adata.var_names.get_loc(gene) for gene in top_n_genes_cluster[key]]], axis = 1)
		data[layer] = data[layer] + list(results)
	data_f = pd.DataFrame(data)
	pplot = sns.violinplot(x = "condition", y = layer, data = data_f,inner = "box", palette = color_dict)
	pplot1 = add_stat_annotation(pplot, data = data_f, x='condition', y = layer, box_pairs = [tuple(Order_plot)], test = 't-test_welch', text_format = "full", loc = 'outside')
	#pplot1[0].set_xticklabels(pplot.get_xticklabels(),rotation = 90)
	pplot1[0].set_title("{}:{}".format(key,str(top_n_genes_cluster[key])))
	pplot1[0].set(ylabel = "mean({})".format(layer))
	pplot1[0].get_figure().savefig("figures/violin/cluster_{}_group_{}_genes_plot.png".format(key,layer),bbox_inches='tight')
	pplot1[0].get_figure().clf()
	pplot.get_figure().clf()

# RNA expression
groupby = "Condition"
cluster = "cluster"
layer  = "expression"
for clust in top_n_genes_cluster.keys():
	dict_of_clusters = {each_cluster : [i for i,x in enumerate(adata.obs[cluster]) if x == each_cluster] for each_cluster in [clust]}
	data = {"condition" : [], layer : [], "cluster" : []}
	for key,value in dict_of_clusters.items():
		data["condition"] = data["condition"] + list(adata.obs[groupby][value])
		data["cluster"] = data["cluster"] + list(adata.obs[cluster][value])
		results = adata.X.toarray()[value,:]
		results = np.nanmean(results[:,[adata.var_names.get_loc(gene) for gene in top_n_genes_cluster[key]]], axis = 1)
		data[layer] = data[layer] + list(results)
	data_f = pd.DataFrame(data)
	pplot = sns.violinplot(x = "condition", y = layer, data = data_f,inner = "box", palette = color_dict)
	pplot1 = add_stat_annotation(pplot, data = data_f, x='condition', y = layer, box_pairs = [tuple(Order_plot)], test = 't-test_welch', text_format = "full", loc = 'outside')
	#pplot1[0].set_xticklabels(pplot.get_xticklabels(),rotation = 90)
	pplot1[0].set_title("{}:{}".format(key,str(top_n_genes_cluster[key])))
	pplot1[0].set(ylabel = "mean({})".format(layer))
	path = "figures/violin/{}/cluster_{}_group_RNAexp_genes_plot.png".format(key,key)
	os.makedirs(os.path.join(os.getcwd(),os.path.dirname(path)),exist_ok = True)
	pplot1[0].get_figure().savefig(path,bbox_inches='tight')
	pplot1[0].get_figure().clf()
	pplot.get_figure().clf()

groupby = "Condition"
cluster = "cluster"
data = {"condition" : [], "velocity" : [], "cluster" : [], "gene" : [], "sample" : [], "Ms" : [], "Mu" : [], "cell" : []}
for clust,genes in top_n_genes_cluster.items():
	dict_of_clusters = {each_cluster : [i for i,x in enumerate(adata.obs[cluster]) if x == each_cluster] for each_cluster in [clust]}
	for key,value in dict_of_clusters.items():
		for gene in genes:
			n_len = len(list(adata.obs[groupby][value]))
			data["condition"] = data["condition"] + list(adata.obs[groupby][value])
			data["cluster"] = data["cluster"] + list(adata.obs[cluster][value])
			data["velocity"] = data["velocity"] + list(adata.layers['velocity'][value,adata.var_names.get_loc(gene)])	
			data["Ms"] = data["Ms"] + list(adata.layers['Ms'][value,adata.var_names.get_loc(gene)])
			data["Mu"] = data["Mu"] + list(adata.layers['Mu'][value,adata.var_names.get_loc(gene)])
			data["gene"] = data["gene"] + [gene for i in range(0,n_len)]
			#data["top_cluster"] = data["top_cluster"] + [top_key for i in range(0,n_len)]
			data["sample"] = data["sample"] + list(adata.obs["orig_ident"][value])
			data["cell"] = data["cell"] + list(adata.obs.index[value])
df = pd.DataFrame(data)


layer = "Mu"
for clust in pd.unique(df["cluster"]):
	results = df[df["cluster"] == clust]
	results = results[~np.isnan(results[layer])]
	pplot = sns.violinplot(x = "gene", y = layer, hue = "condition", data = results,inner = "box", hue_order = Order_plot,palette = color_dict)
	pplot.set_xticklabels(pplot.get_xticklabels(),rotation = 90)
	pplot.set_title("{}:{}".format(clust,str(top_n_genes_cluster[clust])))
	pplot.get_figure().savefig("figures/violin/violin_group_{}_{}_by_genes.png".format(clust,layer),bbox_inches='tight')
	pplot.get_figure().clf()
	
	
scv.tl.paga(adata, groups='cluster',root_key = key_name)	
scv.pl.paga(adata, basis='umap', size=50, alpha=.1,min_edge_width=2, node_size_scale=1.5, save = "test_PAGA.png")
scv.pl.paga(adata, basis = 'umap', vkey = "velocity", color = "cluster", size=50, alpha=.1,min_edge_width=2, node_size_scale=1.5, save = "test_PAGA.png")

HD = adata2[adata2.obs["Condition"].isin(["HD"])]
FPD = adata2[adata2.obs["Condition"].isin(["FPD"])]
FPD_root = return_root(FPD, root_cluster = "HSC", method = "max_x")
HD_root = return_root(HD, root_cluster = "HSC", method = "max_x")

FPD_root_name = FPD.obs.index[FPD_root]
HD_root_name = HD.obs.index[HD_root]
scv.tl.velocity_pseudotime(HD, root_key = HD_root)
scv.tl.velocity_pseudotime(FPD, root_key = FPD_root)
#scv.pl.scatter(adata, color='velocity_pseudotime', color_map='gnuplot', save = "pseudotime.png")
#scv.tl.score_genes_cell_cycle(adata)
HD.uns['neighbors']['distances'] = HD.obsp['distances']
HD.uns['neighbors']['connectivities'] = HD.obsp['connectivities']
FPD.uns['neighbors']['distances'] = FPD.obsp['distances']
FPD.uns['neighbors']['connectivities'] = FPD.obsp['connectivities']
scv.tl.paga(HD, groups='cluster',root_key = HD_root_name)
scv.tl.paga(FPD, groups='cluster',root_key = FPD_root_name)
scv.pl.paga(FPD, basis='umap', size=50, alpha=.1,min_edge_width=2, node_size_scale=1.5, save = "test_PAGA_FPD.png")
scv.pl.paga(HD, basis='umap', size=50, alpha=.1,min_edge_width=2, node_size_scale=1.5, save = "test_PAGA_HD.png")
#scv.pl.velocity_embedding_stream(adata, basis='umap', color = "Condition", save = "test_stream.png")

x, y = scv.utils.get_cell_transitions(FPD, basis='umap', starting_cell=FPD_root)
ax = scv.pl.velocity_graph(FPD, c='lightgrey', edge_width=.05, show=False)
ax = scv.pl.scatter(FPD, x=x, y=y, s=120, c='ascending', cmap='gnuplot', ax=ax)
ax.get_figure().savefig("figures/test_FPD_scatter.png",bbox_inches='tight')



#### Dynamics #####################################

df = adata.var
df = df[(df['fit_likelihood'] > .1) & df['velocity_genes'] == True]

kwargs = kwargs = dict(xscale='log', fontsize=16, show = False)

pplot = scv.pl.hist(df['fit_alpha'], xlabel='transcription rate', **kwargs)
pplot.get_figure().savefig("figures/transcription_hist.png",bbox_inches="tight")
pplot = scv.pl.hist(df['fit_beta'] * df['fit_scaling'], xlabel='splicing rate', xticks=[.1, .4, 1], **kwargs)
pplot.get_figure().savefig("figures/splicing_rate_hist.png",bbox_inches="tight")
pplot = scv.pl.hist(df['fit_gamma'], xlabel='degradation rate', xticks=[.1, .4, 1], **kwargs)
pplot.get_figure().savefig("figures/degradation_rate_hist.png",bbox_inches="tight")

#scv.tl.latent_time(adata, root_key = )
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80, save = "latent_time.png")
scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, color='velocity_pseudotime', color_map='gnuplot', save = "test_pseudotime.png")
adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']
scv.tl.paga(adata, groups='cluster', use_time_prior = "latent_time",root_key = key_name)
adata.obs["dpt_pseudotime"] = adata.obs["latent_time"]
adata.uns['iroot'] = key_root
scanpy.tl.dpt(adata)
scanpy.pl.paga_path(adata, nodes = order, keys = hm_genes, show = False, return_data = False, save = "paga_path.png")

scv.pl.paga(adata, basis='umap', size=50, alpha=.1,min_edge_width=2, node_size_scale=1.5, save = "test_PAGA_latent_time.png")

scanpy.tl.paga(adata, groups = "clusters", use_rna_velocity = True)
adata.uns['connectivities_tree'] = adata.uns["paga"]['connectivities_tree']
adata.uns['connectivities'] = adata.uns["paga"]['connectivities']
scanpy.pl.paga_adjacency(adata,save = "test_paga_adj.png")


scv.pl.scatter(adata, color=['root_cells', 'end_points','latent_time',], color_map='viridis', size=80, dpi = 300,  figsize=(7,5), save = "latent_time2.png", show=False)

scv.tl.velocity(adata, diff_kinetics = True, filter_genes = False, use_latent_time = True)
scv.tl.velocity_graph(adata, mode_neighbors = "connectivities")
scv.tl.terminal_states(adata)
scv.tl.latent_time(adata)
scv.pl.scatter(adata, color='latent_time', color_map='viridis', size=80, dpi = 300,  figsize=(7,5), save = "latent_time0.png", show=False)


####

scv.tl.differential_kinetic_test(adata, var_names = hm_genes, groupby = "clusters",add_key = "heatmap_genes")
kwargs = dict(linewidth=2, add_linfit=True, frameon=False)
scv.pl.scatter(adata, basis=hm_genes[:15], ncols=5, add_outline='heatmap_genes_diff_kinetics', **kwargs, save = "hm_diff_kinetics15-1.png",legend_loc = "best",color = "cluster")

## recomputing velocities with diff kinetics

scv.tl.velocity(adata, diff_kinetics = True, filter_genes = False, use_latent_time = True,mode='dynamical')
#tests
scv.tl.velocity_graph(adata, mode_neighbors = "connectivities", basis = "umap", tkey = "latent_time")
scv.tl.terminal_states(adata)
scv.tl.latent_time(adata)
scv.pl.scatter(adata, color=['root_cells', 'end_points','latent_time',], color_map='gnuplot', size=80, dpi = 300,  figsize=(7,5), save = "latent_time_redo.png", show=False)
scv.tl.velocity_pseudotime(adata, root_key = key_root)
scv.tl.latent_time(adata)
scv.pl.velocity_embedding(adata, dpi=120, arrow_size=2, arrow_length=2, save = "recompute_kinectis1.png")

## recomputing velocities on certain genes

# change out genes in .var["highly_variable"]
adata2 = adata
adata2.var["velocity_genes"] = [x in hm_genes for x in adata2.var["velocity_genes"].index]

scv.tl.velocity(adata2, diff_kinetics = True, filter_genes = False, use_latent_time = True,mode='dynamical')
scv.tl.velocity_graph(adata2, xkey = "Mu", gene_subset = hm_genes,mode_neighbors = "connectivities", tkey = "latent_time")
scv.tl.terminal_states(adata2)
scv.tl.latent_time(adata2)
scv.tl.velocity_pseudotime(adata2, root_key = key_root)
adata2.uns['neighbors']['distances'] = adata2.obsp['distances']
adata2.uns['neighbors']['connectivities'] = adata2.obsp['connectivities']
scv.tl.paga(adata2, groups='cluster',root_key = key_name)
scv.pl.paga(adata2, basis='umap', size=50, alpha=.1,min_edge_width=2, node_size_scale=1.5, save = "PAGA2.png")
scv.pl.heatmap(adata2, var_names = hm_genes, layer = "Mu", sort = False, col_cluster = True, save = "heatmap_test2.png")

## Draw descendents and anscetors coming from a specific cell
x, y = scv.utils.get_cell_transitions(adata, basis='umap', starting_cell=key_root)
ax = scv.pl.velocity_graph(adata, c='lightgrey', edge_width=.05, show=False)
scv.pl.scatter(adata, x=x, y=y, s=120, c='ascending', cmap='gnuplot', ax=ax, save = "draw_descendents.png")

#Do Latent time difference in HSC vs FPD
scv.pl.heatmap(adata2, var_names = hm_genes, layer = "latent_time", sort = False, col_cluster = False, save = "heatmap_test2.png")
scv.pl.heatmap(HD, var_names = hm_genes, layer = "latent_time", sort = False, col_color = "cluster",col_cluster = False, save = "HD_heatmap_test2.png")

scv.pl.heatmap(FPD, var_names = hm_genes, layer = "latent_time", sort = False, col_color = "cluster",col_cluster = False, save = "FPD_heatmap_test2.png")

## Do HSC / Progenitor cells support  FPD divide slower or faster compared to HD
clust = "HSC"

scv.tl.rank_velocity_genes(adata, groupby="cluster", min_corr=.3, n_genes = adata.n_vars, match_with = "Condition")
scv.pl.scatter(adata, basis=adata.uns['rank_velocity_genes']['names']["HSC"][:5], save = "HSC_5_velo_genes.png",legend_loc="best")
adata.uns['rank_velocity_genes']['names']

## HD
scv.tl.rank_velocity_genes(HD, groupby="cluster", min_corr=.3, n_genes = HD.n_vars)
scv.pl.scatter(HD, basis=HD.uns['rank_velocity_genes']['names']["HSC"][:5], save = "HD_HSC_5_velo_genes.png",legend_loc="best")
HD.uns['rank_velocity_genes']['names']

## FPD
scv.tl.rank_velocity_genes(FPD, groupby="cluster", min_corr=.3, n_genes = FPD.n_vars)
scv.pl.scatter(FPD, basis=FPD.uns['rank_velocity_genes']['names']["HSC"][:5], save = "FPD_HSC_5_velo_genes.png",legend_loc="best")
FPD.uns['rank_velocity_genes']['names']

#make clusers into lineages
order = ["HSC","Progenitor","MKP","GMP","Pro-Mono","Mono","pDC","cDC","Early Eryth","Late Eryth","CLP","Pro-B","B","Plasma","CTL","NK","T","Stroma"]

HSC_Progenitors = ["HSC","Progenitor"]
Lymphoid = ["CLP","NK","T","B","Plasma","CTL","Pro-B"]  #Lymphoid
Myeloid = ["GMP","Pro-Mono","Mono"]
Dendritic = ["pDC","cDC","Stroma"]
Erythroid = ["Early Eryth","Late Eryth"]
MKP = ["MKP"]
linages_sets = [HSC_Progenitors,Lymphoid,Myeloid,Dendritic,Erythroid,MKP]

lineages = {"HSC":"HSC_Progenitors", "Progenitor":"HSC_Progenitors","CLP":"Lymphoid","NK":"Lymphoid","T":"Lymphoid","B":"Lymphoid","Plasma":"Lymphoid","CTL":"Lymphoid","Pro-B":"Lymphoid",
"GMP":"Myeloid","Pro-Mono":"Myeloid","Mono":"Myeloid","pDC":"Dendritic","cDC":"Dendritic","Stroma":"Dendritic","Early Eryth":"Erythroid","Late Eryth":"Erythroid",
"MKP":"MKP"}
## top ranking velocity genes
#HD	
scv.tl.rank_velocity_genes(HD, groupby="lineages", min_corr=.3, n_genes = HD.n_vars)
scv.pl.scatter(HD, color = "lineages",basis=HD.uns['rank_velocity_genes']['names']["HSC_Progenitors"][:5], save = "HD_HSC_5_velo_genes.png",legend_loc="best")
#FPD
scv.tl.rank_velocity_genes(FPD, groupby="lineages", min_corr=.3, n_genes = FPD.n_vars)
scv.pl.scatter(FPD, color = "lineages",basis=FPD.uns['rank_velocity_genes']['names']["HSC_Progenitors"][:5], save = "FPD_HSC_5_velo_genes.png",legend_loc="best")
#top_genes["HSC"] + ["MEIS1","EGR1","MSI2","CD34"]

scv.pl.scatter(HD, color = "lineages",basis=top_genes["HSC"] + ["MEIS1","EGR1","MSI2","CD34"], save = "HD_HSC_5_velo_hmgenes.png",legend_loc="best")
scv.pl.scatter(FPD, color = "lineages",basis=top_genes["HSC"] + ["MEIS1","EGR1","MSI2","CD34"], save = "FPD_HSC_5_velo_hmgenes.png",legend_loc="best")

#Hi Anupriya. Thanks for the details and suggestions.  I think a good approach to getting at the HD vs FPD differentiation differences in the HSC and progenitor cells is to broaden the cell types/clusters into lineages.  
#rank velocity genes
def top_velo_genes(adata,name = "scores"):
	RNA_results = adata.uns['rank_velocity_genes']
	groups = RNA_results['names'].dtype.names
	#B_results = adata.uns['rank_velocity_genes']
	#B_groups = B_results['names'].dtype.names
	#groups = [group for group in groups if group in B_groups]
	keys = list(RNA_results.keys())
	keys = keys[:-1]
	all_dicts = {group : pd.DataFrame({"{}".format(key) : RNA_results[key][group] for key in keys}) for group in groups}
	for df in all_dicts.keys():
		#all_dicts[df]["cluster"] = df
		all_dicts[df].columns = ["gene",name]
	return all_dicts
	
HD_df = top_velo_genes(HD)
FPD_df = top_velo_genes(FPD)
all_df = {}
for df in HD_df:
	all_df[df] = pd.merge(HD_df[df],FPD_df[df], on = "gene", suffixes = ["HD","FPD"])

def get_velocity_table(adata):
	"""
	input: a scvelo or scanpy object
	
	Returns a merged pandas data frame of 'rank_velocity_genes'
	
	Must run 'scv.tl.rank_velocity_genes(adata)' before running
	"""
	df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
	df_scores = scv.DataFrame(adata.uns['rank_velocity_genes']['scores'])
	all_dicts = dict()
	for col in df.columns:
		all_dicts[col] = pd.DataFrame()
		all_dicts[col]["genes"] = df[col]
		all_dicts[col]["scores"] = df_scores[col]
	all_list = [value for key,value in all_dicts.items()] 
	df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['genes'],how='inner'), all_list)
	df_merged.columns = ["genes"] + list(df.columns)
	return(df_merged)