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

curr_dir = os.getcwd()
begin_time = datetime.datetime.now().timestamp()

sys.stderr.write("beginning analyze!")

in_object         = snakemake.input.in_object

out_object        = snakemake.output.out_object
out_dir           = os.path.dirname(out_object)

genes_of_interest = snakemake.params.genes

adata = scv.read(in_object)

#################  CUSTOM Colors ################################

Color_hex = ["#54990F","#0F8299","#990F26","#3D0F99"]
Order_plot = ["wt_wk04","mut_wk04","wt_wk36","mut_wk36"]
color_dict = dict(zip(Order_plot,Color_hex))
cluster = snakemake.params.seurat_cluster
condition = snakemake.params.seurat_status
cond_levels = list(pd.unique(adata.obs[condition]))



stats_order = []
for i in range(len(cond_levels)):
	for j in range(i+1,len(cond_levels)):
		stats_order.append((cond_levels[i],cond_levels[j]))
###################################################################


os.chdir(out_dir)
if not os.path.exists("figures"):
	os.makedirs("figures")

#loop through the clusters
#analyze the condition
#adata.obs['cluster_Condition'] = adata.obs.apply(lambda x: "{}_{}".format(x['cluster'],x[condition]), axis=1)

scv.tl.rank_velocity_genes(adata, groupby="cluster", n_genes = adata.n_vars)

def get_velocity_table(adata):
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

#scv.tl.rank_velocity_genes(adata_subset, groupby="kit_version", n_genes = adata_subset.n_vars)
#df = get_velocity_table(adata_subset)
#df.to_csv("version_velocity_genes.tsv",sep="\t")
#scv.tl.rank_velocity_genes(ad_sub, groupby="kit_version", n_genes = ad_sub.n_vars)
#df.to_csv("top_{}_kit_velocity_genes.tsv".format(cont),sep="\t")
#scv.tl.score_genes_cell_cycle(adata_subset)
for cont in pd.unique(adata.obs[condition]):
	adata_subset = adata[adata.obs[condition].isin([str(cont)])]
	scv.tl.rank_velocity_genes(adata_subset, groupby="cluster", n_genes = adata_subset.n_vars)
	df = get_velocity_table(adata_subset)
	df.to_csv("top_{}_cluster_velocity_genes.tsv".format(cont),sep="\t")


unique_clusters = pd.unique(adata.obs["cluster"])

keys = 'velocity_length', 'velocity_confidence'


def return_root(adata, root_cluster = "HSC",location="topleft"):
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
	def top_left_return(potential_umap,max_std,ymax_value,xmin_value,xmean_value,ymean_value,i,j):
		
		putative_root = potential_umap[np.where( (xmin_value <= potential_umap[:,0]) & (potential_umap[:,0] < (xmin_value+(i*xsd_value))) & (potential_umap[:,1] <= ymax_value) & (potential_umap[:,1] > (ymax_value-(j*ysd_value))) & (potential_umap[:,0] >= (xmean_value-(xsd_value*i))) & (potential_umap[:,1] <= (ymean_value+(ysd_value*j))) )]
		if putative_root.shape[0] < min(10,int((len(potential_umap)/100)+1)):
			return top_left_return(potential_umap,max_std,ymax_value,xmin_value,xmean_value,ymean_value,i+0.1,j+0.1)
		elif putative_root.shape[0] > int(potential_umap.shape[0]/2):
			print(i)
			return np.argwhere(adata.obsm["umap"] == putative_root[int(putative_root.shape[0]/2)])[0][0]
		else:
			print(i,j)
			return np.argwhere(adata.obsm["umap"] == putative_root[int(putative_root.shape[0]/2)])[0][0]
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
		
	if type(location) == int:
		if location > len(potential_roots):
			location = len(potential_roots)-1
		if location < 0:
			location = 0
		return potential_roots[location]
	elif location == "topleft":
		return top_left_return(potential_umap,max_std,ymax_value,xmin_value,xmean_value,ymean_value,1,1)
	elif location == "center":
		return center_return(potential_umap,max_std,ymed_value,xmed_value,xmean_value,ymean_value,1,1)
	elif location == "max_x":
		return max_expression(potential_roots)
	elif location == "any":
		return(potential_roots[0])
	elif location == "random":
		return potential_roots[np.random.randint(1,len(potential_roots))-1]
	else:
		return np.argwhere(adata.obsm["umap"] == potential_roots[0])[0][0]

key_root = return_root(adata, root_cluster = "HSC", location = "random")
scv.tl.velocity_pseudotime(adata, root_key = key_root)
scv.pl.scatter(adata, color='velocity_pseudotime', color_map='gnuplot', save = "pseudotime.png")

adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']
scv.tl.paga(adata, groups='cluster', root_key = 'HSC')
df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
df.to_csv("PAGA_cluster_transitions.tsv",sep="\t")
color_pal = dict(zip(adata.obs["cluster"],adata.obs["ClusterID_color"]))
scv.pl.paga(adata, basis='umap', size=50, alpha=.1,min_edge_width=2, node_size_scale=1.5, save = "PAGA.png", palette = color_pal)

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
	genes_of_interest = list(set(genes_of_interest) | more_genes | less_genes )
	for gene in genes_of_interest:
		path_file = "figures/scvelo_scatter_gene_cluster_{}_{}.png".format(clust.replace("/","."),gene)
		if not os.path.isfile(path_file):
			try:
				scv.pl.velocity(adata_subset,str(gene), dpi = 120, figsize = (7,5), color = condition,legend_loc = 'best',palette = color_dict,save = "scatter_gene_cluster_{}_{}.png".format(clust.replace("/","."),gene))
			except:
				sys.stderr.write("{} not included in {}\n".format(gene,clust.replace("/",".")))
		else:
			sys.stderr.write("{} in {} already plotted \n".format(gene,clust.replace("/",".")))
		sys.stderr.write("\n")
	try:
		scv.pl.scatter(adata_subset, c=keys, cmap='coolwarm', perc=[5, 95], save = "scatter_confidence_{}.png".format(clust.replace("/",".")))
		scv.pl.velocity_embedding_stream(adata_subset, basis='umap', color = condition, palette = color_dict,save = "scvelo_stream_{}.png".format(clust.replace("/",".")))
		#scv.pl.velocity_embedding(adata_subset, basis='umap', color = 'Condition', palette = color_dict,arrow_length=0, arrow_size=0, dpi=120, save = "scvelo_embedding_{}.png".format(str(clust)))
		scanpy.pl.violin(adata_subset, keys = "velocity_length", groupby = condition, save = "velocity_length_{}.png".format(clust.replace("/",".")), palette = color_dict)
		scanpy.pl.violin(adata_subset, keys = "velocity_length", groupby = "orig_ident", save = "orig_velocity_length_{}.png".format(clust.replace("/",".")), rotation = 90, palette = color_dict)
	except:
		sys.stderr.write("Error in one of the plots\n")
	
	Count = adata_subset.obs.groupby(condition)['velocity_length'].count()
	Max = adata_subset.obs["velocity_length"].max()
	Means = adata_subset.obs.groupby(condition)["velocity_length"].mean()
	Median = adata_subset.obs.groupby(condition)["velocity_length"].median()
	p_plot = scanpy.pl.violin(adata_subset, keys = "velocity_length", groupby = condition, show = False, inner = "box", size = 0, order = cond_levels, palette = color_dict)
	
	fig = add_stat_annotation(p_plot, data = adata_subset.obs, x=condition, y = "velocity_length", box_pairs = stats_order, test = 't-test_welch', text_format = "full", loc = 'outside')
	
	add_fig = fig[0]
	for i,x in enumerate(Count):
		add_fig.text(cond_levels.index(Count.index.values[i]),Max+1,"{}".format(x))
	add_fig.get_figure().savefig("figures/stats_velocity_length_{}.png".format(str(clust).replace("/",".")),bbox_inches = "tight")
	add_fig.get_figure().clf()
	### additional testing specific to this case
	valid_var = [np.isnan(adata_subset.layers['velocity'][:,x]).any() for x,i in enumerate(adata_subset.var_names)] #any velocities in which genes have NAN ?  False means there is no NAN and has a velocities for the cells in that feature
	valid_gene_names = [adata_subset.var_names[i] for i,x in enumerate(valid_var) if x == False]
	condish = list(adata_subset.obs[condition])
	df2 = {"sum_values" : [], "mean_values" : [], "condition" : []}
	valid_velocities = adata_subset.layers['velocity'][:,[i for i,x in enumerate(valid_var) if x == False]] # array of cells -> array of velocities of valid genes
	condition_dict = {}
	for j,y in enumerate(condish): #for y in in conditions
		#append the mean/sum of the velocities of all valid genes. valid genes being ones that have velocities / spliced/unspliced
		df2["sum_values"].append(sum(valid_velocities[i,:]))
		df2["mean_values"].append(np.mean(valid_velocities[i,:]))
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
	pplot = sns.violinplot(x = 'condition', y = 'sum_values', data = df, inner = "box",order = cond_levels)
	fig2 = add_stat_annotation(pplot, data = df, x='condition', y = "sum_values", box_pairs = stats_order, test = 't-test_welch', text_format = "full", loc = 'outside')
	fig2[0].get_figure().savefig("figures/violin_test.png")
	fig2[0].get_figure().clf()
	pplot.get_figure().clf()
	
	pplot = sns.violinplot(x = 'condition', y = 'mean_values', data = df, inner = "box",order = cond_levels)
	fig2 = add_stat_annotation(pplot, data = df, x='condition', y = "sum_values", box_pairs = stats_order, test = 't-test_welch', text_format = "full", loc = 'outside')
	fig2[0].get_figure().savefig("figures/violin_mean_test.png")
	fig2[0].get_figure().clf()
	pplot.get_figure().clf()
	
	df3 = pd.DataFrame(cond_mean_genes)
	try:
		stat,p = stats.f_oneway(valid_velocities[condition_dict[list(condition_dict.keys())[0]],:],valid_velocities[condition_dict[list(condition_dict.keys())[1]],:])
		df3['anova'] = stat
		df3['p_val'] = p
	except:
		sys.stderr.write("Error in Anova\n")

	# all genes 
	df3.to_csv("anova_{}_mean_velocity_genes.tsv".format(clust.replace("/",".")), sep="\t")
	
	###
	

os.chdir(curr_dir)
adata.obs.to_csv(out_object,sep="\t")