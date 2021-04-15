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

from velocity_plots import get_velocity_table, top_velo_genes


def gene_correlation(adata,gene, layer = "velocity", n = 10):
	"""
    # returns a heatmap of genes correlated to gene according to values in layer
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

def main():

    curr_dir = os.getcwd()
    begin_time = datetime.datetime.now().timestamp()

    sys.stderr.write("beginning analyze")

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


    ############################################## HEATMAPS



    heatmap_genes = ["MEIS1","EGR1","MSI2","CD34","EGFL7","CD38","MPO","ELANE","CTSG","AZU1","LYZ","CEBPD","MNDA","FCER1G","FCN1","CD14","C5AR1","CLEC4A","CLEC10A","FCER1A","CLEC4C","PTPRS","IRF8","TCF4","CSF1","KIT","HBB","HBD","GYPA","CD24","EBF1","MME","CD19","CD79A","MS4A1","BANK1","MZB1","IGLL5","JCHAIN","CD3D","CD3G","IL32","IL7R","CCL5","GZMK","CD8A","KLRB1","KLRD1","GZMB","NCAM1"]
    hm_genes = [gene for gene in heatmap_genes if gene in adata.var_names]

    order = ["HSC","Progenitor","MKP","GMP","Pro-Mono","Mono","pDC","cDC","Early Eryth","Late Eryth","CLP","Pro-B","B","Plasma","CTL","NK","T","Stroma"]

    top_genes = get_top_n_genes_per_cluster(adata, n_genes = 10)
    genes_of_interest = ["RUNX1","CD74","MIF","FOS","CCL2", "PU.1", "TLR4", "TLR2","CD44", "SRC", "PI3K","AKT", "TEN", "JAB1", "CXCL8", "MAPK", "ERK", "SPI1", "STAT3", "STAT5", "NFKb", "CXCR1/2","CXCR1","CXCR2",  "JUN", "GATA1"]
    GOI = list(set([gene for genes in top_genes.values() for gene in genes] + heatmap_genes + genes_of_interest))

    adata.obs['clusters'] = adata.obs['cluster'].cat.reorder_categories(list(order), ordered=True)

    ################
    # HEATMAPS
    ##########
    scanpy.pl.heatmap(adata, var_names = hm_genes, swap_axes = True, var_group_positions = [(1,3)], standard_scale = "var",groupby = "clusters", layer = "Mu", save = "scanpy_Mu_heatmap.png")
    scanpy.pl.heatmap(adata, var_names = hm_genes, swap_axes = True, var_group_positions = [(1,3)], standard_scale = "var",groupby = "clusters", layer = "Ms", save = "scanpy_Ms_heatmap.png")

    FPD = adata[adata.obs["Condition"].isin(["FPD"])]
    HD = adata[adata.obs["Condition"].isin(["HD"])]
    scanpy.pl.heatmap(FPD, var_names = hm_genes, swap_axes = True, var_group_positions = [(1,3)], standard_scale = "var",groupby = "clusters", layer = "Ms", save = "scanpy_FPD_heatmap.png")
    scanpy.pl.heatmap(HD, var_names = hm_genes, swap_axes = True, var_group_positions = [(1,3)], standard_scale = "var",groupby = "clusters", layer = "Ms", save = "scanpy_HD_heatmap.png")
    scanpy.pl.heatmap(FPD, var_names = hm_genes , swap_axes = True, var_group_positions = [(1,3)], standard_scale = "var",groupby = "clusters", layer = "Mu", save = "scanpy_Ms_top_genes_FPD.png", figsize = (14,10), show_gene_labels=True)
    scanpy.pl.heatmap(HD, var_names = hm_genes , swap_axes = True, var_group_positions = [(1,3)], standard_scale = "var",groupby = "clusters", layer = "Mu", save = "scanpy_Ms_top_genes_HD.png", figsize = (14,10), show_gene_labels=True)


    ####################################################




    ############



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
                data["gene"] = data["gene"] + [gene for i in range(0,n_len)] #repeat gene for the length of entries
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

if __name__ == '__main__':
    main()
