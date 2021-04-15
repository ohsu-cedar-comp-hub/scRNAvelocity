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
import cellrank as cr

def plot_gene_velocity(adata,gene,groupby = "Condition",cluster = "cluster", clust = 'HSC',layer = "velocity", dir_base = "violin_plots", stat = False, show = False,**kwargs):
    """
    # seaborn violin plot 
    # if a list of genes: plots the mean of a given layer, x= 'groupy', y = mean(layer), subset by 'clust' in 'cluster' 
    # if a single gene: plots the values of the layer,
    # returns plot if show = True or saves plots in [current working directory]/figures/[dir_base]/[clust]/  
    # 
    # args: 
    #
    # adata  = scanpy or scvelo anndata object
    # gene = gene (str) or list of genes
    # groupby = (str)is the condition splitting the data
    # cluster = (str) obs column of the clusters or other categorical data
    # clust   = (str) a specific entry that is in 'cluster'
    # dir_base = (str) specifies an additional extension for the save location 
    # stat = (bool) perform stats across the groupby (prone to breaking)
    # show = (bool) if true returns plot (default = False)
    # **kwargs to be using in seaborn violinplot
    #returns nothing but saves an png image
    # 
    #example usage:  plot_gene_velocity(adata,gene = "MSI2", groupby = "Condition", cluster = "cluster", layer = "velocity",hue_order = Order_plot, palette = color_dict)
    """
    if isinstance(gene,str):
        # if gene not in var names, return error
        if gene not in adata.var_names:
            sys.stderr.write("{} not included in the object\n".format(gene))
            return None
        #dictionary of idx list of clust in cluster
        dict_of_clusters = {each_cluster : [i for i,x in enumerate(adata.obs[cluster]) if x == each_cluster] for each_cluster in [clust]}
        data = {groupby : [], layer : []}
        # gather data into a dictionary
        for key,value in dict_of_clusters.items():
            data[groupby] = data[groupby] + list(adata.obs[groupby][value])
            data[layer] = data[layer] + list(adata.layers[layer][value,adata.var_names.get_loc(gene)])
        # convert dictionary to pandas dataframe for ease of plot making
        data_f = pd.DataFrame(data)
        sns.set(style="whitegrid")
        # plot data
        pplot = sns.violinplot(x = groupby, y = layer, data = data_f,inner = "box", **kwargs)
        # rotate xlabels, set titles, save figure
        pplot.set_xticklabels(pplot.get_xticklabels(),rotation = 50)
        pplot.set_title(gene)
        if show:
            return pplot
        pplot.get_figure().savefig("figures/violin_genes_{}_{}_{}.png".format(cluster,layer,gene),bbox_inches='tight')
        # clear figure to avoid plotting on the same figure
        pplot.get_figure().clf()
        
    elif isinstance(gene,list):
        genes = gene
        #dictionary of clusters being the keys and list of index locations being the values
        dict_of_clusters = {each_cluster : [i for i,x in enumerate(adata.obs[cluster]) if x == each_cluster] for each_cluster in [clust]}
        data = {groupby : [], layer : []}
        included_genes = []
        #determine genes that are valid (in adata object)
        for gene in genes:
            if gene not in adata.var_names:
                sys.stderr.write("{} not included in the object\n".format(gene))
                continue
            else:
                included_genes.append(gene)
        for key,value in dict_of_clusters.items():
            n_len = len(list(adata.obs[groupby][value]))
            data[groupby] = data[groupby] + list(adata.obs[groupby][value])
            #get layer of the cells to calculate
            results = adata.layers[layer][value,:] 
            #get the mean of each cell excluding not a number (NaN)
            results = np.nanmedian(results[:,[adata.var_names.get_loc(gene) for gene in included_genes]], axis = 1)
            data[layer] = data[layer] + list(results)
        data_f = pd.DataFrame(data)
        sns.set(style="whitegrid")
        pplot = sns.violinplot(x = groupby, y = layer, data = data_f,inner = "box",  **kwargs)
        if (stat):
            #add statistic for groupby. easy to break if box_pairs parameter is not satisfied properly
            pplot = add_stat_annotation(pplot, data = data_f, x=groupby, y = layer, box_pairs = [tuple(Order_plot)], test = 't-test_welch', text_format = "full", loc = 'outside')
            pplot = pplot[0]
            #raise plot title to make room for stats
            pplot.set_title("{}: top {} DE genes".format(clust,len(genes)), y = 1.2)
        else:
            pplot.set_title("{}: top {} DE genes".format(clust,len(genes)))
        pplot.set_xticklabels(pplot.get_xticklabels(),rotation = 50)
        pplot.set(ylabel = "median({})".format(layer))
        base_dir = "figures/{}/{}".format(dir_base,clust)
        os.makedirs(os.path.join(os.getcwd(),base_dir),exist_ok = True)
        if show:
            return pplot
        pplot.get_figure().savefig("{}/violin_genes_{}_{}_{}_multigene.png".format(base_dir,groupby,clust,layer),bbox_inches='tight')
        pplot.get_figure().clf()
    else:
        sys.stderr.write("{} not list nor str\n".format(gene))
        return
        
def velocity_by_sample(adata,genes,groupby = "Condition",cluster = "cluster", clust = 'HSC',layer = "velocity", stat = False, show = False,**kwargs):
    ###############################
    #  SHOW Violin plots by gene ##
    #  y axis is 'layer' of choice of adata object
    #  x axis is gene in list of genes
    #  hue is 'groupby'
    # example kwargs: kwargs = {hue_order = Order_plot,palette = color_dict}
    ###############################

    dict_of_clusters = {each_cluster : [i for i,x in enumerate(adata.obs[cluster]) if x == each_cluster] for each_cluster in [clust]}
    data = {groupby : [], layer : [], "gene": []}
    included_genes = []
    #determine genes that are valid (in adata object)
    for gene in genes:
        # disregard gene if... 
        #gene isn't in the object
        if gene not in adata.var_names:
            sys.stderr.write("{} not included in the object\n".format(gene))
            continue
        #or gene has any 'nan' entries in layer
        elif any(np.isnan(adata.layers[layer][:,adata.var_names.get_loc(gene)])):
            sys.stderr.write("{} not included in the layer: {}\n".format(gene,layer))
            continue
        else:
            #include the gene otherwise
            included_genes.append(gene)
    #loop through the cells of the cluster
    for key,value in dict_of_clusters.items():
        for gene in included_genes:
            n_len = len(list(adata.obs[groupby][value]))
            data[groupby] = data[groupby] + list(adata.obs[groupby][value])
            #get layer of the cells to calculate
            results = adata.layers[layer][value,adata.var_names.get_loc(gene)]
            data[layer] = data[layer] + list(results)
            #repeat gene for the length of entries
            data["gene"] = data["gene"] + [gene for i in range(0,n_len)]
    data_df = pd.DataFrame(data)
    if show == False:
        return data_df
    else:
        #remove any nan (shouldn't be any since we should have caught them in the elif statment above)
        data_df = data_df[~np.isnan(data_df[layer])]
        #plot 
        pplot = sns.violinplot(x = "gene", y = layer, hue = groupby, data = data_df,inner = "box", **kwargs)
        pplot.set_xticklabels(pplot.get_xticklabels(),rotation = 90)
        pplot.set_title("{}".format(clust))
        return pplot

def heatmap_velocity(adata,clust, gene, groupby = "Condition", save_path = None):
    #Function plot Heatmap of velocity by abundance
    ##by gene and cluster
    ## if save_path is None (default) then it will not save. otherwise save path should be the path to the directory to save
    #clust = "HSC"  #choose cluster
    #gene = "FOS" #top_n_genes_cluster[clust][12] #choose gene
    if np.isnan(adata.layers['velocity'][:,adata.var_names.get_loc(gene)]).any():
        return "gene {} has no velocity".format(gene)
    groupby = "Condition"
    cluster = "cluster"
    data = {"condition" : [], "velocity" : [], "cluster" : [], "gene" : [], "sample" : [], "cell": [],"gene_count": []}
    if clust is None:
        n_len = len(list(adata.obs[groupby]))
        #dict_of_clusters = {each_cluster : [i for i,x in enumerate(adata.obs[cluster]) if x == each_cluster] for each_cluster in pd.unique(adata.obs[cluster])}
        data["condition"] = data["condition"] + list(adata.obs[groupby])
        data["cluster"] = data["cluster"] + list(adata.obs[cluster])
        data["velocity"] = data["velocity"] + list(adata.layers['velocity'][:,adata.var_names.get_loc(gene)])
        data["gene"] = data["gene"] + [gene for i in range(0,n_len)]
        data["sample"] = data["sample"] + list(adata.obs["orig_ident"])
        data["cell"] = data["cell"] + list(adata.obs.index)
        data["gene_count"] = data["gene_count"] + list(np.squeeze(np.asarray(adata.X[:,adata.var_names.get_loc(gene)].todense())))
        clust = set(data["cluster"])
    else:
        dict_of_clusters = {each_cluster : [i for i,x in enumerate(adata.obs[cluster]) if x == each_cluster] for each_cluster in [clust]}
    #gather results
        for key,value in dict_of_clusters.items():
            n_len = len(list(adata.obs[groupby][value]))
            data["condition"] = data["condition"] + list(adata.obs[groupby][value])
            data["cluster"] = data["cluster"] + list(adata.obs[cluster][value])
            data["velocity"] = data["velocity"] + list(adata.layers['velocity'][value,adata.var_names.get_loc(gene)])
            data["gene"] = data["gene"] + [gene for i in range(0,n_len)]
            data["sample"] = data["sample"] + list(adata.obs["orig_ident"][value])
            data["cell"] = data["cell"] + list(adata.obs.index[value])
            data["gene_count"] = data["gene_count"] + list(np.squeeze(np.asarray(adata.X[value,adata.var_names.get_loc(gene)].todense())))
    results = pd.DataFrame(data)
    #return results
    
    results = results.dropna(axis=0) #remove NaN

    HD_res = results[results["condition"] == "HD"]
    #HD_res = HD_res.dropna(axis=0)
    FPD_res = results[results["condition"] == "FPD"]
    #FPD_res = FPD_res.dropna(axis=0)
    #genes_count = (results.shape[0]/len(pd.unique(results['cell'])))
    from scipy.stats.kde import gaussian_kde
    #ymin, ymax = min(HD_res['velocity'],FPD_res['velocity']),max(HD_res['velocity'],FPD_res['velocity'])
    #matplotlib setup
    fig, (ax1, ax2) = plt.subplots(1,2, sharey=True)
    #HD counts and velocity
    x, y  = HD_res["gene_count"], HD_res["velocity"]
    nbins = 300 
    # gaussian kernel density estimate
    try:
        k = gaussian_kde(np.vstack([x,y]))
        #scaling
        #root the length of x and y * 1j (to make it imaginary)
        #multi dimensional meshgrid
        #xi, yi = np.mgrid[x.min():x.max():x.size**0.5*1j,y.min():y.max():y.size**0.5*1j]
        ## other ways of scaling is to manually select the number of bins
        xi, yi = np.mgrid[x.min():x.max():nbins*1j,y.min():y.max():nbins*1j]
        zi = k(np.vstack([xi.flatten(), yi.flatten()]))
        #plot HD
        ax1.pcolormesh(xi, yi, zi.reshape(xi.shape))
        ax1.set_title("HD: {} cells".format(len(pd.unique(HD_res['cell']))))
        ax1.set_ylabel('velocity')
    except:
        return "gene {} has no velocity FPD singular matrix".format(gene)
    #FPD counts and velocity
    x, y  = FPD_res["gene_count"], FPD_res["velocity"]
    try:
        k = gaussian_kde(np.vstack([x,y]))
        #xi, yi = np.mgrid[x.min():x.max():x.size**0.5*1j,y.min():y.max():y.size**0.5*1j]
        xi, yi = np.mgrid[x.min():x.max():nbins*1j,y.min():y.max():nbins*1j]
        zi = k(np.vstack([xi.flatten(), yi.flatten()]))
        #plt.subplot(1,2,1)
        ax2.pcolormesh(xi, yi, zi.reshape(xi.shape))
        ax2.set_title("FPD: {} cells".format(len(pd.unique(FPD_res['cell']))))
        ax2.set_xlabel('abundance (log normalized mRNA counts)', loc = 'left')
    except:
        return "gene {} has no velocity FPD singular matrix".format(gene)
    fig.suptitle("{}: {} ".format(str(clust),gene), y = 1.02,x = 0.58)
    
    #plt.colorbar()
    if save_path is not None:
        try:
            if isinstance(clust,set):
                clust = "all_clusters"
            path_to_save = os.path.join("{}/{}/velocity_{}_heatmap.png".format(save_path,clust,gene))
            os.makedirs(os.path.dirname(path_to_save),exist_ok = True)
            fig.savefig(path_to_save,bbox_inches='tight')
        except:
            sys.stderr.write("did not save, check the save parameter be a string path")
    else:
        return fig.tight_layout()
        
def heatmap_velocity_sample(adata,clust, gene, groupby = "Condition", save_path = None):
    from scipy.stats.kde import gaussian_kde
    #
    #Function plot Heatmap of velocity by abundance
    ##by gene and 'clust' from 'cluster' obs. separated by 'groupby'
    ## if save_path is None (default) then it will not save. otherwise save path should be the path to the directory to save
    #clust = "HSC"  #choose cluster 
    #gene = "FOS" #top_n_genes_cluster[clust][12] #choose gene
    if np.isnan(adata.layers['velocity'][:,adata.var_names.get_loc(gene)]).any():
        return "gene {} has no velocity".format(gene)
    
    groupby = "Condition"
    cluster = "cluster"
    data = {"condition" : [], "velocity" : [], "cluster" : [], "gene" : [], "sample" : [], "cell": [],"gene_count": []}
    if clust is None:
        
        n_len = len(list(adata.obs[groupby]))
        #dict_of_clusters = {each_cluster : [i for i,x in enumerate(adata.obs[cluster]) if x == each_cluster] for each_cluster in pd.unique(adata.obs[cluster])}
        data["condition"] = data["condition"] + list(adata.obs[groupby])
        data["cluster"] = data["cluster"] + list(adata.obs[cluster])
        data["velocity"] = data["velocity"] + list(adata.layers['velocity'][:,adata.var_names.get_loc(gene)])
        data["gene"] = data["gene"] + [gene for i in range(0,n_len)]
        data["sample"] = data["sample"] + list(adata.obs["orig_ident"])
        data["cell"] = data["cell"] + list(adata.obs.index)
        data["gene_count"] = data["gene_count"] + list(np.squeeze(np.asarray(adata.X[:,adata.var_names.get_loc(gene)].todense())))
        clust = set(data["cluster"])
    else:
        dict_of_clusters = {each_cluster : [i for i,x in enumerate(adata.obs[cluster]) if x == each_cluster] for each_cluster in [clust]}
    #gather results
        for key,value in dict_of_clusters.items():
            n_len = len(list(adata.obs[groupby][value]))
            data["condition"] = data["condition"] + list(adata.obs[groupby][value])
            data["cluster"] = data["cluster"] + list(adata.obs[cluster][value])
            data["velocity"] = data["velocity"] + list(adata.layers['velocity'][value,adata.var_names.get_loc(gene)])
            data["gene"] = data["gene"] + [gene for i in range(0,n_len)]
            data["sample"] = data["sample"] + list(adata.obs["orig_ident"][value])
            data["cell"] = data["cell"] + list(adata.obs.index[value])
            data["gene_count"] = data["gene_count"] + list(np.squeeze(np.asarray(adata.X[value,adata.var_names.get_loc(gene)].todense())))
    results = pd.DataFrame(data)
    #return results
    
    results = results.dropna(axis=0) #remove NaN

    HD_res = results[results["condition"] == "HD"]
    #HD_res = HD_res.dropna(axis=0)
    FPD_res = results[results["condition"] == "FPD"]
    samples_HD = pd.unique(HD_res['sample'])
    samples_FPD = pd.unique(FPD_res['sample'])


    #matplotlib setup
    fig = plt.figure()
    #HD counts and velocity
    for i,samp in enumerate(samples_HD):
        samp_res = HD_res[HD_res["sample"]==samp]
        x, y  = samp_res["gene_count"], samp_res["velocity"]
        # gaussian kernel density estimate
        try:
            k = gaussian_kde(np.vstack([x,y]))
            #scaling
            #root the length of x and y * 1j (to make it imaginary)
            #multi dimensional meshgrid
            xi, yi = np.mgrid[x.min():x.max():x.size**0.5*1j,y.min():y.max():y.size**0.5*1j]

            ## other ways of scaling is to manually select the number of bins
            nbins = 300 
            xi, yi = np.mgrid[x.min():x.max():nbins*1j,y.min():y.max():nbins*1j]
            #
            zi = k(np.vstack([xi.flatten(), yi.flatten()]))
            #plot HD
            if i == 0:
                ax1 = fig.add_subplot(2,5,i+1)
                ax1.pcolormesh(xi, yi, zi.reshape(xi.shape))
                ax1.set_title("{}: {} cells".format(samp,len(pd.unique(samp_res['cell']))), fontsize = 6)
                ax1.set_ylabel('velocity')
            else:
                ax = fig.add_subplot(2,5,i+1, sharey=ax1, sharex=ax1)
                ax.pcolormesh(xi, yi, zi.reshape(xi.shape))
                ax.set_title("{}: {} cells".format(samp,len(pd.unique(samp_res['cell']))), fontsize = 6)
        except:
            continue

    for j,samp in enumerate(samples_FPD):
        samp_res = FPD_res[FPD_res["sample"]==samp]
        x, y  = samp_res["gene_count"], samp_res["velocity"]
        try:
            # gaussian kernel density estimate
            k = gaussian_kde(np.vstack([x,y]))
            #scaling
            #root the length of x and y * 1j (to make it imaginary)
            #multi dimensional meshgrid
            xi, yi = np.mgrid[x.min():x.max():x.size**0.5*1j,y.min():y.max():y.size**0.5*1j]

            ## other ways of scaling is to manually select the number of bins
            nbins = 300 
            xi, yi = np.mgrid[x.min():x.max():nbins*1j,y.min():y.max():nbins*1j]
            #
            zi = k(np.vstack([xi.flatten(), yi.flatten()]))
            #plot FPD

            if j == 0:
                ax2 = fig.add_subplot(2,5,j+1+i+len(samples_HD), sharex=ax1)
                ax2.pcolormesh(xi, yi, zi.reshape(xi.shape))
                ax2.set_title("{}: {} cells".format(samp,len(pd.unique(samp_res['cell']))), fontsize = 6)
                ax2.set_ylabel('velocity')
            else:
                ax = fig.add_subplot(2,5,j+1+i+len(samples_HD), sharey=ax2, sharex=ax1)
                ax.pcolormesh(xi, yi, zi.reshape(xi.shape))
                ax.set_title("{}: {} cells".format(samp,len(pd.unique(samp_res['cell']))), fontsize = 6)
        except:
            continue
    fig.suptitle("{}: {}\nby sample".format(clust,gene), y = 1.02,x = 0.58)
    #fig.tight_layout()
    if save_path is not None:
        try:
            if isinstance(clust,set):
                clust = "all_clusters"
            path_to_save = os.path.join("{}/{}/velocity_{}_heatmap_by_sample.png".format(save_path,clust,gene))
            os.makedirs(os.path.dirname(path_to_save),exist_ok = True)
            fig.savefig(path_to_save,bbox_inches='tight')
        except:
            sys.stderr.write("did not save, check the save parameter be a string path")
    else:
        return fig.tight_layout()        
        
        
def get_velocity_table(adata, score_type = "velocity_genes"):
    """
    input: a scvelo or scanpy object
    score_type: 'velocity_genes' or 'genes_groups' from running 'scanpy.tl.rank_genes_groups(adata) or 'scv.tl.rank_velocity_genes(adata)'
    Returns a merged pandas data frame of 'rank_velocity_genes'

    Must run 'scv.tl.rank_velocity_genes(adata)' before running 
    """
    df = scv.DataFrame(adata.uns['rank_{}'.format(score_type)]['names'])
    df_scores = scv.DataFrame(adata.uns['rank_{}'.format(score_type)]['scores'])
    all_dicts = dict()
    for col in df.columns:
        all_dicts[col] = pd.DataFrame()
        all_dicts[col]["genes"] = df[col]
        all_dicts[col]["scores"] = df_scores[col]
    all_list = [value for key,value in all_dicts.items()] 
    df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['genes'],how='inner'), all_list)
    df_merged.columns = ["genes"] + list(df.columns)
    return(df_merged)       

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

def main():
    input_path = snakemake.input.cr_object
    adata = cr.read(input_path)
    
    ### Set multiple variables

    ### color variables
    Color_hex = ["#67CAD4","#B0207E"] #purple for FPD, teal for HD
    Order_plot = ["HD","FPD"]

    color_dict        = dict(zip(Order_plot,Color_hex))

    ### number of genes to use
    n = 15

    ### DE cluster markers from seurat FindAllMarkers() ident.1 = cluster, ident.2 = cells
    # must include at minimum column names of 'gene', 'avg_logFC' and 'cluster'
    markers = pd.read_csv("/home/groups/CEDAR/enright/Runx_AML/seuratObj/markers_DE_FindAllMarkers/FindAllMarkers.tsv",sep = "\t")

    #only markers that are in var_names
    markers = markers[markers['gene'].isin(adata.var_names)]

    #make dictionary of the top n genes per cluster
    significance = 0.05
    #top n genes to use
    submark = {clust: markers[markers['cluster'] == clust].nlargest(n,"avg_logFC") for clust in pd.unique(markers["cluster"])}
    # ensure p_val_adj below significance threshold
    top_n_genes_cluster = {clust: list(submark[clust][submark[clust]['p_val_adj']<significance]['gene']) for clust in submark.keys()}

    sample_names = sorted(list(pd.unique(adata.obs["samples"])))

    # set plot style
    scv.set_figure_params('scvelo')
    
    MKP_genes = ["BAP1","PRKAR2B","DDI2","P2RY14","GPR171","CEP95","LINC00152","TPP2","RGS18","FCER1G", "PTK2","VWF","FNIP1","ACTN1","PSMB4","PDLIM5", "LRRC8B"]
    GENES = ["MEIS1", "EGR 1", "MSI2", "CD34", "PROM1", "EGFL7", "CD38",
           "MPO","ELANE","CTSG","AZU1","LYST","LYZ","CEBPD","MNDA",
           "FCER1G", "FCN1", "CD14", "C5AR1", "CLEC4A", "CLEC10A",
           "FCER1A", "CLEC4C", "PTPRS", "IRF8", "TCF4",
           "CSF1", "KIT", "HBB", "HBD", "GYPA",
           "CD24", "EBF1", "MME", "VPREB1", "PAX5", "CD19",
           "CD79A", "MS4A1", "BANK1", "MZB1", "IGLL5", "JCHAIN",
           "CD3D", "CD3G", "IL32", "IL7R", "TCF7", "CCL5",
           "GZMK", "CD8A", "KLRB1", "KLRD1", "GZMB", "NCAM1"]
    # heatmap genes
    heatmap_genes = ["MEIS1","EGR1","MSI2","CD34","EGFL7","CD38","MPO","ELANE","CTSG","AZU1","LYZ","CEBPD","MNDA","FCER1G","FCN1","CD14","C5AR1","CLEC4A","CLEC10A","FCER1A","CLEC4C","PTPRS","IRF8","TCF4","CSF1","KIT","HBB","HBD","GYPA","CD24","EBF1","MME","CD19","CD79A","MS4A1","BANK1","MZB1","IGLL5","JCHAIN","CD3D","CD3G","IL32","IL7R","CCL5","GZMK","CD8A","KLRB1","KLRD1","GZMB","NCAM1"]

    sample_names = pd.unique(adata.obs["samples"])
    sample_dict = dict(zip(sample_names,list(sns.color_palette(None, len(sample_names)))))
    # genes of interest
    genes_of_interest = ["RUNX1","CD74","MIF","FOS","CCL2", "PU.1", "TLR4", "TLR2","CD44", "SRC", "PI3K","AKT", "TEN", "JAB1", "CXCL8", "MAPK", "ERK", "SPI1", "STAT3", "STAT5", "NFKb", "CXCR1/2","CXCR1","CXCR2",  "JUN", "GATA1"]
    
    all_genes = set(MKP_genes + GENES + heatmap_genes + genes_of_interest)
    all_genes = [gene for gene in all_genes if gene in list(adata.var_names)]
    HSC_genes = [x for x in ["MEIS1","EGR1","MSI2","CD34","PROM1","EGFL7"] if x in adata.var_names]
    Prog_genes = [x for x in ["MEIS1","EGR1","MSI2","CD38","PROM1","EGFL7"] if x in adata.var_names]
    MKP_genes = [x for x in ["BAP1","PRKAR2B","DDI2","P2RY14","GPR171","CEP95","LINC00152","TPP2","RGS18","FCER1G", "PTK2","VWF","FNIP1","ACTN1","PSMB4","PDLIM5", "LRRC8B"] if x in adata.var_names]
    GMP_genes = [x for x in ["MPO","ELANE","CTSG","AZU1","LYST"] if x in adata.var_names]
    Mono_genes =  [x for x in["LYZ","CEBPD","MNDA","FCER1G","FCN1","CD13","CSAR1","CLEC4A"] if x in adata.var_names]
    top_genes = { "HSC": HSC_genes,"Progenitor" : Prog_genes,"MKP" : MKP_genes,"GMP" : GMP_genes, "Mono" : Mono_genes}
    #genes of interest
    GOI = list(set(MKP_genes + GENES + heatmap_genes + genes_of_interest))
    ## set input path
    os.chdir(os.path.dirname(input_path))
    
    
    # loop through clusters and plot the mean of the top n genes 
    path_ext = "violin_top_{}_DE_genes".format(n)
    for clust in pd.unique(markers["cluster"]):
        # subset of markers file
        #sub_mark = df_table[markers['cluster'] == clust].nlargest(n,"avg_logFC")
        # gene list made from the top ranking avg_logFC
        #gene_list = list(sub_mark[sub_mark['p_val_adj'] < 0.05]['gene']) #ensure adjusted p-value is significant
        gene_list = list(df_table.nlargest(n,clust)["genes"])
        # plot velocity by Condition
        plot_gene_velocity(adata,gene = gene_list, stat = True,clust = clust, groupby = "Condition", cluster = "cluster", layer = "velocity",dir_base = path_ext,hue_order = Order_plot, palette = color_dict, hue = "Condition")
        # plot Mu (moment of spliced) 30 PCA space by Condition
        plot_gene_velocity(adata,gene = gene_list, stat = True,clust = clust, groupby = "Condition", cluster = "cluster", layer = "Mu",dir_base = path_ext,hue_order = Order_plot, palette = color_dict, hue = "Condition")
        # plot Ms (moment of spliced) 30 PCA space by Condition
        plot_gene_velocity(adata,gene = gene_list, stat = True,clust = clust, groupby = "Condition", cluster = "cluster", layer = "Ms",dir_base = path_ext,hue_order = Order_plot, palette = color_dict, hue = "Condition")
        # plot velocity by samples
        plot_gene_velocity(adata,gene = gene_list, clust = clust, groupby = "samples", cluster = "cluster", layer = "velocity",dir_base = path_ext, order = sample_names)
    
    
    #DE test for expression data using 'wilcoxon' method
    if 'rank_genes_groups' not in adata.uns:
        scanpy.tl.rank_genes_groups(adata, groupby = 'clusters', n_genes=adata.n_vars, method = "wilcoxon")
    
    df_table = get_velocity_table(adata, score_type = "genes_groups")
    top_n_genes_cluster = {clust : list(df_table.nlargest(n,clust)["genes"]) for clust in pd.unique(adata.obs["cluster"]) }
    #make long data frame for all the genes and all clusters from from 'top_n_genes_cluster' dict including of multiple layers and metadata

    groupby = "Condition"
    cluster = "cluster"
    path_ext = "violin_top_{}_DE_genes".format(n)
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
    
    # plot the genes on the x axis and layer on the y split by "condition" for each cluster and each of the following layers
    layers = ["velocity","Ms","Mu"]
    sns.set(rc={'figure.figsize':(n,10)}) #sns.set(rc={'figure.figsize':(11.7,8.27)})
    for clust in pd.unique(df["cluster"]):
        for layer in layers:
            results = df[df["cluster"] == clust]
            results = results[~np.isnan(results[layer])]
            pplot = sns.violinplot(x = "gene", y = layer, hue = "condition", data = results,inner = "box", hue_order = Order_plot,palette = color_dict)
            pplot.set_xticklabels(pplot.get_xticklabels(),rotation = 45, fontsize=16)
            pplot.set_title("{}:{}".format(clust,str(top_n_genes_cluster[clust])))
            pplot.get_figure().savefig("figures/{}/{}/violin_group_{}_{}_by_genes.png".format(path_ext,clust,clust,layer),bbox_inches='tight')
            pplot.get_figure().clf()
    #Heatmaps
    path_file = "{}/figures/velocity_heatmaps".format(input_path)
    clusts = ["HSC","Progenitor","GMP","MKP"]
    for clust in clusts:
        for gene in all_genes:
            heatmap_velocity(adata, clust = clust, gene = gene, save_path = path_file)
            heatmap_velocity_sample(adata, clust = clust, gene = gene, save_path = path_file)
            
    ## Average Expression Score by time variable
    for clust,genes in top_genes.items():
    #sub_mark = markers[markers['cluster'] == clust].nlargest(n,"avg_logFC")  #top n by avg_log FC
    #genes = list(sub_mark[sub_mark['p_val_adj'] < 0.05]['gene'])             #ensure adjusted p-value is significant
        scanpy.tl.score_genes(adata, gene_list = genes, score_name = "{}_DE_{}".format(len(genes),clust))
    
    ## Heatmap
    
 
    hm_genes = [gene for gene in heatmap_genes if gene in adata.var_names]
    order = ["HSC","Progenitor","MKP","GMP","Pro-Mono","Mono","pDC","cDC","Early Eryth","Late Eryth","CLP","Pro-B","B","Plasma","CTL","NK","T","Stroma"]
    scanpy.pl.heatmap(FPD, var_names = ordered_genes , swap_axes = True, var_group_positions = [(1,3)], standard_scale = "var",groupby = "clusters", layer = "Ms", save = "scanpy_Ms_top_genes_FPD.png", figsize = (14,10), show_gene_labels=True)
    scanpy.pl.heatmap(HD, var_names = ordered_genes , swap_axes = True, var_group_positions = [(1,3)], standard_scale = "var",groupby = "clusters", layer = "Ms", save = "scanpy_Ms_top_genes_HD.png", figsize = (14,10), show_gene_labels=True)
       
        
if __name__ == '__main__':
    main()