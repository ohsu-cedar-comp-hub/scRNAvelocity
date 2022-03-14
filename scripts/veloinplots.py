import matplotlib.pyplot as plt
import numpy as np
import os
import datetime
import scvelo as scv
import scanpy
import sys
from velocity_fx import plot_gene_velocity,velocity_by_sample
import itertools
import random
import pandas as pd
from scipy import stats
from statannot import add_stat_annotation
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats.kde import gaussian_kde
from functools import reduce
import seaborn as sns
from pathlib import Path

in_object = snakemake.input.in_object
output_file = snakemake.output.output_file
markers_file = snakemake.input.markers_file
genes_of_interest = snakemake.params.genes

markers_dir = snakemake.params.markers_dir
out_dir     = snakemake.params.out_dir
condition = snakemake.params.seurat_status
clusters_of_interest = snakemake.params.clusters_of_interest
random.seed(123)

cur_dir = os.getcwd()
# read inputs
adata = scv.read(in_object)
markers = pd.read_csv(markers_file,sep = "\t")


if not any(["." in col for col in adata.obs.columns]):
    condition = condition.replace(".", "_")

# Global Variables
# set plot style to Arial publication ready
sns.set_style('whitegrid', {'font.family': 'standard', 'font.serif': 'Arial'})

# setup condition color for figures
Order_plot, Color_hex = None, None

if (snakemake.params.color_hex == "None"):
    Color_hex = sns.color_palette(palette=None, n_colors=len(pd.unique(adata.obs[condition])))
else:
    Color_hex = snakemake.params.color_hex
if (snakemake.params.order_plot == "None"):
    Order_plot = list(pd.unique(adata.obs[condition]))
else:
    Order_plot = snakemake.params.order_plot

color_dict = dict(zip(Order_plot, Color_hex))
adata.uns['Condition_colors'] = np.array(Color_hex)
adata.uns['Condition_color_dict'] = color_dict

Order_cluster, Cluster_hex = None, None

if (snakemake.params.cluster_hex == "None"):
    Cluster_hex = sns.color_palette(palette=None, n_colors=len(pd.unique(adata.obs['clusters'])))
else:
    Cluster_hex = snakemake.params.cluster_hex
if (snakemake.params.order_cluster == "None"):
    Order_cluster = list(pd.unique(adata.obs['clusters']))
    adata.obs['clusters'] = adata.obs['cluster']
else:
    Order_cluster = snakemake.params.order_cluster
    adata.obs['clusters'] = adata.obs['cluster'].cat.reorder_categories(list(Order_cluster), ordered=True)
cluster_colors = dict(zip(Order_cluster, Cluster_hex))

adata.uns['clusters_colors'] = np.array(Cluster_hex)
adata.uns['clusters_colors_dict'] = cluster_colors



# change output dir
os.chdir(out_dir)
violin_files = os.path.join(os.getcwd(),"violin_files")
dir_violin = os.path.join(violin_files ,"DE_between_clusters")
os.makedirs(dir_violin ,exist_ok = True)
os.chdir(dir_violin)

#specific genes:
HSC_genes = [x for x in ["MEIS1","EGR1","MSI2","CD34","PROM1","EGFL7"] if x in adata.var_names]
Prog_genes = [x for x in ["MEIS1","EGR1","MSI2","CD38","PROM1","EGFL7"] if x in adata.var_names]
MKP_genes = [x for x in ["BAP1","PRKAR2B","DDI2","P2RY14","GPR171","CEP95","LINC00152","TPP2","RGS18","FCER1G", "PTK2","VWF","FNIP1","ACTN1","PSMB4","PDLIM5", "LRRC8B"] if x in adata.var_names]
GMP_genes = [x for x in ["MPO","ELANE","CTSG","AZU1","LYST"] if x in adata.var_names]
Mono_genes =  [x for x in["LYZ","CEBPD","MNDA","FCER1G","FCN1","CD13","CSAR1","CLEC4A"] if x in adata.var_names]
top_genes = { "HSC": HSC_genes,"Progenitor" : Prog_genes,"MKP" : MKP_genes,"GMP" : GMP_genes, "Mono" : Mono_genes}

# only significant
significance = 0.005
markers = markers[markers['p_val_adj' ] <significance]
#only markers that are in var_names
markers = markers[markers['gene'].isin(adata.var_names)]
# only markers that have velocities
valid_genes = [gene for gene in list(markers['gene']) if all(~np.isnan(adata.layers["velocity"][: ,adata.var_names.get_loc(gene)]))]
markers = markers[markers['gene'].isin(valid_genes)]

Numb_genes = [5,30,50]
show = False
for i ,clust in enumerate(clusters_of_interest):
    clust_markers = markers[markers['cluster'].isin([clust])]
    # the differences between these two is that 'valid' refers to eliminating genes without velocities first to give a number of genes closer to the target n
    # whereas the others have the genes validated after picking n genes.
    for n in Numb_genes:

        submark = clust_markers.nlargest(n ,"avg_logFC")
        s_submark = clust_markers.nsmallest(n ,"avg_logFC")

        submark_genes = list(submark['gene'])
        s_submark_genes = list(s_submark['gene'])

        # get genes that are up regulated in FPD
        if submark.shape[0] > 0:
            plot_gene_velocity(adata ,gene = submark_genes, groupby = condition, cluster = "cluster", clust = clust
                               ,layer = "velocity" ,median = False ,stat = True ,show = show ,dir_base = "FPD_up"
                               ,hue_order = Order_plot, palette = color_dict)
            velocity_by_sample(adata ,genes = submark_genes, groupby = condition, cluster = "cluster", clust = clust
                               ,layer = "velocity", show = False ,dir_base = "FPD_up_by_gene", stat = True ,hue_order = Order_plot, palette = color_dict)
        # get genes that are down regulated in FPD
        if s_submark.shape[0] > 0:
            plot_gene_velocity(adata ,gene = s_submark_genes ,groupby = condition, cluster = "cluster", clust = clust
                               ,layer = "velocity" ,median = False ,stat = True ,show = show ,dir_base = "FPD_down"
                               ,hue_order = Order_plot, palette = color_dict)
            velocity_by_sample(adata ,genes = s_submark_genes, groupby = condition, cluster = "cluster", clust = clust
                               ,layer = "velocity", show = False ,dir_base = "FPD_down_by_gene", stat = True
                               ,hue_order = Order_plot, palette = color_dict)


# change output dir
next_violin_files = os.path.join(violin_files ,"canonical_genes")
os.makedirs(next_violin_files ,exist_ok = True)
os.chdir(next_violin_files)


for clust in clusters_of_interest:
    #get genes to plot
    cur_genes = top_genes[clust]
    valid_genes = [gene for gene in cur_genes if all(~np.isnan(adata.layers["velocity"][: ,adata.var_names.get_loc(gene)]))]

    plot_gene_velocity(adata,gene = valid_genes, groupby = "Condition", cluster = "cluster", clust = clust,layer = "velocity",median = False,stat = True,show = show,hue_order = Order_plot, palette = color_dict)
    velocity_by_sample(adata ,genes = valid_genes, groupby = condition, cluster = "cluster", clust = clust ,layer = "velocity", show = False ,dir_base = "by_gene", stat = True
                               ,hue_order = Order_plot, palette = color_dict)

os.chdir(cur_dir)
Path(output_file).touch()
