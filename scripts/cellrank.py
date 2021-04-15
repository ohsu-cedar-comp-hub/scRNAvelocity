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
#from statannot import add_stat_annotation
from functools import reduce
import seaborn as sns
import statsmodels.api as sm
import cellrank as cr
import plotly

import mpi4py





def main():
    velocity_adata = snakemake.input.velocity_adata
    adata = scv.read(velocity_adata)

    #### project specific variables ####
    heatmap_genes = ["MEIS1","EGR1","MSI2","CD34","EGFL7","CD38","MPO","ELANE","CTSG","AZU1","LYZ","CEBPD","MNDA","FCER1G","FCN1","CD14","C5AR1","CLEC4A","CLEC10A","FCER1A","CLEC4C","PTPRS","IRF8","TCF4","CSF1","KIT","HBB","HBD","GYPA","CD24","EBF1","MME","CD19","CD79A","MS4A1","BANK1","MZB1","IGLL5","JCHAIN","CD3D","CD3G","IL32","IL7R","CCL5","GZMK","CD8A","KLRB1","KLRD1","GZMB","NCAM1"]
    # color for project
    Color_hex = ["#67CAD4","#B0207E"] #purple for FPD, teal for HD
    Order_plot = ["HD","FPD"]

    color_dict  = dict(zip(Order_plot,Color_hex))

    hm_genes = [gene for gene in heatmap_genes if gene in adata.var_names]

    order = ["HSC","Progenitor","MKP","GMP","Pro-Mono","Mono","pDC","cDC","Early Eryth","Late Eryth","CLP","Pro-B","B","Plasma","CTL","NK","T","Stroma"]
    adata.obs['cluster'] = adata.obs['cluster'].astype('category')
    adata.obs['clusters'] = adata.obs['cluster'].cat.reorder_categories(list(order), ordered=True)
    lineages = {"HSC":"HSC_Progenitors", "Progenitor":"HSC_Progenitors","CLP":"Lymphoid","NK":"Lymphoid","T":"Lymphoid","B":"Lymphoid","Plasma":"Lymphoid","CTL":"Lymphoid","Pro-B":"Lymphoid",
    "GMP":"Myeloid","Pro-Mono":"Myeloid","Mono":"Myeloid","pDC":"Dendritic","cDC":"Dendritic","Stroma":"Dendritic","Early Eryth":"Erythroid","Late Eryth":"Erythroid",
    "MKP":"MKP"}
    adata.obs["lineages"] = [lineages[key] for key in adata.obs["cluster"]]
    adata.obs["cluster_cond"] = ["{}_{}".format(x,y) for x,y in zip(adata.obs["cluster"],adata.obs["Condition"])] #combine cluster and conditions
    
    
    ###determine terminal and initial states
    cr.tl.terminal_states(adata, cluster_key='clusters', n_states = 12, weight_connectivities=0.2)
    
    cr.tl.initial_states(adata, cluster_key='clusters')
    ## determine lineages
    cr.tl.lineages(adata)
    #recover latent time
    scv.tl.recover_latent_time(adata, root_key='initial_states_probs', end_key='terminal_states_probs')
    #determine lineage drivers
    cr.tl.lineage_drivers(adata, cluster_key = "cluster", seed = 0)
    
    #diffusion pseudotime
    root_idx = np.where(adata.obs['cluster'] == 'HSC')[0][0]

    scanpy.tl.diffmap(adata)
    adata.uns['iroot'] = root_idx
    scanpy.tl.dpt(adata)
    # save object
    adata.write_h5ad(snakemake.output.output_file)
    
if __name__ == '__main__':
    main()