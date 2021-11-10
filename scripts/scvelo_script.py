import scanpy
import matplotlib.pyplot as plt
import numpy as np
from sklearn.manifold import TSNE
import seaborn as sns                     
import os
import sys

import scvelo as scv

curr_dir = os.getcwd()

velocity_loom = snakemake.input.velocity_loom
seurat_loom = snakemake.input.seurat_loom

out_object = snakemake.output.out_object
n_cores = snakemake.params.n_cores
out_dir = os.path.dirname(out_object)
#walkthrough
#https://colab.research.google.com/github/theislab/scvelo_notebooks/blob/master/VelocityBasics.ipynb#scrollTo=iHl8jdCUd1j8

#scvelo documentation
#https://readthedocs.org/projects/scvelo/downloads/pdf/latest/





#read loom 
adata = scv.read(velocity_loom)

#remove cell id duplicates
adata.obs_names_make_unique("-")
non_duplicates = [x for x in adata.obs_names if "-" not in x]
adata = adata[adata.obs_names.isin(non_duplicates)]

#make gene names unique
adata.var_names_make_unique("-")

os.chdir(out_dir)

#proportions of spliced/unspliced counts 
prop_plot = scv.pl.proportions(adata, show = False)

sns.set_style('whitegrid', {'font.family':'standard', 'font.serif':'Arial'})
#filter genes, normalize per cell, filter genes dispersion, and scale (log1p)
#normalize and run velocity analysis
scv.pp.filter_and_normalize(adata, n_top_genes=None, min_shared_counts = 20)
#PCA
scanpy.tl.pca(adata)
#find knn 
scanpy.pp.neighbors(adata, n_neighbors = 30, n_pcs = 30)
#first and second order moments (means and uncentered variances) computed among nearest neighbors in PCA space, computes: pca and neighbors
scv.pp.moments(adata, n_pcs = 30, n_neighbors = 30)


#default mode for velocity is stochastic,  mode = 'dynamical' and mode = "deterministic" are also available.   see https://scvelo.readthedocs.io/about.html
scv.tl.recover_dynamics(adata,n_jobs = n_cores)
scv.tl.velocity(adata,mode = 'dynamical')
#transition probabilties calculated by cosine correlation between the potential cell-to-cell transitions

#project velocities onto umap
scv.tl.velocity_graph(adata, basis = 'umap')
scv.pl.velocity_embedding_stream(adata, basis='umap', color = 'cluster', save = "scvelo_stream.png")
scv.pl.velocity_embedding(adata, basis='umap', color = 'cluster', arrow_length=3, arrow_size=2, dpi=120, save = "scvelo_embedding.png")
scv.pl.velocity_embedding_grid(adata,basis='umap', color = 'cluster', save = "scvelo_grid.png")


#get velocity confidence
scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95], save = "scatter_confidence.png")


df = adata.obs.groupby('samples')[keys].mean().T
#df.style.background_gradient(cmap='coolwarm', axis=1)
df.to_csv("velo_confidence_samples.tsv",sep="\t")

df = adata.obs.groupby('cluster')[keys].mean().T
#df.style.background_gradient(cmap='coolwarm', axis=1)
df.to_csv("velo_confidence_cluster.tsv",sep="\t")

 



#save plot proportions
fig = prop_plot.get_figure()
fig.savefig('figures/proportions.png')

#meta_data = scv.get_df(adata, keys=None, layer=None)
os.chdir(curr_dir)

adata.write_h5ad(out_object)


