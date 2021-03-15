import scanpy
import matplotlib.pyplot as plt
import numpy as np
from sklearn.manifold import TSNE
import os
import sys
import datetime
import scvelo as scv

curr_dir = os.getcwd()
begin_time = datetime.datetime.now().timestamp()
velocity_loom = snakemake.input.velocity_loom
genes_of_interest = snakemake.params.genes
out_object = snakemake.output.out_object
out_dir = os.path.dirname(out_object)
#walkthrough
#https://colab.research.google.com/github/theislab/scvelo_notebooks/blob/master/VelocityBasics.ipynb#scrollTo=iHl8jdCUd1j8

#scvelo documentation
#https://readthedocs.org/projects/scvelo/downloads/pdf/latest/
condition = snakemake.params.seurat_status

#ds = loompy.connect(seurat_loom,mode = "r") #seurat object to loom ... in r 'as.loom(seuratObj, filename = "seuratObj.loom") 
adata = scv.read(velocity_loom)

#remove cell id duplicates
adata.obs_names_make_unique("-")
non_duplicates = [x for x in adata.obs_names if "-" not in x]
adata = adata[adata.obs_names.isin(non_duplicates)]

#make gene names unique
adata.var_names_make_unique("-")

os.chdir(out_dir)
#matplotlib settings to 'upgraded' images
scv.set_figure_params('scvelo')
scv.logging.print_version()
#proportions of spliced/unspliced counts 
prop_plot = scv.pl.proportions(adata, show = False)


#filter genes, normalize per cell, filter genes dispersion, and scale (log1p)
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=None)
#first and second order moments (means and uncentered variances) computed among nearest neighbors in PCA space, computes: pca and neighbors
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
#default mode for velocity is stochastic,  mode = 'dynamical' and mode = "deterministic" are also available.   see https://scvelo.readthedocs.io/about.html
scv.tl.velocity(adata)
#transition probabilties calculated by cosine correlation between the potential cell-to-cell transitions

#project velocities onto umap
scv.tl.velocity_graph(adata)
color_pal = dict(zip(adata.obs["cluster"],adata.obs["ClusterID_color"]))
scv.pl.velocity_embedding_stream(adata, basis='umap', color = 'cluster', save = "scvelo_stream.png", palette = color_pal)
scv.pl.velocity_embedding(adata, basis='umap', color = 'cluster', arrow_length=3, arrow_size=2, dpi=120, save = "scvelo_embedding.png")
scv.pl.velocity_embedding_grid(adata,basis='umap', color = 'cluster', save = "scvelo_grid.png")

# timestamp
plots_time = datetime.datetime.now().timestamp()
sys.stderr.write("finished plots: " + str(round((plots_time-begin_time)/60/60,2)) + " hours")


scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95], save = "scatter_confidence.png")


df = adata.obs.groupby('samples')[keys].mean().T
#df.style.background_gradient(cmap='coolwarm', axis=1)
df.to_csv("velo_confidence_samples.tsv",sep="\t")

df = adata.obs.groupby('cluster')[keys].mean().T
#df.style.background_gradient(cmap='coolwarm', axis=1)
df.to_csv("velo_confidence_cluster.tsv",sep="\t")

#scv.tl.velocity_pseudotime(adata)
#scv.tl.velocity_clusters(adata, match_with = "cluster")
#scv.pl.scatter(adata, color='velocity_clusters', save = "scatter_velo.png")

#scv.tl.rank_velocity_genes(adata, groupby='cluster', min_corr=.3)
#df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
#df.to_csv("rank_velocity_genes_by_cluster.tsv",sep="\t")

#scv.tl.velocity_pseudotime(adata) #, groupby = "cluster")

#scv.pl.scatter(adata, color='velocity_pseudotime', color_map='gnuplot', save = "pseudotime.png")
#scv.tl.rank_velocity_genes(adata, groupby='samples', min_corr=.3)
#df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
#df.to_csv("rank_velocity_genes_by_samples.tsv",sep="\t")

#We can visualize the velocity graph to portray all velocity-inferred cell-to-cell connections/transitions
#for i in np.arange(0.1,1,step=.1):
#		i = np.round(i,1)
#		scv.pl.velocity_graph(adata, color = "cluster", palette = color_pal, threshold=i, save = "velo_graph_analysis{}.png".format(int(i*10)))
#Further, the graph can be used to draw descendents/anscestors coming from a specified cell. Here, a pre-endocrine cell is traced to its potential fate.
#x, y = scv.utils.get_cell_transitions(adata, basis='umap', starting_cell=70)
#ax = scv.pl.velocity_graph(adata, c='lightgrey', edge_width=.05, show=False)
#ax = scv.pl.scatter(adata, x=x, y=y, s=120, c='ascending', cmap='gnuplot', ax=ax)

#

# Likelihood ratio test for differential kinetics to detect clusters/lineages that display kinetic behavior that cannot be sufficiently explained by a single model for the overall dynamics
#scv.tl.differential_kinetic_test(adata, var_names = 'velocity_genes', groupby='cluster')


for gene in genes_of_interest:
	try:
		scv.pl.velocity(adata,str(gene), dpi = 120, figsize = (7,5), color = "cluster",legend_loc = 'best',save = "scatter_gene_cluster_{}.png".format(gene))
		scv.pl.velocity(adata,str(gene), dpi = 120, figsize = (7,5), color = condition,legend_loc = 'best',save = "scatter_gene_condition_{}.png".format(gene))
		#scv.pl.velocity(adata,str(gene), dpi = 120, figsize = (7,5), color = "celltype_Condition",legend_loc = 'best',save = "scatter_gene_celltype_{}.png".format(gene))
	except:
		sys.stderr.write("{} not included".format(gene))

almost_time = datetime.datetime.now().timestamp()
sys.stderr.write("almost finished in: " + str(round((almost_time-begin_time)/60/60,2)) + " hours")

#save plot proportions
fig = prop_plot.get_figure()
fig.savefig('figures/proportions.png')

#meta_data = scv.get_df(adata, keys=None, layer=None)
os.chdir(curr_dir)

adata.write_h5ad(out_object)




#completed timestamp
end_time = datetime.datetime.now().timestamp()
sys.stderr.write("finished in: " + str(round((end_time - begin_time)/60/60,2)) + " hours")