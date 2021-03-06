import scanpy
import matplotlib.pyplot as plt
import numpy as np
from sklearn.manifold import TSNE
import os
import datetime
import scvelo as scv
import sys

curr_dir = os.getcwd()
begin_time = datetime.datetime.now().timestamp()
sys.stderr.write("beginning scvelo!")
velocity_loom = snakemake.input.velocity_loom
subset_CB = snakemake.params.subset_CB

out_object = snakemake.output.out_object
out_dir = os.path.dirname(out_object)
#walkthrough
#https://colab.research.google.com/github/theislab/scvelo_notebooks/blob/master/VelocityBasics.ipynb#scrollTo=iHl8jdCUd1j8

#scvelo documentation
#https://readthedocs.org/projects/scvelo/downloads/pdf/latest/


lookup = dict()
with open('input/paths.tsv') as file_in:
	 for i,line in enumerate(file_in):
			if i == 0:
				continue #skip header line
			else:
				location = line.split("\t")[0]
				uniq	=	str(line.split("\t")[1]).rstrip()
				basename = location.split("/")[-1]
				lookup[basename] = uniq


base_number = lookup[subset_CB]



#ds = loompy.connect(seurat_loom,mode = "r") #seurat object to loom ... in r 'as.loom(seuratObj, filename = "seuratObj.loom") 
adata = scv.read(velocity_loom)

#remove cell id duplicates
adata.obs_names_make_unique("-")

non_duplicates = [x for x in adata.obs_names if "-" not in x]

adata = adata[adata.obs_names.isin(non_duplicates)]

#make gene names unique
adata.var_names_make_unique("-")

subset_CB = [x for x in adata.obs_names if base_number == x.split("_")[1]]

#subset the data
adata = adata[adata.obs_names.isin(subset_CB)]

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
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, basis='umap', color = 'cluster', save = "scvelo_stream.png")
scv.pl.velocity_embedding(adata, basis='umap', color = 'cluster', arrow_length=3, arrow_size=2, dpi=120, save = "scvelo_embedding.png")
scv.pl.velocity_embedding_grid(adata,basis='umap', color = 'cluster', save = "scvelo_grid.png")

# timestamp
plots_time = datetime.datetime.now().timestamp()
sys.stderr.write("finished plots: " + str(round((plots_time-begin_time)/60/60,2)) + " hours")

scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95], save = "scatter_confidence.png")


df = adata.obs.groupby('cluster')[keys].mean().T
#df.style.background_gradient(cmap='coolwarm', axis=1)
df.to_csv("velo_confidence_cluster.tsv",sep="\t")


#scv.tl.velocity_pseudotime(adata)
#scv.tl.velocity_clusters(adata, match_with = "cluster")
#scv.pl.scatter(adata, color='velocity_clusters', save = "scatter_velo.png")

#scv.tl.rank_velocity_genes(adata, groupby='cluster', min_corr=.3)
#df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
#df.to_csv("rank_velocity_genes_by_cluster.tsv",sep="\t")

#genes_of_interest = ["RUNX1","EEF2","HIST1H1D","RPS29","HSP90AB1","RPL3","HIST1H1E","NCL","RPS27","RPS28","H3F3B","NUCB2","CD74","RPS21","PRDX1","RPL36","FOS","RPS4X","RPL36A","JUN","NRIP1","SLC25A6","IFI44L","RPS4Y1","HIST1H1C","HSP90AA1","IER2","RPS26","SPINK2","ACTB","RPL17","ITGA4","HIST1H4C","S100A10","RPL41","LARS"]
genes_of_interest = ["RUNX1", "CD74", "MIF", "FOS", "CCL2", "PU.1", "TLR4", "TLR2"]
for gene in genes_of_interest:
	try:
		scv.pl.velocity(adata,str(gene), dpi = 120, figsize = (7,5),legend_loc = 'best',save = "scatter_gene_{}.png".format(gene))
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

adata.obs.to_csv("scvelo_obs.tsv",sep="\t")

#completed timestamp
end_time = datetime.datetime.now().timestamp()
sys.stderr.write("finished in: " + str(round((end_time - begin_time)/60/60,2)) + " hours")