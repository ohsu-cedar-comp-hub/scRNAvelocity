import matplotlib.pyplot as plt
import numpy as np
from sklearn.manifold import TSNE
import os
import datetime
import scvelo as scv
import scanpy
import sys
import pandas as pd

curr_dir = os.getcwd()
begin_time = datetime.datetime.now().timestamp()

sys.stderr.write("beginning analyze!")

in_object         = snakemake.input.in_object

out_object        = snakemake.output.out_object
out_dir           = os.path.dirname(out_object)

genes_of_interest = snakemake.params.genes

#cluster = snakemake.params.cluster
#condition = snakemake.params.condition

adata = scv.read(in_object)
os.chdir(out_dir)

#loop through the clusters
#analyze the condition

unique_clusters = pd.unique(adata.obs["cluster"])

#genes_of_interest = ["RUNX1","CD74","MIF","FOS","CCL2", "PU.1", "TLR4", "TLR2","CD44", "SRC", "PI3K","AKT", "TEN", "JAB1", "CXCL8", "MAPK", "ERK", "SPI1", "STAT3", "STAT5", "NFKb", "CXCR1/2","CXCR1","CXCR2",  "JUN", "GATA1"]
keys = 'velocity_length', 'velocity_confidence'

for clust in unique_clusters:
	adata_subset = adata[adata.obs["cluster"].isin([str(clust)])]
	for gene in genes_of_interest:
		try:
			scv.pl.velocity(adata_subset,str(gene), dpi = 120, figsize = (7,5), color = "Condition",legend_loc = 'best',save = "scatter_gene_cluster_{}_{}.png".format(str(clust),gene))
		except:
			sys.stderr.write("{} not included in {}, ".format(gene,str(clust)))
	sys.stderr.write("\n")
	try:
		scv.pl.scatter(adata_subset, c=keys, cmap='coolwarm', perc=[5, 95], save = "scatter_confidence_{}.png".format(str(clust)))
		scv.pl.velocity_embedding_stream(adata_subset, basis='umap', color = 'Condition', save = "scvelo_stream_{}.png".format(str(clust)))
		scv.pl.velocity_embedding(adata_subset, basis='umap', color = 'Condition', arrow_length=3, arrow_size=2, dpi=120, save = "scvelo_embedding_{}.png".format(str(clust)))
		scanpy.pl.violin(adata_subset, keys = "velocity_length", groupby = "Condition", save = "velocity_length_{}.png".format(str(clust)))
	except:
		sys.stderr.write("Error in one of the plots\n")
	

os.chdir(curr_dir)
adata.obs.to_csv(out_object,sep="\t")