import scanpy
import matplotlib.pyplot as plt
import numpy as np
from sklearn.manifold import TSNE
import os
import datetime
import scvelo as scv
import sys
import loompy
import time


curr_dir = os.getcwd()
begin_time = datetime.datetime.now().timestamp()
sys.stderr.write("beginning scvelo!")
seurat_loom = snakemake.input.seurat_loom
sample_loom = snakemake.input.subset_CB
subset_CB = snakemake.params.subset_CB
genes_of_interest = snakemake.params.genes
out_object = snakemake.output.out_object
out_dir = os.path.dirname(out_object)
index = snakemake.params.indices
cluster_identity = snakemake.params.seurat_cluster
cluster_identity = cluster_identity.replace(".","_") #remove any periods and replace with underscores
sample_id = snakemake.params.seurat_batch
sample_id = sample_id.replace(".","_") #remove any periods and replace with underscores
status_id = snakemake.params.seurat_status
status_id = status_id.replace(".","_") #remove any periods and replace with underscores
seurat_connection = snakemake.params.seurat_cb_correction 
#walkthrough
#https://colab.research.google.com/github/theislab/scvelo_notebooks/blob/master/VelocityBasics.ipynb#scrollTo=iHl8jdCUd1j8

#scvelo documentation
#https://readthedocs.org/projects/scvelo/downloads/pdf/latest/


curr_cwd = os.getcwd()
output_file = os.path.join(out_dir,"sample.loom")
if not os.path.exists(output_file):

	#load loom objects
	ds = loompy.connect(seurat_loom, mode ='r') #seurat object to loom ... in r 'as.loom(seuratObj, filename = "seuratObj.loom") 
	merged = loompy.connect(sample_loom, mode = 'r')



	#seurat loom meta data
	seurat_cells = ds.ca["CellID"]						#Cell ID (barcodes)
	umap_coord = ds.ca["umap_cell_embeddings"]		 #umap coordinates  
	cluster_ID = ds.ca[cluster_identity]			#cluster id  such as "seurat_clusters" or "integrated_snn_res.0.5" -> "integrated_snn_res_0_5"
	sample_ids = ds.ca[sample_id]
	status_ids = ds.ca[status_id]


	#make copy of merged loom file
	view = merged.view[:,:]
	#close and tidy

	merged.close()

	#make corrected cell barcodes to match seurat objects
	#view.ca['CellID'] = np.array([s.split(":")[1].replace("x","_") + lookup[s.split(":")[0]] for s in view.ca['CellID']])
	view.ca['CellID'] = np.array([s.split(":")[1].replace("x","_") + index for s in view.ca['CellID']])

	## added to correct for funky seurat cells
	seurat_cells = np.array([s.replace(seurat_connection,"_") for s in seurat_cells])
	seurat_cells = np.array([cell for cell in seurat_cells if cell.split("_")[1] == index])
	#filter to keep seurat cells
	view = view.view[:, np.isin(view.ca.CellID, seurat_cells)]

	#add all keys
	for ca in ds.ca.keys():
		if ca == "CellID":
			continue
		curr_dict = {}
		for i,cb in enumerate(seurat_cells):
			curr_dict[cb] = ds.ca[ca][i]
		view.ca[ca] = np.array([curr_dict[j] for j in view.ca['CellID']])
		


	cell_cluster_dict = {}
	for i in range(len(seurat_cells)):
		cell_cluster_dict[seurat_cells[i]] = cluster_ID[i]

	# dictionary of cell id to umap coord
	cell_umap_dict = {}
	for i in range(len(seurat_cells)):
		cell_umap_dict[seurat_cells[i]] = umap_coord[i]

	# dictionary of cell id to umap coord
	cell_sample_dict = {}
	for i in range(len(seurat_cells)):
		cell_sample_dict[seurat_cells[i]] = sample_ids[i]
		
		
	# dictionary of cells to be split by a condition such as WT and MUT
	cell_status_dict = {}
	for i in range(len(seurat_cells)):
		cell_status_dict[seurat_cells[i]] = status_ids[i]	
		
	ds.close()
	#add meta data array to velocyto loom.
	# every cell gets assigned an cluster
	view.ca['cluster'] = np.array([cell_cluster_dict[cb] for cb in view.ca['CellID']])
	# every cell gets assigned umap_coordinates
	view.ca['umap_seurat'] = np.array([cell_umap_dict[cb] for cb in view.ca['CellID']])
	# every cell gets assigned a sample name
	view.ca['samples'] = np.array([cell_sample_dict[cb] for cb in view.ca['CellID']])
	# every cell gets assigned a status name
	view.ca['status'] = np.array([cell_status_dict[cb] for cb in view.ca['CellID']])

	#ds = loompy.connect(seurat_loom,mode = "r") #seurat object to loom ... in r 'as.loom(seuratObj, filename = "seuratObj.loom") 
	#adata = scv.read(velocity_loom)
	#create filtered loom object

	loompy.create(output_file,view.layers,view.ra,view.ca)



	counter = 0
	kill_it = 100
	while not os.path.exists(output_file):
		time.sleep(10)
		counter += 1
		if counter > kill_it:
			break

	if not os.path.isfile(output_file):
		raise ValueError("{} isn't a file!  Didn't save properly".format(output_file))

adata = scv.read(output_file)
#remove cell id duplicates
adata.obs_names_make_unique("-")

non_duplicates = [x for x in adata.obs_names if "-" not in x]

adata = adata[adata.obs_names.isin(non_duplicates)]

#make gene names unique
adata.var_names_make_unique("-")

subset_CB = [x for x in adata.obs_names if index == x.split("_")[1]]

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
scv.pp.moments(adata, n_neighbors=30,n_pcs=30)
color_pal = dict(zip(adata.obs["cluster"],adata.obs["ClusterID_color"]))

## run paga to get initial position for umap embedding
scanpy.tl.paga(adata,groups = 'cluster')
scv.pl.paga(adata, layout = "fr", size=50, alpha=.1,min_edge_width=2, node_size_scale=1.5, save = "PAGA_fr.png", palette = color_pal)
scv.tl.umap(adata, init_pos = 'paga')

#default mode for velocity is stochastic,  mode = 'dynamical' and mode = "deterministic" are also available.   see https://scvelo.readthedocs.io/about.html
scv.tl.velocity(adata)
#transition probabilties calculated by cosine correlation between the potential cell-to-cell transitions
scv.tl.velocity_graph(adata, basis = 'umap')

scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'



df = adata.obs.groupby('cluster')[keys].mean().T
#df.style.background_gradient(cmap='coolwarm', axis=1)
df.to_csv("velo_confidence_cluster.tsv",sep="\t")

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

	
key_root = return_root(adata, root_cluster = "HSC", location = "any")
scv.tl.velocity_pseudotime(adata, root_key = key_root)
scv.pl.scatter(adata, color='velocity_pseudotime', color_map='gnuplot', save = "pseudotime.png")

adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']
scv.tl.paga(adata, groups='cluster', root_key = 'HSC')
#scv.pl.paga(adata, layout = "fr", size=50, alpha=.1,min_edge_width=2, node_size_scale=1.5, save = "PAGA_fr.png", palette = color_pal)

scv.pl.paga(adata, basis = "umap", size=50, alpha=.1,min_edge_width=2, node_size_scale=1, save = "PAGA_umap.png", palette = color_pal)
df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
df.to_csv("PAGA_cluster_transitions.tsv",sep="\t")
scv.pl.velocity_embedding_stream(adata, basis='umap', color = 'cluster', save = "scvelo_stream.png", palette = color_pal)
scv.pl.velocity_embedding(adata, basis='umap', color = 'cluster', arrow_length=3, arrow_size=2, dpi=120, save = "scvelo_embedding.png")
scv.pl.velocity_embedding_grid(adata,basis='umap', color = 'cluster', save = "scvelo_grid.png")
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95], save = "scatter_confidence.png")
# timestamp
plots_time = datetime.datetime.now().timestamp()
sys.stderr.write("finished plots: " + str(round((plots_time-begin_time)/60/60,2)) + " hours")




#scv.tl.velocity_pseudotime(adata)
#scv.tl.velocity_clusters(adata, match_with = "cluster")
#scv.pl.scatter(adata, color='velocity_clusters', save = "scatter_velo.png")

#scv.tl.rank_velocity_genes(adata, groupby='cluster', min_corr=.3)
#df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
#df.to_csv("rank_velocity_genes_by_cluster.tsv",sep="\t")


for gene in genes_of_interest:
	try:
		scv.pl.velocity(adata,str(gene), dpi = 120, figsize = (7,5),color = 'cluster',legend_loc = 'best',save = "scatter_gene_{}.png".format(gene))
	except:
		sys.stderr.write("{} not included".format(gene))

almost_time = datetime.datetime.now().timestamp()
sys.stderr.write("almost finished in: " + str(round((almost_time-begin_time)/60/60,2)) + " hours")

#save plot proportions
fig = prop_plot.get_figure()
fig.savefig('figures/proportions.png')

#meta_data = scv.get_df(adata, keys=None, layer=None)

#output_file = os.path.join(out_dir,"sample.loom")
adata.write_h5ad(os.path.basename(out_object))

adata.obs.to_csv("scvelo_obs.tsv",sep="\t")

#completed timestamp
end_time = datetime.datetime.now().timestamp()
sys.stderr.write("finished in: " + str(round((end_time - begin_time)/60/60,2)) + " hours")