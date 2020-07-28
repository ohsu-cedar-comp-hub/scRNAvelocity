
# coding: utf-8
##### MOST COMMON ERROR:  cell barcodes do not match seurat with velocity

import velocyto as vcy
import loompy
import matplotlib.pyplot as plt
import numpy as np
from sklearn.manifold import TSNE
import os



#import sys

## consult the following github source
#https://github.com/velocyto-team/velocyto.py/blob/master/velocyto/analysis.py

#logging.basicConfig(stream=sys.stdout, format='%(asctime)s - %(levelname)s - %(message)s', level=logging.DEBUG)
#%matplotlib inline
#plt.rcParams['pdf.fonttype'] = 42


velocity_loom = snakemake.input.velocity_loom
seurat_loom = snakemake.input.seurat_loom
cluster.identy = snakemake.params.cluster
out_plot = snakemake.output.out_plot

out_dir = os.path.dirname(out_plot)

vlm = vcy.VelocytoLoom(velocity_loom)		  #velocity loom file
ds = loompy.connect(seurat_loom) #seurat object to loom ... in r 'as.loom(seuratObj, filename = "seuratObj.loom") 

umap_coord = ds.ca["umap_cell_embeddings"]		 #umap coordinates  
cluster_ID = ds.ca[cluster.identy]			#cluster id  such as "seurat_clusters" or "integrated_snn_res.0.5" found in seurat meta data
cell_ID = ds.ca["CellID"]						#Cell ID (barcodes)

v_cells = vlm.ca["CellID"]						#Cell ID (must match)

# Hardcoded to enforce that the cellbarcode IDs match the seurat objects and the merged loom file.
# setup to change the basename attached to the cellbarcodes with a colon in velocity such as:
# '[basename(_#)]:[barcode]x' to seurat format '[barcode]_#'
# 'CEL0024_A_1:TTTACGTCAACCAATCx'  -> 'TTTACGTCAACCAATC_1'
v_cells = [x.split("_")[-1] for x in v_cells] #takes only the number after the last "_" before the ":"
v_cells = [x.replace("x","_") for x in v_cells] #replace the 'x' at the end with "_"
v_cells = [x.split(":")[1] + x.split(":")[0] for x in v_cells] 
#v_cells = [x.replace("pool","") for x in v_cells]             ## no need for this; remove this
vlm.ca["CellID"] = np.asarray(v_cells) #replace the array in the velocity file

#keep cells that are in cell_ID from seurat object; filter out the rest
vlm.filter_cells(bool_array=np.in1d(vlm.ca['CellID'],cell_ID))

#dictionary of cell id to cluster id
cell_cluster_dict = {}
for i in range(len(cell_ID)):
	cell_cluster_dict[cell_ID[i]] = cluster_ID[i]

# dictionary of cell id to umap coord
cell_umap_dict = {}
for i in range(len(cell_ID)):
	cell_umap_dict[cell_ID[i]] = umap_coord[i]

#add meta data array to velocyto loom.
# every cell gets assigned an cluster
vlm.ca['cluster'] = np.array([cell_cluster_dict[i] for i in vlm.ca['CellID']])
# every cell gets assigned umap_coordinates
vlm.ca['umap'] = np.array([cell_umap_dict[i] for i in vlm.ca['CellID']])

vlm.normalize("S", size=True, log=True)

vlm.initial_Ucell_size = vlm.U.sum(axis=0)
#print(len(vlm.U[0,:]))
vlm.filter_cells(bool_array=vlm.initial_Ucell_size > np.percentile(vlm.initial_Ucell_size, 0.5)) #filter out cells with exteremely low unspliced detection
#print(len(vlm.U[0,:]))
#vlm.plot_fractions()

#make velocyto aware of cluster annotation
vlm.set_clusters(vlm.ca["cluster"])

#now using clustering annotation select the genes that are expressed above a threshold of total number of molecules in any of the clusters
vlm.score_detection_levels(min_expr_counts=40, min_cells_express=35)
vlm.filter_genes(by_detection_levels=True)

#feature selection
vlm.score_cv_vs_mean(1000, plot=True, max_expr_avg=35)
vlm.filter_genes(by_cv_vs_mean=True)


# normalize our data by size (total molecule count)
vlm._normalize_S(relative_size=vlm.S.sum(0),
			 target_size=vlm.S.sum(0).mean())
vlm._normalize_U(relative_size=vlm.U.sum(0),
			 target_size=vlm.U.sum(0).mean())



#For the preparation of the gamma fit we smooth the data using a kNN neighbors pooling approach. kNN neighbors can be calculated directly in gene expression space or reduced PCA space, using either correlation distance or euclidean distance. One example of set of parameters is provided below
vlm.perform_PCA()
vlm.knn_imputation(n_pca_dims=20, balanced=True, b_sight=1000, b_maxl=500, k=500, n_jobs=16)

# fit gamma to every gene that survived the filtering step run
vlm.fit_gammas()

#The fit can be visualized by calling plot_phase_portraits and listing the gene names
#vlm.plot_phase_portraits(["Igfbpl1", "Pdgfra"])

#calculate velocity and extrapolate the future state of the cells
vlm.predict_U()
vlm.calculate_velocity()
vlm.calculate_shift(assumption="constant_velocity")
vlm.extrapolate_cell_at_t(delta_t=1.)

#set up for projection onto UMAP

bh_tsne = TSNE()
vlm.ts = vlm.ca['umap']


cluster_list = vlm.ca["cluster"]


#Use correlation to estimate transition probabilities for every cells to its embedding neighborhood
vlm.estimate_transition_prob(hidim="Sx_sz", embed="ts", transform="sqrt", psc=1,
							 n_neighbors=2000, knn_random=True, sampled_fraction=0.5)

#Use the transition probability to project the velocity direction on the embedding							 
vlm.calculate_embedding_shift(sigma_corr = 0.1, expression_scaling=True)

#plot umap with cluster labels
plt.figure(figsize=(10,10))
vcy.scatter_viz(vlm.ts[:,0], vlm.ts[:,1], c=vlm.colorandum, s=2)
for i in range(len(cluster_list)):
	ts_m = np.median(vlm.ts[vlm.ca["cluster"] == cluster_list[i], :], 0)
	plt.text(ts_m[0], ts_m[1], str(vlm.cluster_labels[vlm.ca["cluster"] == cluster_list[i]][0]),
			 fontsize=13, bbox={"facecolor":"w", "alpha":0.6})

plt.savefig(os.path.join(out_dir,"velocity_label_plot.png"))



#Calculate the velocity using a points on a regular grid and a gaussian kernel
vlm.calculate_grid_arrows(smooth=0.5, steps=(40, 40), n_neighbors=100)
#vlm.flow = vlm.flow*50
#plt.figure(None,(20,10))
plt.figure(figsize=(10,10))
vlm.plot_grid_arrows(quiver_scale=0.1,
					scatter_kwargs_dict={"alpha":0.35, "lw":0.35, "edgecolor":"0.4", "s":38, "rasterized":True}, min_mass=1, angles='xy', scale_units='xy',
					headaxislength=2.75, headlength=5, headwidth=4.8, minlength=0.5,plot_dots = True,
					plot_random=False, scale_type = "absolute")


plt.savefig(out_plot)
vlm.to_hdf5(os.path.join(out_dir,"saved_velocyto.hdf5")) #can be reloaded with load_velocyto_hdf5(filename)

ds.close()
