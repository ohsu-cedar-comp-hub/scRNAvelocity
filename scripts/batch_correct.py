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
import plotly
from combat.pycombat import pycombat
import scipy

"""
https://www.hansenlab.org/velocity_batch

1) Construct M = S + U, as well as R = S/U

2) Batch correct M, for example by using ComBat, to get a corrected matrix Mb

3) Back transform into corrected matrices Sb, Ub by Sb = Mb * R, Ub = Mb * (1-R)

"In our experience, this approach works best on log transformed expression measures. Specifically, we correct S and U for library size, and form M=U+S. Then we log transform M as log(M+1), use ComBat and invert the log transformation. This is what we feed to scVelo and friends. Note that the sum happens on the original scale."
"""
velocity_loom = snakemake.input.loom_file
output_file = snakemake.output.out_file

adata = scv.read(velocity_loom)
adata.obs_names_make_unique("-")

non_duplicates = [x for x in adata.obs_names if "-" not in x]

adata = adata[adata.obs_names.isin(non_duplicates)]

#make gene names unique
adata.var_names_make_unique("-")
adata.obsm["X_pca"] = adata.obsm["pca_cell_embeddings"]
scv.set_figure_params('scvelo')

row_sum_splice = adata.layers["spliced"].sum(axis = 0)
row_sum_unsplice = adata.layers["unspliced"].sum(axis = 0)


spliced_indicies = list(np.where(row_sum_splice.A1 == 0)[0])
unspliced_indicies = list(np.where(row_sum_unsplice.A1 == 0)[0])


#remove rows with zero spliced counts
mask = np.ones(len(row_sum_splice.A1),dtype = bool)
mask[spliced_indicies] = False
spliced = adata.layers["spliced"][:,mask]

unspliced_indicies = list(np.where(row_sum_unsplice.A1 == 0)[0])
unspliced = adata.layers["unspliced"][:,mask]
#correct each for library size
scale_factor = 10000
unspliced_normalized = adata.layers["unspliced"]/row_sum_unsplice*scale_factor
spliced_normalize = adata.layers["spliced"]/row_sum_splice*scale_factor

spliced_normalize = np.nan_to_num(spliced_normalize,nan=0.)
unspliced_normalize = np.nan_to_num(unspliced_normalized,nan=0.)

unspliced_normalize = scipy.sparse.csr_matrix(unspliced_normalize)
spliced_normalize = scipy.sparse.csr_matrix(spliced_normalize)

M = unspliced_normalize + spliced_normalize

#locate the cells in each sample
lookup_dict = {x:i for i,x in enumerate(pd.unique(adata.obs["samples"]))}
lookup_dict_cond = {x:i for i,x in enumerate(pd.unique(adata.obs["Condition"]))}
samples_list = [lookup_dict[x] for x in adata.obs["samples"]]


#remove which genes have zero variance (all zeros>?)
row_sums = M.sum(axis = 0)
positions_zero = np.where(row_sums.A1 == 0)[0]

#remove genes that have zero variance/zero counts
indices = list(positions_zero)
mask = np.ones(M.shape[1],dtype = bool)
mask[indices] = False

M = M[:,mask]
#log normalize
M = np.log1p(M)

var_names = adata.var_names
corrected_vars = var_names[mask]

M_df = pd.DataFrame(M.transpose().todense(),columns = adata.obs_names,index = corrected_vars)



#R is ratio of spliced / unspliced

R = spliced_normalize[:,mask] / unspliced_normalize[:,mask]
R = np.nan_to_num(R,nan=0.)
R = scipy.sparse.csr_matrix(R)

#combat batch correction
df_corrected = pycombat(M_df, batch=samples_list)

#convert to sparse matrix
df_corrected = scipy.sparse.csr_matrix(df_corrected)

#inverse log1p
df_corrected = np.expm1(df_corrected)

#transpose
df_corrected = df_corrected.transpose()

#Sb = corrected_data * R
Sb = df_corrected.multiply(R)

#Ub = corrected_data * (1-R)
Ub = df_corrected.multiply((scipy.sparse.csr_matrix(np.ones(R.shape))-R))

#add to adata object
adata.layers["Ub"] = Ub
adata.layers["Sb"] = Sb
adata.layers["M"] = df_corrected

#save object
#out_object = os.path.join(os.path.dirname(velocity_loom),"batch_corrected.h5ad")
adata.write_h5ad(output_file)
