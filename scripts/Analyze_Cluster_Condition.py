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
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats.kde import gaussian_kde
from functools import reduce
import seaborn as sns
import statsmodels.api as sm
from scipy.stats import pearsonr
from scipy import spatial
from velocity_fx import umap_plot,get_top_n_genes_per_cluster,cos_sim,add_entropy,scatter_velocity,heatmap_velocity,plot_gene_velocity,get_velocity_table,velocity_by_sample
import itertools
import random

random.seed(123)


def main():

    # set parameters and load objects
    curr_dir = os.getcwd()
    begin_time = datetime.datetime.now().timestamp()
    
    sys.stderr.write("beginning analyze\n")
    sys.stderr.write("loading parameters and variables\n")
    in_object         = snakemake.input.in_object
    out_object        = snakemake.output.out_object
    

    adata_out         = snakemake.output.adata_out
    genes_of_interest = snakemake.params.genes
    
    condition = snakemake.params.seurat_status
    clusters_of_interest = snakemake.params.clusters_of_interest

    
    markers_dir = snakemake.params.markers_dir
    #if out directory doesn't exist:
    out_dir           = snakemake.params.out_dir
    
    if not os.path.exists(out_dir):
        sys.stderr.write("creating out path: {}".format(out_dir))
        os.mkdir(out_dir)
    
    adata = scv.read(in_object)
    
    
    if not any(["." in col for col in adata.obs.columns]):
        condition = condition.replace(".","_")
    
    #Global Variables
    # set plot style to Arial publication ready
    sns.set_style('whitegrid', {'font.family':'standard', 'font.serif':'Arial'})
    
    # setup condition color for figures
    Order_plot,Color_hex = None,None
    
    if (snakemake.params.color_hex == "None"):
        Color_hex = sns.color_palette(palette= None,n_colors = len(pd.unique(adata.obs[condition])))
    else:
        Color_hex         = snakemake.params.color_hex    
    if (snakemake.params.order_plot == "None"):
        Order_plot = list(pd.unique(adata.obs[condition]))
    else:
        Order_plot        = snakemake.params.order_plot

    color_dict        = dict(zip(Order_plot,Color_hex))    
    adata.uns['Condition_colors'] = np.array(Color_hex)
    adata.uns['Condition_color_dict'] = color_dict
    sys.stderr.write("\n color hex: {}\n".format(color_dict))
    
    # set up cluster color for figures
    adata.obs['cluster'] = adata.obs['cluster'].astype('category')
    
    
    Order_cluster,Cluster_hex = None,None
    
    if (snakemake.params.cluster_hex == "None"):
        Cluster_hex = sns.color_palette(palette= None,n_colors = len(pd.unique(adata.obs['clusters'])))
    else:
        Cluster_hex         = snakemake.params.cluster_hex    
    if (snakemake.params.order_cluster == "None"):
        Order_cluster = list(pd.unique(adata.obs['clusters']))
        adata.obs['clusters'] = adata.obs['cluster']
    else:
        Order_cluster        = snakemake.params.order_cluster
        adata.obs['clusters'] = adata.obs['cluster'].cat.reorder_categories(list(Order_cluster), ordered=True)
    cluster_colors        = dict(zip(Order_cluster,Cluster_hex))
    sys.stderr.write("\n cluster hex: {}\n".format(cluster_colors))
    
    adata.uns['clusters_colors'] = np.array(Cluster_hex)
    adata.uns['clusters_colors_dict'] = cluster_colors
    os.chdir(out_dir)
    #sample names
    sample_names = sorted(pd.unique(adata.obs['samples']))[::-1]
    #make directory for putting 'tsv' files
    #analysis
    input_des = os.path.splitext(os.path.basename(in_object))[0]
    #out_file_path = os.path.join(os.getcwd(),"{}_analysis".format(input_des))
    out_file_path = os.path.join(os.getcwd())
    
    tsv_file_path = os.path.join(os.getcwd(),"tsv_files")
    sys.stderr.write("creating tsv files path: {}".format(tsv_file_path))
    os.makedirs(tsv_file_path,exist_ok = True)
    heatmap_files = os.path.join(os.getcwd(),"heatmap_files")
    os.makedirs(heatmap_files,exist_ok = True)
    violin_files = os.path.join(os.getcwd(),"violin_files")
    os.makedirs(violin_files,exist_ok = True)
    pseudotime_files = os.path.join(os.getcwd(),"pseudotime_files")
    os.makedirs(pseudotime_files,exist_ok = True)
    Velo_heat_files = os.path.join(os.getcwd(),"Velo_heat_files")
    os.makedirs(Velo_heat_files,exist_ok = True)
    
 
    os.chdir(out_file_path)
    
    #get the unique clusters
    unique_clusters = pd.unique(adata.obs["clusters"])
    
    if clusters_of_interest == 'None':
        clusters_of_interest = unique_clusters
    
    n = 15

    sys.stderr.write("\n Applied condition:{}\n".format(condition))
    
    #get genes of interest and make sure they are in the object
    genes_of_interest = [gene for gene in genes_of_interest if gene in adata.var_names]
    
    #color for phase plot    
    phase_colors = {'G1':'red', 'G2M':'blue', 'S':'green'}        
    
    #specific genes:
    HSC_genes = [x for x in ["MEIS1","EGR1","MSI2","CD34","PROM1","EGFL7"] if x in adata.var_names]
    Prog_genes = [x for x in ["MEIS1","EGR1","MSI2","CD38","PROM1","EGFL7"] if x in adata.var_names]
    MKP_genes = [x for x in ["BAP1","PRKAR2B","DDI2","P2RY14","GPR171","CEP95","LINC00152","TPP2","RGS18","FCER1G", "PTK2","VWF","FNIP1","ACTN1","PSMB4","PDLIM5", "LRRC8B"] if x in adata.var_names]
    GMP_genes = [x for x in ["MPO","ELANE","CTSG","AZU1","LYST"] if x in adata.var_names]
    Mono_genes =  [x for x in["LYZ","CEBPD","MNDA","FCER1G","FCN1","CD13","CSAR1","CLEC4A"] if x in adata.var_names]
    top_genes = { "HSC": HSC_genes,"Progenitor" : Prog_genes,"MKP" : MKP_genes,"GMP" : GMP_genes, "Mono" : Mono_genes}
    
    #expand dictionary for all clusters
    blank_dict = {clust:[] for clust in unique_clusters}
    for clust,value in top_genes.items():
        blank_dict[clust] = value 
    top_genes = blank_dict
    ### Wilcoxon rank 
    file_path ="{}/cluster_wilcoxon_velocity.tsv".format(tsv_file_path)
    if not os.path.exists(file_path):
        #DE test for expression data using 'wilcoxon' method
        if 'rank_genes_groups' not in adata.uns:
            scanpy.tl.rank_genes_groups(adata, groupby = 'cluster', n_genes=adata.n_vars, method = "wilcoxon")
        df_table = get_velocity_table(adata, groupby = 'cluster')
        df_table.to_csv(file_path,sep="\t")
    else:
        df_table = pd.read_csv(file_path, sep='\t', index_col = 0)
        
    top_n_genes_cluster = {clust : list(set(list(df_table.nlargest(n,clust)["genes"]) + top_genes[clust])) for clust in unique_clusters }  
    
    ######
    # scvelo plots
    #####
    scv.pl.velocity_embedding_stream(adata, basis='umap', color = 'cluster', save = "scvelo_stream.png")
    scv.pl.velocity_embedding(adata, basis='umap', color = 'cluster', arrow_length=3, arrow_size=2, dpi=120, save = "scvelo_embedding.png")
    scv.pl.velocity_embedding_grid(adata,basis='umap', color = 'cluster', save = "scvelo_grid.png")
    keys = 'velocity_length', 'velocity_confidence'
    scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95], save = "scatter_confidence.png")
    
    
    
    ############################################## 
    # HEATMAPS SETUP
    ##############################################
    os.chdir(heatmap_files)
    heatmap_genes = ["MEIS1","EGR1","MSI2","CD34","EGFL7","CD38","MPO","ELANE","CTSG","AZU1","LYZ","CEBPD","MNDA","FCER1G","FCN1","CD14","C5AR1","CLEC4A","CLEC10A","FCER1A","CLEC4C","PTPRS","IRF8","TCF4","CSF1","KIT","HBB","HBD","GYPA","CD24","EBF1","MME","CD19","CD79A","MS4A1","BANK1","MZB1","IGLL5","JCHAIN","CD3D","CD3G","IL32","IL7R","CCL5","GZMK","CD8A","KLRB1","KLRD1","GZMB","NCAM1"]
    hm_genes = [gene for gene in heatmap_genes if gene in adata.var_names]
    top_genes = [gene for cluster,gene in top_n_genes_cluster.items()] #list of genes to plot for heatmap
    top_genes = [gene for group in top_genes for gene in group] #flatten list
    
    ################
    #   HEATMAPS   #
    ################ 
    scanpy.pl.heatmap(adata, var_names = top_genes, swap_axes = True, var_group_positions = [(1,3)], standard_scale = "var",groupby = "clusters", layer = "Mu", save = "scanpy_Mu_heatmap.pdf")
    scanpy.pl.heatmap(adata, var_names = top_genes, swap_axes = True, var_group_positions = [(1,3)], standard_scale = "var",groupby = "clusters", layer = "Ms", save = "scanpy_Ms_heatmap.pdf")
    
    # heatmap genes
    
    cmap_pal = "viridis"#sns.color_palette("viridis", as_cmap=True)
    scanpy.pl.heatmap(adata, var_names = hm_genes, swap_axes = True, standard_scale = "var",groupby = "clusters", layer = "Mu", cmap= cmap_pal,save = "scanpy_Mu_heatmap.pdf")
    scanpy.pl.heatmap(adata, var_names = hm_genes, swap_axes = True, standard_scale = "var",groupby = "clusters", layer = "Ms", cmap= cmap_pal,save = "scanpy_Ms_heatmap.pdf")
    # by condition
    for cond in pd.unique(adata.obs[condition]):
        cur_cond = adata[adata.obs[condition].isin([cond])]

        scanpy.pl.heatmap(cur_cond, var_names = hm_genes , swap_axes = True,  standard_scale = "var",groupby = "clusters", layer = "Mu", save = "scanpy_Mu_top_genes_{}.pdf".format(cond), cmap= cmap_pal,figsize = (14,10), show_gene_labels=True)
        scanpy.pl.heatmap(cur_cond, var_names = hm_genes , swap_axes = True,  standard_scale = "var",groupby = "clusters", layer = "Ms", save = "scanpy_Ms_top_genes_{}.pdf".format(cond), cmap= cmap_pal,figsize = (14,10), show_gene_labels=True)
    

    os.chdir(out_file_path)
    # latent time
    scv.tl.latent_time(adata)
    scv.pl.scatter(adata, color = 'latent_time',legend_loc = 'on data', color_map = "magma", save = "latent_time.pdf")
    
    
    
    root_indices = {}
    end_indices = {}
    # find terminal states per cluster and per condition
    for clust in clusters_of_interest:
        subset_adata =  adata[adata.obs['clusters'].isin([clust])]
        for cond in pd.unique(adata.obs[condition]):
            cur_cond = subset_adata[subset_adata.obs[condition].isin([cond])]
            scv.pp.neighbors(cur_cond)
            scv.tl.terminal_states(cur_cond)
            x = cur_cond.obs['root_cells']
            m_list = list(cur_cond.obs.index[np.concatenate(np.where(x == x.max())).tolist()])
            root_indices["{}_{}".format(cond,clust.replace(".","_").replace("-","_").replace(" ","_"))] = m_list
            x = cur_cond.obs['end_points']
            m_list = list(cur_cond.obs.index[np.concatenate(np.where(x == x.max())).tolist()])
            end_indices["{}_{}".format(cond,clust.replace(".","_").replace("-","_").replace(" ","_"))] = m_list
            
    # pdt_pseudotime    
    scanpy.tl.diffmap(adata)
    # run diffusion pseudotime for each root cell in each cluster of interest
    # average the root cell
    run_avg = {}
    for key,value in root_indices.items():
        dpt_cur = {}
        for cell_name in value:
            #look up cell name index
            root_idx = adata.obs_names.get_loc(cell_name)
            adata.uns['iroot'] = root_idx
            scanpy.tl.dpt(adata)
            dpt_cur[cell_name] = adata.obs['dpt_pseudotime']
        #assign the dpt_pseudotime column to be equal to the mean pseudotime starting from all HSC cells    
        mean_all = np.mean([np.array(values) for key,values in dpt_cur.items()], axis = 0)    
        adata.obs['{}_dpt'.format(key)] = mean_all
        run_avg[key] = mean_all
    #averages of the averages for each starting root    
    x = np.mean([np.array(values) for key,values in run_avg.items()], axis = 0)
    normalized = (x-min(x))/(max(x)-min(x))
    adata.obs['dpt_pseudotime_mean'] = normalized
    
    #traditional method for dpt_pseudotime
    adata.uns['iroot'] = np.where((adata.obs['cluster'] == 'HSC'))[0][0]
    scanpy.tl.dpt(adata)
    
    os.chdir(pseudotime_files)
    #latent time and pseudotime plots
    scv.pl.scatter(adata, color=['dpt_pseudotime','latent_time'], cmap='gnuplot', save = "pseudo_latent_time.pdf")
    
    #umap dpt pseudotime
    for cond in pd.unique(adata.obs[condition]):
        cur_cond = adata[adata.obs[condition].isin([cond])]
        umap_plot(cur_cond,color = 'dpt_pseudotime_mean', color_map = 'magma', title = cond, save = True, remove_color_bar = False)
        umap_plot(cur_cond,color = 'dpt_pseudotime_mean', color_map = 'magma', title = cond, save = True, remove_color_bar = True)
        scv.pl.scatter(cur_cond, color = ['dpt_pseudotime_mean'],legend_loc = None, color_map = "magma", colorbar = False,title = cond, save = "{}_dpt_pseudotime_rescaled_nobar.pdf".format(cond),rescale_color = [0,1])
        scv.pl.scatter(cur_cond, color = ['dpt_pseudotime_mean'],legend_loc = None, color_map = "magma", colorbar = True ,title = cond, save = "{}_dpt_pseudotime_rescaled_bar.pdf".format(cond),rescale_color = [0,1])
        umap_plot(cur_cond,color = 'dpt_pseudotime', color_map = 'magma', title = cond, save = True, remove_color_bar = False)
        umap_plot(cur_cond,color = 'dpt_pseudotime', color_map = 'magma', title = cond, save = True, remove_color_bar = True)
        scv.pl.scatter(cur_cond, color = ['dpt_pseudotime'],legend_loc = None, color_map = "magma", colorbar = False,title = cond, save = "{}_dpt_pseudotime_HSCroot_nobar.pdf".format(cond),rescale_color = [0,1])
        scv.pl.scatter(cur_cond, color = ['dpt_pseudotime'],legend_loc = None, color_map = "magma", colorbar = True ,title = cond, save = "{}_dpt_pseudotime_HSCroot_bar.pdf".format(cond),rescale_color = [0,1])
    
    os.chdir(out_file_path)
    #PAGA
    #directed graph 
    adata.uns['neighbors']['distances'] = adata.obsp['distances']
    adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']
    scv.tl.paga(adata, groups='clusters')
    scv.pl.paga(adata, basis='umap', size=50, alpha=.1,min_edge_width=2, node_size_scale=1.5, save = "PAGA.pdf")
    

    

    
           
    if condition == "None":
        os.chdir(curr_dir)
        sys.stderr.write("not starting analyses since there is no status to analyze! check your config file\n")
        adata.obs.to_csv(out_object,sep="\t") 
        
        adata.write_h5ad(adata_out)
        adata.obs.to_csv(out_object,sep="\t")
    else:
        #combine cluster and condition into one meta data column
        adata.obs["cluster_cond"] = ["{}_{}".format(x,y) for x,y in zip(adata.obs["cluster"],adata.obs[condition])]
        os.chdir(curr_dir)
            
        #save objects
        adata.write_h5ad(adata_out)
        adata.obs.to_csv(out_object,sep="\t")
        os.chdir(out_file_path)
        file_path ="{}/cluster_condition_velocity.tsv".format(tsv_file_path)
        if not os.path.exists(file_path):
            dataf = get_velocity_table(adata,groupby = 'cluster_cond')
            dataf.to_csv(file_path,sep="\t")
        
        
        #cmap magna
        cmap = "magma"
        #loop through each cluster
        sys.stderr.write("entering outer loop for cluster\n")
        for clust in unique_clusters:
            sys.stderr.write("Analyzing {}".format(clust))
            #make a dir for clust in figures and analyses
            os.makedirs(os.path.join(os.getcwd(),"figures/scvelo_{}".format(str(clust).replace("/","."))),exist_ok = True)
            os.makedirs(os.path.join(tsv_file_path,"{}".format(str(clust).replace("/","."))),exist_ok=True)
            os.makedirs(os.path.join(os.getcwd(),"figures/{}".format(str(clust).replace("/","."))),exist_ok = True)
            os.makedirs(os.path.join(violin_files,str(clust).replace("/",".")),exist_ok = True)
            #subset the data to current cluster
            adata_subset = adata[adata.obs["cluster"].isin([str(clust)])]
            
            #for each gene of interest make plots
            for gene in genes_of_interest:
                sys.stderr.write("On marker: {}".format(gene))
                
                #phase plot
                sys.stderr.write("\n attempting phase plot\n")
                os.chdir(out_file_path)
                scv.pl.velocity(adata_subset,str(gene), dpi = 120, figsize = (7,5), color = condition,palette = color_dict, legend_loc = 'best',save = "{}/scatter_gene_{}_cluster.pdf".format(str(clust).replace("/","."),gene))
                
                #heatmap velocity plot
                os.chdir(Velo_heat_files)
                sys.stderr.write("\n attempting heatmap velocity\n")
                heatmap_velocity(adata_subset, clust = str(clust), groupby = condition, sample = condition,gene = gene, nbins = 300, plot_max = True, order = Order_plot,save_path = "./figures" ,cmap = cmap)
                
                #scatter velocity plot
                # velocity or 'layer' on y-axis separated by groupby on x-axis 
                #sys.stderr.write("\n attempting scatter velocity plot\n")
                #add_args = { 'palette': color_dict}
                
                #if any(~np.isnan(adata_subset.layers["velocity"][:,adata_subset.var_names.get_loc(gene)])):
                #    os.chdir(violin_files)
                #    plot_gene_velocity(adata,gene = str(gene), groupby = condition, cluster = "cluster", clust = str(clust), layer = "velocity", stat = True, show = False,Order_plot =Order_plot,dir_base = "single_marker",hue_order = Order_plot)
            
            
            #plots for cluster
            sys.stderr.write("\n attempting scatter confidence\n")
            os.chdir(out_file_path)
            
            scv.pl.scatter(adata_subset, color=['velocity_length', 'velocity_confidence'], cmap='coolwarm', perc=[5, 95], save = "scatter_confidence_{}.pdf".format(str(clust).replace("/",".")))
            
            scv.pl.velocity_embedding_stream(adata_subset, basis='umap', color = condition, palette = color_dict,save = "scvelo_stream_{}.pdf".format(str(clust).replace("/",".")))
            
            scanpy.pl.violin(adata_subset, keys = "velocity_length", groupby = condition, save = "velocity_length_{}.pdf".format(str(clust).replace("/",".")), palette = color_dict)
                
            sys.stderr.write("\n attempting velocity length\n")        
            # calculate velocity length (not really interesting but why not)
            Count = adata_subset.obs.groupby(condition)['velocity_length'].count()
            Max = adata_subset.obs["velocity_length"].max()
            Means = adata_subset.obs.groupby(condition)["velocity_length"].mean()
            Median = adata_subset.obs.groupby(condition)["velocity_length"].median()
            
            
            p_plot = scanpy.pl.violin(adata_subset, keys = "velocity_length", groupby = condition, show = False, inner = "box", size = 0, order = Order_plot, palette = color_dict)
        
            fig = add_stat_annotation(p_plot, data = adata_subset.obs, x=condition, y = "velocity_length", box_pairs = [tuple(Order_plot)], test = 't-test_welch', text_format = "full", loc = 'outside')
        
            add_fig = fig[0]
            for i,x in enumerate(Count):
                add_fig.text(Order_plot.index(Count.index.values[i]),Max+1,"{}".format(x))
            add_fig.get_figure().savefig("{}/{}/stats_velocity_length_{}.pdf".format(violin_files,str(clust).replace("/","."),str(clust).replace("/",".")))
    
           
            #RANK velocity genes
            # tsv files
            # DE on velocity expression values
            sys.stderr.write("\n attempting DE velocity genes\n")
            scv.tl.rank_velocity_genes(adata_subset, groupby=condition, min_corr=.3,n_genes = adata_subset.n_vars)
            
            df = get_velocity_table(adata_subset, groupby = condition)
            #save in a table
            df.to_csv("{}/{}/top_velocity_{}_genes.tsv".format(tsv_file_path,clust.replace("/","."),clust.replace("/",".")), sep="\t")
            
            
   
            ### additional testing specific to two conditions (ie. wild-type vs mutant)
            #adata_subset.layers['velocity'] is the velocity residual for each gene per cell. axis 0 = gene, axis 1 = cell
            valid_var = [np.isnan(adata_subset.layers['velocity'][:,x]).any() for x,i in enumerate(adata_subset.var_names)] #any velocities in which genes have NAN ?  False means there is no NAN and has a velocities for the cells in that feature
            valid_gene_names = [adata_subset.var_names[i] for i,x in enumerate(valid_var) if x == False]
            valid_index = [i for i,x in enumerate(valid_var) if x == False]
            valid_gene_index = dict(zip(valid_gene_names,valid_index)) 
            #list of conditions
            condish = list(adata_subset.obs[condition])
            #make a dictionary to hold values
            df2 = {"sum_values" : [], "mean_values" : [], "condition" : []}
            valid_velocities = adata_subset.layers['velocity'][:,[i for i,x in enumerate(valid_var) if x == False]] # -> [array of cells : array of velocities of valid genes]

            condition_dict = {}
            #fill in the dictionary
            for j,y in enumerate(condish): #for y in in conditions
                #append the mean/sum of the velocities of all valid genes. valid genes being ones that have velocities / spliced/unspliced
                df2["sum_values"].append(sum(valid_velocities[j,:]))
                df2["mean_values"].append(np.mean(valid_velocities[j,:]))
                df2["condition"].append(y)
                if y not in condition_dict:
                    condition_dict[y] = [j]
                else:
                    condition_dict[y].append(j)
            cond_mean_genes = {}
            
            for key,value in condition_dict.items():
                cond_mean_genes[key] = np.mean(valid_velocities[value,:], axis = 0)  #mean of each genes velocity across the  cells of this cluster
            cond_mean_genes["gene"] = valid_gene_names
            #turn ditioncary into a dataframe
            df = pd.DataFrame(df2)
            pplot = sns.violinplot(x = 'condition', y = 'sum_values', data = df, inner = "box")
            fig2 = add_stat_annotation(pplot, data = df, x='condition', y = "sum_values", box_pairs = [tuple(Order_plot)], test = 't-test_welch', text_format = "full", loc = 'outside')
            fig2[0].get_figure().savefig("{}/{}/violin_sum_{}_plot.pdf".format(violin_files,str(clust).replace("/","."),str(clust).replace("/",".")))
            fig2[0].get_figure().clf()
            pplot.get_figure().clf()
            
            pplot = sns.violinplot(x = 'condition', y = 'mean_values', data = df, inner = "box", palette = color_dict,order = Order_plot)
            fig2 = add_stat_annotation(pplot, data = df, x='condition', y = "sum_values", box_pairs = [tuple(Order_plot)], test = 't-test_welch', text_format = "full", loc = 'outside')
            fig2[0].get_figure().savefig("{}/{}/violin_mean_{}_plot.pdf".format(violin_files,str(clust).replace("/","."),str(clust).replace("/",".")))
            fig2[0].get_figure().clf()
            pplot.get_figure().clf()
            
            df3 = pd.DataFrame(cond_mean_genes)

            try:
                stat,p = stats.f_oneway(valid_velocities[condition_dict[list(condition_dict.keys())[0]],:],valid_velocities[condition_dict[list(condition_dict.keys())[1]],:])
                df3['anova'] = stat
                df3['p_val'] = p
                df3['p_val_Adjusted'] = sm.stats.multipletests(p,method = "fdr_bh")[1]
                
            except:
                sys.stderr.write("Error in Anova\n")
            
            #wilcoxon rank-sum test on velocity values
            scanpy.tl.rank_genes_groups(adata_subset, groupby = condition,layer = 'velocity', method = "wilcoxon", use_raw = False)
            results = adata_subset.uns['rank_genes_groups']
            groups = results['names'].dtype.names #get names of groupby (condition)
            keys = list(results.keys()) #params, names (genes), scores, p_val, ...
            keys = keys[1:] #remove params
            #results for DE for each condition
            results_readable = {group : pd.DataFrame({key:results[key][group] for key in keys}) for group in groups}
            for group in groups:
                results_readable[group].to_csv("{}/{}/wilcoxon_{}_vs_rest_velocity.tsv".format(tsv_file_path,clust.replace("/","."),group), sep="\t")
            
            #wilcoxon rank-sum test on total RNA
            scanpy.tl.rank_genes_groups(adata_subset, groupby = condition, method = "wilcoxon")
            RNA_results = adata_subset.uns['rank_genes_groups']
            RNA_results_readable = {group : pd.DataFrame({key:RNA_results[key][group] for key in keys}) for group in groups}
            
            # merge velocity and total RNA tests
            for group in groups:
                RNA_results_readable[group] = pd.merge(results_readable[group],RNA_results_readable[group],on = ['names'])
                RNA_results_readable[group].columns = [k.replace("x","Vel").replace("y","RNA") for k in list(RNA_results_readable[group].columns)]
                RNA_results_readable[group].to_csv("{}/{}/wilcoxon_{}_RNA_velocity_merged.tsv".format(tsv_file_path,clust.replace("/","."),group), sep="\t")
            # all genes 
            df3.to_csv("{}/{}/anova_mean_velocity_genes.tsv".format(tsv_file_path,clust.replace("/",".")), sep="\t")
            
    # loop through clusters and plot the mean of the top n genes 
    #path_ext = "violin_top_{}_DE_genes".format(n)
    os.chdir(violin_files)
    for i,clust in enumerate(clusters_of_interest):
        DE_markers = pd.read_csv("{}/DE.{}_FPD.{}_Healthy.Markers.txt".format(markers_dir,clust,clust), sep = "\t", index_col = 0)
        
        DE_markers["genes"] = DE_markers.index
        
        #only significant
        significance = 0.05
        markers = DE_markers[DE_markers['p_val_adj']<significance]
        valid_genes = [gene for gene in adata.var_names]
        markers = markers[markers['genes'].isin(valid_genes)]
        #only markers that have velocities
        valid_genes = [gene for gene in list(markers['genes']) if any(~np.isnan(adata.layers["velocity"][:,adata.var_names.get_loc(gene)]))]
        markers = markers[markers['genes'].isin(valid_genes)]
        Numb_genes = [5,30,50]
        show = False
        for n in Numb_genes:
             #get genes that are up regulated in FPD
            submark = markers.nlargest(n,"avg_logFC")
            #get genes that are down regulated in FPD
            submark_s = markers.nsmallest(n,"avg_logFC")
            
            
            #plot median 
            #plot up regulated
            plot_gene_velocity(adata,gene = list(submark.index), groupby = condition, cluster = "cluster", clust = clust,layer = "velocity",median = True,stat = True,show = show,dir_base = "FPD_up",hue_order = Order_plot, palette = color_dict)
            #plot down regulated
            plot_gene_velocity(adata,gene = list(submark_s.index),groupby = condition, cluster = "cluster", clust = clust,layer = "velocity",median = True,stat = True,show = show,dir_base = "FPD_down",hue_order = Order_plot, palette = color_dict)
            #plot means
            plot_gene_velocity(adata,gene = list(submark.index), groupby = condition, cluster = "cluster", clust = clust,layer = "velocity",median = False,stat = True,show = show,dir_base = "FPD_up",hue_order = Order_plot, palette = color_dict)
            plot_gene_velocity(adata,gene = list(submark_s.index),groupby = condition, cluster = "cluster", clust = clust,layer = "velocity",median = False,stat = True,show = show,dir_base = "FPD_down",hue_order = Order_plot, palette = color_dict)
            
            velocity_by_sample(adata,genes = list(submark["genes"]), groupby = condition, cluster = "cluster", clust = clust,layer = "velocity", show = False,dir_base='FPD_up_by_gene', stat = True,hue_order = Order_plot, palette = color_dict)
            velocity_by_sample(adata,genes = list(submark_s["genes"]), groupby = condition, cluster = "cluster", clust = clust,layer = "velocity", show = False,dir_base='FPD_down_by_gene', stat = True,hue_order = Order_plot, palette = color_dict)

            plot_gene_velocity(adata,gene = list(submark['genes']), clust = clust, groupby = "samples",Order_plot = sample_names, cluster = "cluster", layer = "velocity",dir_base = 'samples_FPD_up')
            plot_gene_velocity(adata,gene = list(submark_s['genes']), clust = clust, groupby = "samples",Order_plot = sample_names, cluster = "cluster", layer = "velocity",dir_base = 'samples_FPD_down')
        
    #change output dir
    dir_violin = os.path.join(violin_files,"cross_referenced_DE")
    os.makedirs(dir_violin,exist_ok = True)
    os.chdir(dir_violin)
    for i,clust in enumerate(clusters_of_interest):    
        DE_markers = pd.read_csv("{}/DE.{}_FPD.{}_Healthy.Markers.txt".format(markers_dir,clust,clust), sep = "\t", index_col = 0)
        
        DE_markers["genes"] = DE_markers.index
        valid_genes = [gene for gene in adata.var_names]
        markers = markers[markers['genes'].isin(valid_genes)]
        #only significant
        significance = 0.05
        markers = DE_markers[DE_markers['p_val_adj']<significance]
        
        Rank_tsv = pd.read_csv("{}/{}/top_velocity_{}_genes.tsv".format(tsv_file_path,clust.replace("/","."),clust.replace("/",".")), sep = "\t", index_col = 0)
        markers = DE_markers.merge(Rank_tsv,left_index = True, right_on = 'genes')
        #only markers that have velocities
        valid_genes = [gene for gene in list(markers['genes']) if any(~np.isnan(adata.layers["velocity"][:,adata.var_names.get_loc(gene)]))]
        markers = markers[markers['genes'].isin(valid_genes)]
        #add difference between FPD and HD 
        markers['delta'] = markers['{}'.format(Order_plot[1])] - markers['{}'.format(Order_plot[0])]
        #add agreement if difference matches the sign of the logFC (upregulated in FPD will match a positive sign)
        markers['agreement'] = np.sign(markers['avg_logFC']) == np.sign(markers['delta'])
        show = False
        for n in Numb_genes:
            #get genes that are up regulated in FPD
            submark = markers.nlargest(n,"avg_logFC")
            #plot means
            if submark.shape[0] > 0:
                plot_gene_velocity(adata,gene = list(submark["genes"]), groupby = condition, cluster = "cluster", clust = clust,layer = "velocity",median = False,stat = True,show = show,dir_base = "FPD_up",hue_order = Order_plot, palette = color_dict)
                velocity_by_sample(adata,genes = list(submark["genes"]), groupby = condition, cluster = "cluster", clust = clust,layer = "velocity", show = False,dir_base = "FPD_up_by_gene", stat = True,hue_order = Order_plot, palette = color_dict)
            #get genes that are down regulated in FPD
            submark = markers.nsmallest(n,"avg_logFC")
            if submark.shape[0] > 0:
                plot_gene_velocity(adata,gene = list(submark["genes"]),groupby = condition, cluster = "cluster", clust = clust,layer = "velocity",median = False,stat = True,show = show,dir_base = "FPD_down",hue_order = Order_plot, palette = color_dict)
                velocity_by_sample(adata,genes = list(submark["genes"]), groupby = condition, cluster = "cluster", clust = clust,layer = "velocity", show = False,dir_base = "FPD_down_by_gene", stat = True,hue_order = Order_plot, palette = color_dict)
            #get genes that are up regulated in FPD and down in velocity
            submark = markers[markers['agreement'] == False].nlargest(n,"avg_logFC")
            if submark.shape[0] > 0:
                plot_gene_velocity(adata,gene = list(submark["genes"]), groupby = condition, cluster = "cluster", clust = clust,layer = "velocity",median = False,stat = True,show = show,dir_base = "FPD_up_velo_down",hue_order = Order_plot, palette = color_dict)
                velocity_by_sample(adata,genes = list(submark["genes"]), groupby = condition, cluster = "cluster", clust = clust,layer = "velocity", show = False,dir_base = "FPD_up_velo_down_by_gene", stat = True,hue_order = Order_plot, palette = color_dict)
            #get genes that are both up regulated in FPD and up in velocity
            submark = markers[markers['agreement'] == True].nlargest(n,"avg_logFC")
            if submark.shape[0] > 0:
                plot_gene_velocity(adata,gene = list(submark["genes"]), groupby = condition, cluster = "cluster", clust = clust,layer = "velocity",median = False,stat = True,show = show,dir_base = "FPD_up_velo_up",hue_order = Order_plot, palette = color_dict)
                velocity_by_sample(adata,genes = list(submark["genes"]), groupby = condition, cluster = "cluster", clust = clust,layer = "velocity", show = False,dir_base = "FPD_up_velo_up_by_gene", stat = True,hue_order = Order_plot, palette = color_dict)

            #get genes that are down regulated in FPD and up in velocity
            submark = markers[markers['agreement'] == False].nsmallest(n,"avg_logFC")
            if submark.shape[0] > 0:
                plot_gene_velocity(adata,gene = list(submark["genes"]), groupby = condition, cluster = "cluster", clust = clust,layer = "velocity",median = False,stat = True,show = show,dir_base = "FPD_down_velo_up",hue_order = Order_plot, palette = color_dict)
                velocity_by_sample(adata,genes = list(submark["genes"]), groupby = condition, cluster = "cluster", clust = clust,layer = "velocity", show = False,dir_base = "FPD_down_velo_up_by_gene", stat = True,hue_order = Order_plot, palette = color_dict)
            #get genes that are both up regulated in FPD and up in velocity
            submark = markers[markers['agreement'] == True].nsmallest(n,"avg_logFC")
            if submark.shape[0] > 0:
                plot_gene_velocity(adata,gene = list(submark["genes"]), groupby = condition, cluster = "cluster", clust = clust,layer = "velocity",median = False,stat = True,show = show,dir_base = "FPD_down_velo_down",hue_order = Order_plot, palette = color_dict)
                velocity_by_sample(adata,genes = list(submark["genes"]), groupby = condition, cluster = "cluster", clust = clust,layer = "velocity", show = False,dir_base = "FPD_down_velo_down_by_gene", stat = True,hue_order = Order_plot, palette = color_dict)
        
if __name__ == '__main__':
    main()