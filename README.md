# RNAVelocity

Runs velocyto on multiple 10x outputs from cellranger to project trajectory arrows onto an UMAP from a partially processsed merged/integrated seurat object.

[Velocyto Reference](https://velocyto.org/velocyto.py/tutorial/cli.html#run10x-run-on-10x-chromium-samples)

[Velocyto preprint: RNA velocity in single cells](https://www.biorxiv.org/content/10.1101/206052v1)

## Setup


Add 10x paths into the "[name without underscores]_locations.tsv" file in the input directory.

SAMPLE:
```
location	name
/PATH/TO/SAMPLES/10x_1 1
/PATH/TO/SAMPLES/10x_2 2
```

The "outs" directory from 10x cellranger would be a sub-directory of the "/PATH/TO/SAMPLES/10x_1".  The name column references the postfix of the cell barcode used in the seurat object.


Add the paths in the *omic_config.yaml* for the GTF files.
As well as the path to the seurat object with an *.rds* extension.


Repeat masker downloaded from [UCSC](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=611454127_NtvlaW6xBSIRYJEBI0iRDEWisITa&clade=mammal&org=&db=hg38&hgta_group=allTracks&hgta_track=rmsk&hgta_table=rmsk&hgta_regionType=genome&position=&hgta_outputType=gff&hgta_outFileName=hg19_repeatmask.gtf). 
Repat masker is optional to velocyto but it is required here.


Reference gtf annotation files from 10x.
[Human 10x gtf](http://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-hg19-3.0.0.tar.gz)

The links are to hg19.

Example omics file.

*omic_config.yaml*
```
GTF_ref:
    /path/to/10x.gtf
repeat_mask:
    /path/to/repeat_masker.gtf
seuratObj:
    /path/to/seurat_object.rds
seurat_cluster:
    RNA_snn_res.0.8 # seurat meta data for naming clusters
seurat_status:
    Condition #seurat meta data for contrast to be analyzed
seurat_batch:
    orig.ident #seurat meta data for batch correction to be applied
contrast:
    FPD  #value from seurat_status to be used as contrast
sam_mem:
    50000 #memory allocation for sam tools
n_cores:
    8 #number of cores to provide certain algorithms
cr_environment:
    /path/to/environment/cellrank #path to cell rank environment

order_plot: #names of values in 'seurat_status' to be used for plotting
    - HD
    - FPD
color_hex:  #colors to match with 'order_plot' entries
    - '#67CAD4' 
    - '#B0207E'
order_cluster: #names of values in 'seurat_cluster' entry
    - HSC
    - Progenitor
    - MKP
    - GMP
    - Pro-Mono
    - Mono
    - pDC
    - cDC
    - Early Eryth
    - Late Eryth
    - CLP
    - Pro-B
    - B
    - Plasma
    - CTL
    - NK
    - T
    - Stroma
cluster_hex: #colors to match with 'order_cluster' entrie
    - '#F8766D'
    - '#E88526'
    - '#D39200'
    - '#B79F00'
    - '#93AA00'
    - '#5EB300'
    - '#00BA38'
    - '#00BF74'
    - '#00C19F'
    - '#00BFC4'
    - '#00B9E3'
    - '#00ADFA'
    - '#619CFF'
    - '#AE87FF'
    - '#DB72FB'
    - '#F564E3'
    - '#FF61C3'
    - '#FF699C'
clusters_of_interest: #clusters of interest. enteries are a subset of 'order_cluster'
    - HSC
    - Progenitor
    - MKP
    - GMP
    - Mono    
```


## Marker Directory

Has output `tsv` files from Seurats `FindMarkers()` function. Format of the file names are as follows: 'DE.{}_FPD.{}_Healthy.Markers.txt' where the {} are place holders for the cluster/cell type from the `seurat_cluster` field.

## Cell Rank

Cell rank is installed as a python virtual environment separately and then referenced the path to the environment in the rule. A Mesasage passing interface (MPI) are required to leverage parallel computing.  Iie: `module load  /path/to/software/modules/openmpi/3.1.6-openib`.


### cellrank installation: 

```
cd envs
python3 -m venv cellrank
cd cellrank
source bin/activate
pip install --upgrade pip
pip install -r ../cellrank_requirements.txt
```

## Test Build

Activate snakemake environment and run the following:

`$snakemake --use-conda --cores=4 --printshellcmds --reason -np`

## Run on cluster

`$sbatch submit_snakemake.sh`

## Common Errors

Velocyto uses the basename of the path to the 10x location directory and this can cause problems. The script *correct_CB.py* corrects for this by using the "paths.tsv" file as lookup table to connect  Velocyto output with seurat object cell barcodes. 

## Reference
[Velocyto Analysis](https://velocyto.org/velocyto.py/tutorial/analysis.html)

[scVelo](https://scvelo.readthedocs.io)

