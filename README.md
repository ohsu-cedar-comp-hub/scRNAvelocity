# RNAVelocity

Runs velocyto on multiple 10x outputs from cellranger to project trajectory arrows onto an UMAP from a partially processsed merged/integrated seurat object.

[Velocyto Reference](https://velocyto.org/velocyto.py/tutorial/cli.html#run10x-run-on-10x-chromium-samples)

[Velocyto preprint: RNA velocity in single cells](https://www.biorxiv.org/content/10.1101/206052v1)

## Setup

Put the seurat object rds file in the **input** directory with *.rds* extension.

Add 10x paths into the "paths.tsv" file in the input directory.

SAMPLE:
```
location	name
/PATH/TO/SAMPLES/10x_1 1
/PATH/TO/SAMPLES/10x_2 2
```

The "outs" directory from 10x cellranger would be a sub-directory of the "/PATH/TO/SAMPLES/10x_1".  The name column references the postfix of the cell barcode used in the seurat object.


Add the paths in the *omic_config.yaml* for the GTF files.


Repeat masker downloaded from [UCSC](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=611454127_NtvlaW6xBSIRYJEBI0iRDEWisITa&clade=mammal&org=&db=hg38&hgta_group=allTracks&hgta_track=rmsk&hgta_table=rmsk&hgta_regionType=genome&position=&hgta_outputType=gff&hgta_outFileName=hg19_repeatmask.gtf). 
Repat masker is optional to velocyto but it is required here.


Reference gtf annotation files from 10x.
[Human 10x gtf](http://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-hg19-3.0.0.tar.gz)

The links are to hg19.

*omic_config.yaml*
```
GTF_ref:
    /path/to/10x.gtf
repeat_mask:
    /path/to/repeat_masker.gtf
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

