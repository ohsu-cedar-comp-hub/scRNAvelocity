__author__ = "Joey Estabrook"
__email__ = "estabroj@ohsu.edu"
__license__ = "MIT"

"""scRNA velocity plot pipeline"""
#from snakemake.io import expand

import numpy as np
import datetime
import sys
import os
import pandas as pd
import json
from pathlib import Path


#samples
PATHS = []
UNIQ = []

PATH_by_base = {}

with open('input/paths.tsv') as file_in:
	 for i,line in enumerate(file_in):
			if i == 0:
				continue #skip header line
			else:
				location = line.split("\t")[0]
				uniq	=	line.split("\t")[1].rstrip()
				PATHS.append(location)
				UNIQ.append(uniq)
				PATH_by_base[os.path.basename(location)] = location
	  
SEURAT,=glob_wildcards('input/{seurat}.rds')

def get_contrast_global(wildcards):
	"""Return each contrast provided in the configuration file"""
	return config["diffexp"]["global_contrasts"][wildcards.contrast_DE]

def get_contrast_local(wildcards):
	"""Return each contrast provided in the configuration file"""
	return config["diffexp"]["local_contrasts"][wildcards.contrast_DE]
	
timestamp = ('{:%Y-%m-%d_%H:%M:%S}'.format(datetime.datetime.now()))

configfile:"omic_config.yaml"

with open('cluster.json') as json_file:
	json_dict = json.load(json_file)

rule_dirs = list(json_dict.keys())

for rule in rule_dirs:
	if not os.path.exists(os.path.join(os.getcwd(),'logs',rule)):
		log_out = os.path.join(os.getcwd(), 'logs', rule)
		os.makedirs(log_out)
		print(log_out)


def message(mes):
	sys.stderr.write("|--- " + mes + "\n")

for sample in PATHS:
	message("10x files in " + sample + " will be processed")


	
rule all:
	input:
		expand("{wave}/velocyto/{base}.loom",zip, wave = PATHS, base = [os.path.basename(x) for x in PATHS]),
		"results/looms/sorted_merged.loom",
		expand(["results/{seurat}/scvelo_object_batch.h5ad"],seurat=SEURAT),
		expand(["results/{seurat}/scvelo_object.h5ad"],seurat=SEURAT),
		expand("results/ind/{seurat_sample}/{sample_name}/scvelo_object.h5ad", seurat_sample = SEURAT, sample_name = [os.path.basename(x) for x in PATHS])
		
		
		
include: "rules/velocyto.smk"
