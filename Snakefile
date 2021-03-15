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


#genes_of_interest = ["RUNX1","CD74","MIF","FOS","CCL2", "PU.1", "TLR4", "TLR2","CD44", "SRC", "PI3K","AKT", "TEN", "JAB1", "CXCL8", "MAPK", "ERK", "SPI1", "STAT3", "STAT5", "NFKb", "CXCR1/2","CXCR1","CXCR2",  "JUN", "GATA1"]
seurat_cb_correction = '_'
genes_of_interest = ["MYC","CEBPA","CEBPE","IRF8","GFI1","Prom1","Itgam","E2F3","Ly6g","E2F3","CDK6","CCND2","KIT","Rpl37","Eef2"]
genes_of_interest = [x[0].upper() + x[1:].lower() for x in genes_of_interest] #mouseify gene names

PATH_by_wave = {}
Index_by_wave = {}
PATH_by_base = {}
Base_by_PATH = {} 
Lookup_index = {}
Index_by_base = {}
WAVES,=glob_wildcards('input/{wave}_locations.tsv')

for wave in WAVES:
	with open('input/{wave}_locations.tsv'.format(wave = wave)) as file_in:
		Locations = []
		Indexes = []
		lookup = {}
		for i,line in enumerate(file_in):
			if i == 0:
				continue #skip header line
			else:
				location = line.split("\t")[0]
				uniq	=	line.split("\t")[1].rstrip()
				PATHS.append(location)
				UNIQ.append(uniq)
				Locations.append(location)
				Indexes.append(uniq)
				PATH_by_base[os.path.basename(location)] = location #assuming basenames are unique
				Base_by_PATH[location] = os.path.basename(location)
				lookup[os.path.basename(location)] = uniq
				Index_by_base[os.path.basename(location)] = uniq #assuming basenames are unique
		PATH_by_wave[wave] = Locations
		Index_by_wave[wave] = Indexes
		Lookup_index[wave] = lookup


PATHS = list(set(PATHS)) #get unique paths to 10x outputs
  

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

for wave in WAVES:
	message("10x files in " + wave + " will be processed")

#sample_name = [os.path.basename(locations) for wave in WAVES for locations in PATH_by_wave[wave]]
#seurat = os.path.splitext(os.path.basename(config["seuratObj"]))[0]	
rule all:
	input:
		#expand("{sample}/velocyto/{base}.loom",zip, sample = PATHS, base = [os.path.basename(locations) for locations in PATHS]),
		#expand("results/{wave}/looms/sorted_merged.loom",wave = WAVES),
		#expand(["results/{wave}/scvelo_object_batch.h5ad"],wave = WAVES),
		#expand(["results/{wave}/scvelo_object.h5ad"],wave = WAVES),
		expand("results/ind/{sample_name}/ind_scvelo_object.h5ad", sample_name = [os.path.basename(locations) for locations in PATHS]),
		expand("results/{wave}/Analyze/scvelo_obs.tsv",wave = WAVES),
		expand("results/{wave}/Analyze/scvelo_analysis.html", wave = WAVES)

include: "rules/velocyto.smk"
