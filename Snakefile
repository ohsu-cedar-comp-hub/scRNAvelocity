__author__ = "Joey Estabrook"
__email__ = "estabroj@ohsu.edu"
__license__ = "MIT"

"""scRNA velocity plot pipeline"""

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


with open('input/paths.txt') as file_in:
	 for i,line in enumerate(file_in):
			if i == 0:
				continue #skip header line
			else:
				location = line.split("\t")[0]
				uniq	=	line.split("\t")[1]
				PATHS.append(location)
				UNIQ.append(uniq)
				
	  
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

for wave in PATHS:
	message("10x files in " + wave + " will be processed")

rule all:
	input:
		expand(["{wave}/outs/cellsorted_possorted_genome_bam.bam"],wave=PATHS),
		["{wave}/velocyto/{base}.loom".format(wave = wave, base = os.path.basename(wave)) for wave in PATHS],
		"results/looms/sorted_merged.loom",
		expand(["results/{seurat}/velocity_plot.png"],seurat=SEURAT)
		
		
		
include: "rules/velocyto.smk"
