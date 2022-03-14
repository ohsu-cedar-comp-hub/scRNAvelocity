#!/usr/bin/env python3

import pandas as pd
import os

file_list = snakemake.input
output = snakemake.output.output_file
file_extention = snakemake.params.ext


#read files
files = [pd.read_csv(x, sep = file_extention) for x in file_list]
#concat all files
res = pd.concat(files, axis =0, ignore_index=True)

#write to output

res.to_csv(output,sep = file_extention)