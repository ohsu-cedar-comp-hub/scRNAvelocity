#!/usr/bin/bash
#SBATCH --time 35:00:00
#SBATCH --partition exacloud
#SBATCH --job-name workflow_submission
#SBATCH --output=logs/workflow_submission_%j.log
#SBATCH --mail-type=END,FAIL


snakemake -j 100 --use-conda --rerun-incomplete  --printshellcmds --cluster-config cluster.json --cluster "sbatch -p {cluster.partition} -N {cluster.N}  -t {cluster.t} -o {cluster.o} -e {cluster.e} -J {cluster.J} -c {cluster.c} --mem {cluster.mem}" -s Snakefile --latency-wait 120 


#test build
#snakemake --use-conda --cores=4 -np
#
#snakemake --use-conda --cores=4 --printshellcmds --reason -np