#!/usr/bin/env python
"""
Author: Trevor Enright

Input a tab separated meta file with a header row.  first column is location to 10x output "../outs/filtered_feature_bc_matrix". 2nd column is a unique name for that sample without spaces

output a merged bam file with all 

"""
import os
import sys
import gzip
import glob
import subprocess
import time

######################
#  GLOBAL VARIABLES  #
######################

input_file = snakemake.input.meta
output_file = snakemake.output.file
output_dir = snakemake.params.outdir
cores = snakemake.params.cores


jobname = "velocyto_run"
output = "log_slurm"
mem = "20G"
time_lim = "6:00:00"
gtf_ref = "/home/groups/precepts/enright/experiments/single_cell/GTF_files/refdata-cellranger-hg19-3.0.0/genes/genes.gtf"
repeat_mask = "/home/groups/precepts/enright/experiments/single_cell/GTF_files/hg38_repeatmask.gtf"
cores = "4"
sam_mem = "{}".format(int(int(mem.split("G")[0])/int(cores)*1000))
charcater_encoding = 'utf-8'


def mkdir_path(dir):
	'''make a directory (dir) if it doesn't exist already'''
	if not os.path.exists(dir):
		os.mkdir(dir)


mkdir_path(output_dir)

def submit_with_meta_file():
	"""
	submits multiple jobs to slurm running velocity on bam files
	returns the job ids
	"""
	jobid = []
	with open(input_file,"r") as fh:
		for i,line in enumerate(fh):
			if i == 0:
				continue #skip header line
			else:
				print("processing {} bam file".format(i))
				location=line.split("\t")[0]
				uniq=line.split("\t")[1]
				bam = "{}/outs/possorted_genome_bam.bam".format(location.split("/outs/")[0])
				barcodes = glob.glob("{}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz".format(location.split("/outs/")[0])) #zipped
				if not barcodes:
					barcodes = "{}/outs/filtered_feature_bc_matrix/barcodes.tsv".format(location.split("/outs/")[0]) #not zipped
				#else:
				#	with gzip.open(barcodes,"rb") as gbcode:
				
				#write sbatch file
				
				job_file = "{}/run_{}.sh".format(output_dir,uniq)
				sbatch_skeleton = "#!/bin/bash" + "\n" + "#SBATCH --job-name={}_{}".format(jobname,uniq) + "\n" + "#SBATCH --output={}_{}".format()
				with open(job_file,"w") as sf:
					sf.writelines("#!/bin/bash\n")
					sf.writelines("#SBATCH --job-name={}_{}\n".format(jobname,uniq))
					sf.writelines("#SBATCH --output={}.out\n".format(uniq))
					sf.writelines("#SBATCH --time={}\n".format(time_lim))
					sf.writelines("#SBATCH --mem={}\n".format(mem))
					sf.writelines("#SBATCH --cpus-per-task={}\n".format(cores))
					sf.writelines("#SBATCH --ntasks-per-node=1\n")
					sf.writelines("velocyto run -e {uniq} -b {barcode} -o {out_dir} -m {repeat_mask} -@ {cores} --samtools-memory {sam} {input_bam} {gtf}\n".format(gtf=gtf_ref,out_dir = output_dir,input_bam=bam,barcode = barcodes[0],uniq=uniq,repeat_mask = repeat_mask,sam=sam_mem))
					#sf.writelines("echo 'hello world'; sleep 120")
					
				p = subprocess.run(["sbatch", "{}".format(job_file)], capture_output = True)
				job_id = p.stdout.split()[-1].decode(charcater_encoding)
				jobid.append(job_id)
				
	return(jobid)
	
def check_for_completion(job_id_list, boolean = True):
	"""check status of job ids that was submitted to slurm
	returns True or False if all job id's have "COMPLETED" state and boolean is set to True
	returns the dictionary of job id's and values of State if boolean is not True
	"""
	check_status = {}
	for job in job_id_list:
		p = subprocess.run(["sacct","-j", "{}".format(job_file),"--format","State"], capture_output = True)
		p_string = str(p.stdout)[4] #State ie: COMPLETED, FAILED, RUNNING, PENDING
		check_status[job] = p_string
	if (boolean):
		return all(elem == "COMPLETED" for elem in check_status.values())
	else:
		return check_status
		
def write_output_file(jobid_dictionary,filename):
	with open(filename,"w") as fh:
		fh.write(json.dumps(jobid_dictionary))
		
def main():
	job_list = submit_with_meta_file()
	t = time.time()
	time_allowed = (int(time_lim.split(":")[0])+1)*60*60 #time allowed to wait for jobs is the hours of time_lim + 1 hour
	some_time = 60 #check every minute
	while True:
		# Break if this takes more than some_limit
		if time.time() - t > time_allowed:
			#check to see if there should be more time
			job_dict = check_for_completion(job_id_list,False)
			if [k
			if [k for k in job_dict.values() if k == "PENDING" or k == "RUNNING"]:
				time_allowed += 360 # add an hour if pending or running
				continue
			else:
				base = os.path.dirname(output_file)
				file = "failed_{}".format(os.path.basename(output_file))
				full_file = os.path.join(base,file)
				write_output_file(job_dict,full_file) #failed
				return 1
		# Check if the jobs are done. This could be done by
		# grep'ing squeue for your username and some tags
		# that you name your jobs
		if(check_for_completion(job_list,True)):
			write_output_file(check_for_completion(job_id_list,False),output_file) #success
			return 0
		# Sleep for a while depending on the estimated completion time of the jobs
		time.sleep(some_time)

main()