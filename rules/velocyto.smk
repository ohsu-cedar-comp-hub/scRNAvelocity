rule make_html:
	input:
		in_object="results/{seurat}/Analyze/scvelo_obs.tsv"
	output: 
		html="results/{seurat}/Analyze/scvelo_analysis.html"
	params:
		base_directory= lambda wildcards: "./results/{seurat}".format(seurat = wildcards.seurat),
		script="scripts/Analysis.Rmd",
		seurat=lambda wildcards: "{}".format(wildcards.seurat),
		genes=genes_of_interest
	conda:
		"../envs/seurat.yaml"
	shell: 
		"""
		Rscript -e 'rmarkdown::render(\"./{params.script}\", output_file = \"../{output.html}\", params=list(inputobs = \"../{input.in_object}\", directory = \"../{params.base_directory}\", seurat=\"{params.seurat}\",genes=\"{params.genes}\"))'
		"""

rule analyze_scvelo:
	input:
		in_object="results/{seurat}/scvelo_object.h5ad"
	output: 
		out_object="results/{seurat}/Analyze/scvelo_obs.tsv"
	params:
		genes=genes_of_interest
	conda:
		"../envs/scvelo.yml"
	script:
		"../scripts/Analyze_Cluster_Condition.py"

rule scvelo_batch:
	input:
		velocity_loom = "results/{seurat}/sorted_merged_filtered.loom",
		seurat_loom = "results/{seurat}/{seurat}.loom"
	output: 
		out_object="results/{seurat}/scvelo_object_batch.h5ad"
	params:
		seurat_cluster=config["seurat_cluster"],
		seurat_sample = config["seurat_sample"],
		genes=genes_of_interest
	conda:
		"../envs/scvelo.yml"
	script:
		"../scripts/scvelo.py"


rule scvelo:
	input:
		velocity_loom = "results/{seurat}/sorted_merged_filtered.loom",
		seurat_loom = "results/{seurat}/{seurat}.loom"
	output: 
		out_object="results/{seurat}/scvelo_object.h5ad"
	params:
		genes=genes_of_interest
	conda:
		"../envs/scvelo.yml"
	script:
		"../scripts/scvelo.py"

rule scvelo_ind:
	input:
		subset_CB=lambda wc:"{sample}/velocyto/{sample_name}.loom".format(sample = PATH_by_base[wc.sample_name], sample_name = wc.sample_name),
		velocity_loom = lambda wc:"results/{seurat}/sorted_merged_filtered.loom".format(seurat = wc.seurat_sample)
	output:
		out_object="results/ind/{seurat_sample}/{sample_name}/scvelo_object.h5ad"
	params:
		subset_CB=lambda wc:"{}".format(wc.sample_name),
		genes=genes_of_interest
	conda:
		"../envs/scvelo.yml"
	script:
		"../scripts/scvelo_ind.py"
	
rule plot_velocity:
	input:
		velocity_loom = "results/{seurat}/sorted_merged_filtered.loom",
		seurat_loom = "results/{seurat}/{seurat}.loom"
	output: 
		out_plot="results/{seurat}/velocity_plot.png"
	params:
		seurat_cluster=config["seurat_cluster"]
	conda:
		"../envs/velocity.yml"
	script:
		"../scripts/plot_velocity.py"
		
rule correct_CB:
	input:
		velocity_loom = "results/looms/sorted_merged.loom",
		seurat_loom = "results/{seurat}/{seurat}.loom"
	output:
		out_file="results/{seurat}/sorted_merged_filtered.loom"
	params:
		seurat_cluster=config["seurat_cluster"],
		seurat_sample = config["seurat_sample"]
	conda:
		"../envs/velocity.yml"
	script:
		"../scripts/correct_CB.py"

rule seurat_to_loom:
	input:
		seurat_file="input/{seurat}.rds"
	output: 
		out_loom="results/{seurat}/{seurat}.loom"
	conda:
		"../envs/seurat.yaml"
	script:
		"../scripts/seurat2loom.R"

rule loom_merge:
	input:  expand("{wave}/velocyto/{base}.loom",zip, wave = PATHS, base = [os.path.basename(x) for x in PATHS])
	output: "results/looms/sorted_merged.loom"
	conda:
		"../envs/velocity.yml"
	script:
		"../scripts/merge_looms.py"


rule RNAvelocity:
	input:
		input_bam="{sample}/outs/cellsorted_possorted_genome_bam.bam"
	output:
		"{sample}/velocyto/{base}.loom"
	params:
		name = lambda wildcards: "{sample}".format(sample = wildcards.sample),
		n_cores=config["n_cores"],
		gtf = config["GTF_ref"],
		repeat_mask = config["repeat_mask"],
		mem=config["sam_mem"]
	conda:
		"../envs/velocity.yml"
	shell:
		"""
		velocyto run10x -m {params.repeat_mask} -@ {params.n_cores} {params.name} {params.gtf}
		"""
		  
rule samtools_sort:
	input:
		"{sample}/outs/possorted_genome_bam.bam"
	output:
		"{sample}/outs/cellsorted_possorted_genome_bam.bam"
	conda:
		"../envs/samtools.yml"
	shell:
		"""
		samtools sort -t CB -O BAM -o {output} {input}
		"""
			 
