rule make_html:
	input:
		in_object="results/{wave}/Analyze/scvelo_obs.tsv"
	output: 
		html="results/{wave}/Analyze/scvelo_analysis.html"
	params:
		script="scripts/Analysis.Rmd",
		seurat=config["seuratObj"],
		seurat_status=config["seurat_status"],
		seurat_cluster=config["seurat_cluster"],
		genes=genes_of_interest,
		wave = lambda wc: "{}".format(wc.wave),
		out_dir="results/{wave}/Analyze"
	conda:
		"../envs/seurat.yml"
	shell: 
		"""
		Rscript -e 'rmarkdown::render(\"./{params.script}\", output_file = \"../{output.html}\", params=list(inputobs = \"../{input.in_object}\", out_dir = \"../{params.out_dir}\", seurat=\"{params.seurat}\",contrast = \"{params.seurat_status}\",cluster = \"{params.seurat_cluster}\",genes=\"{params.genes}\",wave=\"{params.wave}\"))'
		"""

rule analyze_scvelo:
	input:
		in_object="results/{wave}/scvelo_object.h5ad"
	output: 
		out_object="results/{wave}/Analyze/scvelo_obs.tsv"
	params:
		genes=genes_of_interest,
		seurat_status=config["seurat_status"],
		seurat_cluster=config["seurat_cluster"]
	conda:
		"../envs/scvelo_env.yaml"
	script:
		"../scripts/Analyze_Cluster_Condition.py"

rule scvelo_batch:
	input:
		velocity_loom = "results/{wave}/looms/sorted_merged_filtered.loom"
	output: 
		out_object="results/{wave}/scvelo_object_batch.h5ad"
	params:
		seurat_cluster=config["seurat_cluster"],
		genes=genes_of_interest
	conda:
		"../envs/scvelo.yaml"
	script:
		"../scripts/scvelo.py"


rule scvelo:
	input:
		velocity_loom = ancient("results/{wave}/looms/sorted_merged_filtered.loom")
	output: 
		out_object="results/{wave}/scvelo_object.h5ad"
	params:
		genes=genes_of_interest,
		seurat_status=config["seurat_status"],
		seurat_cluster=config["seurat_cluster"]
	conda:
		"../envs/scvelo.yaml"
	script:
		"../scripts/scvelo.py"

rule scvelo_ind:
	input:
		subset_CB=lambda wc:ancient("{sample}/velocyto/{sample_name}.loom".format(sample = PATH_by_base[wc.sample_name], sample_name = wc.sample_name)),
		seurat_loom=ancient("results/seurat.loom")
	output:
		out_object="results/ind/{sample_name}/ind_scvelo_object.h5ad",
		obser = "results/ind/{sample_name}/scvelo_obs.tsv"
	params:
		indices = lambda wc: Index_by_base[wc.sample_name],
		subset_CB=lambda wc:"{}".format(wc.sample_name),
		genes=genes_of_interest,
		seurat_cluster=config["seurat_cluster"],
		seurat_batch = config["seurat_batch"],
		seurat_status = config["seurat_status"],
		seurat_cb_correction = seurat_cb_correction
	conda:
		"../envs/scvelo.yaml"
	script:
		"../scripts/scvelo_ind.py"
	

		
rule correct_CB:
	input:
		velocity_loom = ancient("results/{wave}/looms/sorted_merged.loom"),
		seurat_loom = ancient("results/seurat.loom")
	output:
		out_file="results/{wave}/looms/sorted_merged_filtered.loom"
	params:
		indices = lambda wc: Lookup_index[wc.wave],
		seurat_cluster=config["seurat_cluster"],
		seurat_batch = config["seurat_batch"],
		seurat_status = config["seurat_status"],
		seurat_cb_correction = seurat_cb_correction
	conda:
		"../envs/velocity.yml"
	script:
		"../scripts/correct_CB.py"

rule seurat_to_loom:
	input:
		seurat_file=config["seuratObj"]
	output:
		out_loom="results/seurat.loom"
	conda:
		"../envs/seurat.yaml"
	script:
		"../scripts/seurat2loom.R"

rule loom_merge:
	input:
		input_list = lambda wc: [ancient(x) for x in expand("{locs}/velocyto/{base}.loom",zip, locs = PATH_by_wave[wc.wave], base = [os.path.basename(x) for x in PATH_by_wave[wc.wave]])],
		scvelo_ind = lambda wc: [ancient(x) for x in expand("results/ind/{base}/ind_scvelo_object.h5ad", base = [os.path.basename(x) for x in PATH_by_wave[wc.wave]])]
	output: "results/{wave}/looms/sorted_merged.loom"
	conda:
		"../envs/velocity.yml"
	script:
		"../scripts/merge_looms.py"


rule RNAvelocity:
	input:
		"{sample}/outs/cellsorted_possorted_genome_bam.bam"
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
			 
