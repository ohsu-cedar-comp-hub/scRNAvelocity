

	
rule plot_velocity:
	input:
		velocity_loom = "results/looms/sorted_merged.loom",
		seurat_loom = "results/seurat/{seurat}.loom"
	output: 
		out_plot="results/{seurat}/velocity_plot.png"
	params:
		cluster=config["seurat_cluster"]
	conda:
		"../envs/velocity.yml"
	script:
		"../scripts/plot_velocity.py"

rule seurat_to_loom:
	input:
		seurat_file="input/{seurat}.rds"
	output: 
		out_loom="results/seurat/{seurat}.loom"
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
			 
