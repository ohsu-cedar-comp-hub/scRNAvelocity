rule make_html:
    input:
        "results/{wave}/scvelo_obs.tsv",
        "results/{wave}/cluster_colors.tsv",
        "results/{wave}/scvelo_object_analysis/pseudotime_files/figures/pseudotime_dpt_r.pdf",
        "results/{wave}/scvelo_object_analysis/pseudotime_files/figures/bar_pseudotime_dpt_lineplot.pdf"
    output: 
        html="results/{wave}/scvelo_analysis.html"
    params:
        script="scripts/Analysis.Rmd",
        seurat=config["seuratObj"],
        seurat_status=config["seurat_status"],
        seurat_cluster=config["seurat_cluster"],
        genes=genes_of_interest,
        wave = lambda wc:"{wave}".format(wave = wc.wave),
        out_dir=lambda wc:"results/{wave}/scvelo_object_analysis".format(wave = wc.wave),
        color_names = config["order_plot"],
        color_hex = config["color_hex"],
        clusters_of_interest = config["clusters_of_interest"]
    conda:
        "../envs/seurat.yaml"
    shell: 
        """
        Rscript -e 'rmarkdown::render(\"./{params.script}\", output_file = \"../{output.html}\", params=list(inputobs = \"../{input[0]}\", out_dir = \"../{params.out_dir}\", seurat=\"{params.seurat}\",contrast = \"{params.seurat_status}\",cluster = \"{params.seurat_cluster}\",genes=\"{params.genes}\",wave=\"{params.wave}\",color_tsv=\"../{input[2]}\",clusters_of_interest=\"{params.clusters_of_interest}\",color_hex=\"{params.color_hex}\",color_names=\"{params.color_names}\"))'
        """

rule seurat_downsample_DE:
    input:
        seurat_file=config["seuratObj"]
    output:
        output_file="results/DE_dir/Cluster_ds_markers.tsv"
    params:
        sample_set = config["seurat_downsample"]
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_downsample_DE.R"

rule veloin_plots:
    input:
        in_object="results/{wave}/cellrank_object.h5ad",
        markers_file="results/DE_dir/Cluster_markers_all.tsv"
    output:
        output_file="results/{wave}/scvelo_object_analysis/veloin_plots.txt"
    params:
        seurat_status = config["seurat_status"],
        color_hex = config["color_hex"],
        order_plot = config["order_plot"],
        cluster_hex = config["cluster_hex"],
        order_cluster = config["order_cluster"],
        markers_dir = config['markers_dir'],
        clusters_of_interest = config["clusters_of_interest"],
        genes = genes_of_interest,
        out_dir = lambda wc: "results/{wave}/scvelo_object_analysis".format(wave = wc.wave)
    conda:
        "../envs/scvelo.yaml"
    script:
        "../scripts/veloinplots.py"

rule combine_DE:
    input:
        ["results/DE_dir/{}_markers.tsv".format(cluster) for cluster in clusters]
    output:
       output_file = "results/DE_dir/Cluster_markers_all.tsv"
    params:
        ext = "\t"
    conda:
        "../envs/scvelo.yaml"
    script:
        "../scripts/concat_files.py"

rule seurat_DE:
    input:
        seurat_file=config["seuratObj"]
    output:
        output_file="results/DE_dir/{cluster}_markers.tsv"
    params:
        cluster = lambda wc: "{}".format(wc.cluster),
        output_dir = "results/DE_dir"
    conda:
        "../envs/seurat.yaml"
    script:
        "../scripts/seurat_DE.R"

        
rule pseudotime_r:
    input:
        input_tsv = "results/{wave}/scvelo_obs.tsv",
        color_tsv = "results/{wave}/cluster_colors.tsv" #use cluster colors made from 'pseudotime_py' rule
    output:
        output_file = "results/{wave}/scvelo_object_analysis/pseudotime_files/figures/pseudotime_dpt_r.pdf"
    params:
        clusters_of_interest = config["clusters_of_interest"],
        seurat_status = config["seurat_status"],
        color_names = config["order_plot"],
        color_hex = config["color_hex"]
    conda:
        "../envs/r_note.yaml"
    log:
        notebook="results/{wave}/R_line_processed.ipynb"
    notebook:
        "../notebooks/R_line_plot.r.ipynb"


rule pseudotime_py:
    input:
        velocity_adata = "results/{wave}/scvelo_analysis.h5ad"
    output:
        output_file = "results/{wave}/scvelo_object_analysis/pseudotime_files/figures/bar_pseudotime_dpt_lineplot.pdf",
        cluster_colors = "results/{wave}/cluster_colors.tsv" #make cluster colors not so close to one another
    params:
        clusters_of_interest = config["clusters_of_interest"],
        seurat_status = config["seurat_status"]
        
    conda:
        "../envs/scvelo_note.yaml"
    notebook:
        "../notebooks/Velocity_line_plot_pseudotime.py.ipynb"

rule analyze_scvelo:
    input:
        in_object="results/{wave}/cellrank_object.h5ad"
    output:
        out_object="results/{wave}/scvelo_obs.tsv",
        adata_out ="results/{wave}/scvelo_analysis.h5ad"
    params:
        seurat_status = config["seurat_status"],
        color_hex = config["color_hex"],
        order_plot = config["order_plot"],
        cluster_hex = config["cluster_hex"],
        order_cluster = config["order_cluster"],
        markers_dir = config['markers_dir'],
        clusters_of_interest = config["clusters_of_interest"],
        genes = genes_of_interest,
        out_dir = lambda wc: "results/{wave}/scvelo_object_analysis".format(wave = wc.wave)
    conda:
        "../envs/scvelo.yaml"
    log:
        "logs/velocity/{wave}/analyze_scvelo.log"
    script:
        "../scripts/Analyze_Cluster_Condition.py"    

rule cellrank:
    input:
        velocity_adata = "results/{wave}/scvelo_object.h5ad"
    output:
        output_file="results/{wave}/cellrank_object.h5ad"
    params:
        cell_environment = "{}/bin/activate".format(config['cr_environment']),
        script_run = "scripts/cellrank_script.py"
    shell:
        """
        source {params.cell_environment};
        python3 {params.script_run} {input.velocity_adata} {output.output_file};
        """
      
rule scvelo_batch:
    input:
        velocity_loom = "results/{wave}/sorted_merged_filtered.loom",
        seurat_loom = "results/seurat.loom"
    output: 
        out_object="results/{wave}/scvelo_object_batch.h5ad"
    params:
        seurat_cluster=config["seurat_cluster"],
        genes=genes_of_interest
    conda:
        "../envs/scvelo.yaml"
    script:
        "../scripts/scvelo_script.py"



rule scvelo:
    input:
        velocity_loom = ancient("results/{wave}/sorted_merged_filtered.loom"),
        seurat_loom = "results/seurat.loom"
    output: 
        out_object="results/{wave}/scvelo_object.h5ad"
    params:
        n_cores = config["n_cores"]
    conda:
        "../envs/scvelo.yaml"
    script:
        "../scripts/scvelo_script.py"

rule scvelo_ind:
    input:
        subset_CB=ancient(lambda wc:"{sample}/velocyto/{sample_name}.loom".format(sample = PATH_by_base[wc.sample_name], sample_name = wc.sample_name)),
        seurat_loom="results/seurat.loom"
    output:
        out_object="results/ind/{sample_name}/ind_scvelo_object.h5ad"
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
        seurat_loom = "results/seurat.loom"
    output:
        out_file="results/{wave}/sorted_merged_filtered.loom"
    params:
        indices = lambda wc: Lookup_index[wc.wave],
        seurat_cluster=config["seurat_cluster"],
        seurat_batch = config["seurat_batch"],
        seurat_status = config["seurat_status"],
        seurat_cb_correction = seurat_cb_correction
    conda:
        "../envs/velocity.yaml"
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
        input_list = ancient(lambda wc: expand("{locs}/velocyto/{base}.loom",zip, locs = PATH_by_wave[wc.wave], base = [os.path.basename(x) for x in PATH_by_wave[wc.wave]])),
        scvelo_ind = ancient(lambda wc: expand("results/ind/{base}/ind_scvelo_object.h5ad", base = [os.path.basename(x) for x in PATH_by_wave[wc.wave]]))
    output: "results/{wave}/looms/sorted_merged.loom"
    conda:
        "../envs/velocity.yaml"
    script:
        "../scripts/merge_looms.py"



rule RNAvelocity:
    input:
        ancient("{sample}/outs/cellsorted_possorted_genome_bam.bam")
    output:
        "{sample}/velocyto/{base}.loom"
    params:
        name = lambda wildcards: "{sample}".format(sample = wildcards.sample),
        n_cores=config["n_cores"],
        gtf = config["GTF_ref"],
        repeat_mask = config["repeat_mask"],
        mem=config["sam_mem"]
    conda:
        "../envs/velocity.yaml"
    shell:
        """
        velocyto run10x -m {params.repeat_mask} -@ {params.n_cores} {params.name} {params.gtf}
        """
          
rule samtools_sort:
    input:
        ancient("{sample}/outs/possorted_genome_bam.bam")
    output:
        "{sample}/outs/cellsorted_possorted_genome_bam.bam"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools sort -t CB -O BAM -o {output} {input}
        """
             
