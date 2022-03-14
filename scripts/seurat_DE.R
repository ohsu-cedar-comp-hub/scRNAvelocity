#!/usr/bin/env Rscript

#load or install packages
packages_to_work <- c("Seurat","dplyr","ggplot2","tidyr")

is_installed <- function(mypkg){ is.element(mypkg, installed.packages()[,1])}
load_or_install<-function(package_names){
    r = getOption("repos")
    r["CRAN"] = "http://cran.us.r-project.org"
    options(repos = r)
	#https://www.r-bloggers.com/2012/05/loading-andor-installing-packages-programmatically/
	#quick install or load packages
  for(package_name in package_names)
  {
    if(!is_installed(package_name))
    {
       message(sprintf("installing package: %s",package_name))
       install.packages(package_name,repos = "http://cran.us.r-project.org")
    }
    message(sprintf("loading package: %s",package_name))
    library(package_name,character.only=TRUE,quietly=TRUE,verbose=FALSE)
  }
}

load_or_install(packages_to_work)


output_file <- file.path(getwd(),snakemake@output[["output_file"]])
seurat_file <- snakemake@input[["seurat_file"]]
output_dir <- file.path(getwd(),snakemake@params[["output_dir"]])

seuratObj <- readRDS(seurat_file)

cluster_call = snakemake@params[['cluster']]#gsub("_markers.tsv","",basename(output_file))

dir.create(output_dir, showWarnings = FALSE)

Idents(seuratObj) <- "predicted"
DefaultAssay(seuratObj) <- 'RNA'
type = "RNA"
#cluster_markers <- FindAllMarkers(seuratObj, assay = type,only.pos = F, logfc.threshold=0.0, min.pct = 0)

#possibleIdents <- unique(seuratObj$predicted) #cell types
message(sprintf('Calculating %s cluster',cluster_call))
#cluster_file <- file.path(output_dir,sprintf('%s_markers.tsv',cluster_call))
cluster_markers <- FindMarkers(seuratObj, assay = type, ident.1 = cluster_call, logfc.threshold=0.1, min.pct = 0.1 ) #Seurat Default
cluster_markers['cluster'] <- rep(cluster_call,dim(cluster_markers)[1])
cluster_markers['gene'] <- rownames(cluster_markers)
write.table(cluster_markers,output_file,sep='\t', quote= F, row.names = F, col.names = T)
