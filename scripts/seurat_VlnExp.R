#!/usr/bin/env Rscript


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
    else{
        message(sprintf("loading package: %s",package_name))
        library(package_name,character.only=TRUE,quietly=TRUE,verbose=FALSE)
    }
  }
}


load_or_install(c("Seurat","DT","grid","plyr","dplyr","ggplot2","plotly","knitr","future","patchwork","gridExtra","RColorBrewer","data.table","extrafont"))


args <- commandArgs(trailingOnly = TRUE)

seurat_file <- args[1]
output_dir <- args[2]

seuratObj <- readRDS(seurat_file)

dir.create(output_dir, showWarnings = FALSE)

equalSampleSizeDownsample <- function(Object, downsampleFactor, DEFoctorsOfInterest, DEGroupBy, downsampleCellsNumber = NULL, seed = 8) {
  set.seed(seed)
  
  ## Cells To Subset
  cellBarcodes <- c()
  
  ## DE cells of interest
  DEFactors  <- Object[[DEGroupBy, drop=TRUE]]
  DECellsOfInterest <- names(grep(paste(paste0("^", DEFoctorsOfInterest, "$"), collapse = "|"), DEFactors, value = T))
  #print(DECellsOfInterest)
  ## DE cells 
  downsampleFactors <- Object[[downsampleFactor, drop=TRUE]]
  downsampleFactorsUniq <- unique(downsampleFactors)
  
  for (i in 1:length(downsampleFactorsUniq)) {
    downsampleCellsByFactor <- names(grep(paste0("^", downsampleFactorsUniq[i], "$"), downsampleFactors, value = T))
    
    ## Intersect the barcodes to actually sample
    barcodesToDownsample <- intersect(downsampleCellsByFactor, DECellsOfInterest)
    if (is.null(downsampleCellsNumber)){
        downsampleCellsNumber = length(barcodesToDownsample)
    }
    ## Skip samples that have no cells of that type 
    if (length(barcodesToDownsample) > 0) {
      if(length(barcodesToDownsample) > downsampleCellsNumber) {
        sampledBarcodes<- sample(barcodesToDownsample, size = downsampleCellsNumber)
        cellBarcodes <- c(cellBarcodes, sampledBarcodes)
      } else {
        cellBarcodes <- c(cellBarcodes, barcodesToDownsample)
      }
    }
  }
  Object <- subset(Object, cells = cellBarcodes)
  return(Object)
}


### clean up seurat file



SampleNames <- names(table(seuratObj$orig.ident))
Donor <- gsub("(FPD[0-9]+)_.*", "\\1", SampleNames)
Donor <- gsub("(HD[0-9]+)_.*", "\\1", Donor)
Donor <- gsub("Healthy_", "", Donor)
Donor <- gsub("HD_4_2", "HD4", Donor)
names(Donor) <- SampleNames
seuratObj$Donor <- plyr::revalue(seuratObj$orig.ident, Donor)

seuratObj$celltype.Condition <- paste(seuratObj$predicted, seuratObj$Condition, sep = "_")

Idents(seuratObj) <- 'predicted'
seuratObj$Condition <- factor(x = seuratObj$Condition, levels = c('HD','FPD'))
###

### list of cell types and markers to make vln plots ###
List_cell_genes <- list()
List_cell_genes[['Progenitor']] <- c("CD34","MSI2","CD38")
List_cell_genes[['HSC']] <- c("CD34","MSI2","CD38")
List_cell_genes[['MKP']] <- c("LYZ", "AZU1", "MPO")
List_cell_genes[['GMP']] <- c("LYZ", "AZU1", "MPO")
List_cell_genes[['Pro-Mono']] <- c("LYZ", "MPO", "ELANE", "FCER1G", "CD14")
List_cell_genes[['Mono']] <- c("LYZ", "MPO", "ELANE", "FCER1G", "CD14")
List_cell_genes[['Late Eryth']] <- c("HBB", "HBD")
List_cell_genes[['Early Eyth']] <-  c("HBB", "HBD")
List_cell_genes[['B']] <- c( "CD79A")
List_cell_genes[['T']] <- c("IL7R", "IL32")



# make plots and save

for (name in names(List_cell_genes)){
    markers <- (List_cell_genes[[name]])
    down_seurat <- equalSampleSizeDownsample(Object = seuratObj, 
                                        downsampleFactor = "Donor", 
                                        DEGroupBy = "celltype.Condition", 
                                        DEFoctorsOfInterest = c(sprintf("%s_HD",name), sprintf("%s_FPD",name)),downsampleCellsNumber = 150)
    for (marker in markers){
        file_name <- sprintf("%s/%s/%s_vlnplot.pdf", output_dir, name, marker)
        dir.create(dirname(file_name), showWarnings = FALSE)
        p_plot <- VlnPlot(down_seurat,marker,group.by = 'Condition',cols = c('#67CAD4','#B0207E'))
        ggsave(filename = basename(file_name), plot = p_plot, path = dirname(file_name), device = "pdf")
        
        file_name <- sprintf("%s/%s/%s_no_dots.pdf", output_dir, name, marker)
        p_plot <- VlnPlot(down_seurat,marker,group.by = 'Condition',cols = c('#67CAD4','#B0207E'),pt.size =0)
        ggsave(filename = basename(file_name), plot = p_plot, path = dirname(file_name), device = "pdf")
    }
}