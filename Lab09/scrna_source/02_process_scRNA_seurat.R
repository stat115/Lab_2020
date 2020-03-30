library(Seurat)
library(dplyr)
library(Matrix)
library(data.table)
library(cowplot)
library(ggplot2)

options(future.globals.maxSize = 4000 * 1024^2)

# Import
import_scRNAseq <- function(dir_base, name){
  data.dir <- paste0("", dir_base)
  raw <- Read10X(data.dir = data.dir); colnames(raw) <- paste0(name, "-", colnames(raw))
  
  # import scrublet results
  singlets <- fread(paste0("scrublet_10x/", dir_base, ".scrub.tsv")) %>%
    data.frame() %>% dplyr::filter(!called) %>% pull(barcode)
  
  # Filter for singlet genes adn non-mitos
  raw <- raw[!grepl("^MT", rownames(raw)),paste0(name, "-", substr(singlets,1,16))]
  raw <- CreateSeuratObject(counts = raw,  min.cells = 3, min.features = 2000); raw$source <- name
  raw <- NormalizeData(raw)
  raw <- FindVariableFeatures(raw)
  raw
  
}

# Import data while filtering out doublets from scrublet
pbmc <- import_scRNAseq("healthy_pbmc_5k_nextgem", "Healthy")

# Do dimension reduction
pbmc <- ScaleData(pbmc, verbose = FALSE)

pbmc <- RunPCA(pbmc, verbose = FALSE)
pbmc <- RunUMAP(pbmc, dims = 1:30, verbose = FALSE)

pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 0.5)

DimPlot(pbmc, label = TRUE)


# Verify that there isn't a clsuter due to cell cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
pbmc <- CellCycleScoring(pbmc, s.features = s.genes, g2m.features = g2m.genes)

FeaturePlot(pbmc, features = "G2M.Score")
FeaturePlot(pbmc, features = "nFeature_RNA")


new.cluster.ids <- c(
  '0'='Naive_CD4Tcell',
  '1'='CD14_Monocytes',
  '2'='Memory_CD4Tcell',
  '3'='Effector_CD8Tcell',
  '4'='NK_cell',
  '5'='Activated_Bcell',
  '6'="Memory_Bcell",
  '7'='CD14_Monocytes',
  '8'='Memory_CD8Tcell',
  '9' = "Dendritic_cells",
  '10' = "NK_cell",
  '11' = 'pDC',
  '12' = 'Progenitors'
)

Idents(pbmc) <- new.cluster.ids[as.character(pbmc@meta.data$seurat_clusters)]
pbmc$celltype <- Idents(pbmc)
saveRDS(pbmc, "output/PBMC5k_scRNAseq-for-integration.rds")
