library(data.table)
library(dplyr)
library(Seurat)
library(Signac)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)

# Accessory function to import the peaks x cells count matrix
import_counts_mat_10x <- function(directory){
  peaks <- fread(paste0(directory, "/peaks.bed.gz"), header = FALSE, col.names = c("chr", "start", "end"))
  peak_character <- paste0(peaks$chr, ":", as.character(peaks$start), "-", as.character(peaks$end))
  barcodes <- fread(paste0(directory, "/barcodes.tsv.gz"), header = FALSE)[[1]]
  mtx <- fread(paste0(directory, "/matrix.mtx.gz")[[1]], skip = 3, header = FALSE)
  
  # Assemble binary matrix
  mat <- Matrix::sparseMatrix(i = c(mtx[[1]],length(peaks)), 
                              j = c(mtx[[2]],length(barcodes)),
                              x = c(as.numeric(mtx[[3]]) > 0, 0))
  colnames(mat) <- barcodes
  rownames(mat) <- peak_character
  mat
}

counts <- import_counts_mat_10x("atac_pbmc_500_data")

# Import single-cell meta data
metadata <- read.csv(
  file = "atac_pbmc_500_data/atac_pbmc_500_nextgem_singlecell.csv.gz",
  header = TRUE,
  row.names = 1
)

# Create a Seurat object
pbmc <- CreateSeuratObject(
  counts = counts,
  assay = 'peaks',
  project = 'ATAC',
  min.cells = 1,
  meta.data = metadata
)

# Set fragments path
fragment.path <- "atac_pbmc_500_data/fragments.tsv.gz"
pbmc <- SetFragments(
  object = pbmc,
  file = fragment.path
)

# There are many different QC metrics to filter cells; here, we will only do the total abundance of peaks
# For more options, check out https://satijalab.org/signac/articles/pbmc_vignette.html

# The passed_filters is the number of fragments that is HQ mapped with no duplicates, etc. 
pbmc <- subset(pbmc, subset = passed_filters > 5000)


#---- 
# Create gene activity scores
#----
gene.coords <- genes(EnsDb.Hsapiens.v75, filter = ~ gene_biotype == "protein_coding")
seqlevelsStyle(gene.coords) <- 'UCSC'
genebody.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')
genebodyandpromoter.coords <- Extend(x = gene.coords, upstream = 2000, downstream = 0)

# create a gene by cell matrix-- this step takes ~ 2 minutes
gene.activities <- FeatureMatrix(
  fragments = fragment.path,
  features = genebodyandpromoter.coords,
  cells = colnames(pbmc),
  chunk = 20
)

# convert rownames from chromsomal coordinates into gene names
gene.key <- genebodyandpromoter.coords$gene_name
names(gene.key) <- GRangesToString(grange = genebodyandpromoter.coords)
rownames(gene.activities) <- gene.key[rownames(gene.activities)]

# add the gene activity matrix to the Seurat object as a new assay, and normalize it
pbmc[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(
  object = pbmc,
  assay = 'ACTIVITY',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc$nCount_ACTIVITY)
)

#---- 
# Do standard scATAC-seq dimension reduction
#----
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(
  object = pbmc,
  assay = 'peaks',
  reduction.key = 'LSI_',
  reduction.name = 'lsi', n = 20
)

# Visualize correlation between first component and number of detected peaks
FeatureScatter(pbmc, 'LSI_1', 'nCount_peaks')
FeatureScatter(pbmc, 'LSI_2', 'nCount_peaks')
FeatureScatter(pbmc, 'LSI_3', 'nCount_peaks')

# Finalize dimension reduction
pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:20)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:20)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)
DimPlot(object = pbmc, label = TRUE) + NoLegend()

# Visualize gene scores
DefaultAssay(pbmc) <- 'ACTIVITY'

# Visualize marker genes
FeaturePlot(
  object = pbmc,
  features = c('MS4A1', 'CD3D', 'NKG7', 'TREM1'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)


# Load the pre-processed scRNA-seq data for PBMCs
pbmc_rna <- readRDS("scrna_source/output/PBMC5k_scRNAseq-for-integration.rds")

# Perform label transfer using the scRNA as a reference
transfer.anchors <- FindTransferAnchors(
  reference = pbmc_rna,
  query = pbmc,
  reduction = 'cca'
)

# Predict cell type labels in ATAC space
predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = pbmc_rna$celltype,
  weight.reduction = pbmc[['lsi']],
  dims = 2:20
)

pbmc <- AddMetaData(object = pbmc, metadata = predicted.labels)
DimPlot(
  object = pbmc_rna,
  group.by = 'celltype',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

DimPlot(
  object = pbmc,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')



#---------------------------------------------------------------------------------------------

#---- 
# For graaduate students
#----


# Update identities with transfered labels
Idents(pbmc) <- "predicted.id"

# Find differentially accessible peaks between B-cells and T-cells
DefaultAssay(pbmc) <- 'peaks'

da_peaks <- FindMarkers(
  object = pbmc,
  ident.1 = "Memory_Bcell",
  ident.2 = "Naive_CD4Tcell",
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)

# Visualize chromatin accessibility tracks
CoveragePlot(
  object = pbmc,
  region = rownames(da_peaks)[c(1,2)],
  sep = c(":", "-"),
  peaks = StringToGRanges(rownames(pbmc), sep = c(":", "-")),
  annotation = EnsDb.Hsapiens.v75,
  extend.upstream = 20000,
  extend.downstream = 20000,
  ncol = 1
)

library(JASPAR2018)  # BiocManager::install("JASPAR2018")
library(TFBSTools) # BiocManager::install("TFBSTools")
library(motifmatchr) # BiocManager::install("motifmatchr")
library(BSgenome.Hsapiens.UCSC.hg19) # BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2018,
  opts = list(species = "Homo sapiens", all_versions = FALSE)
)

# Great a motif matrix
motif.matrix <- CreateMotifMatrix(
  features = StringToGRanges(rownames(pbmc), sep = c(":", "-")),
  pwm = pfm,
  genome = 'hg19',
  sep = c(":", "-"),
  use.counts = FALSE
)

# Integrate it into the Seurat object
motif <- CreateMotifObject(
  data = motif.matrix,
  pwm = pfm
)

# Add the Motif object to the assay
pbmc[['peaks']] <- AddMotifObject(
  object = pbmc[['peaks']],
  motif.object = motif
)

pbmc <- RegionStats(
  object = pbmc,
  genome = BSgenome.Hsapiens.UCSC.hg19,
  sep = c(":", "-")
)

# Find overall enriched motifs
enriched.motifs <- FindMotifs(
  object = pbmc,
  features = head(rownames(da_peaks), 1000)
)

# Here's the underlying DNA sequence for the top motifs
MotifPlot(
  object = pbmc,
  motifs = head(rownames(enriched.motifs))
)

# Score single cells
pbmc <- RunChromVAR(
  object = pbmc,
  genome = BSgenome.Hsapiens.UCSC.hg19
)

DefaultAssay(pbmc) <- 'chromvar'

# look at the activity of SPI1, the top hit from the overrepresentation testing
FeaturePlot(
  object = pbmc,
  features = rownames(enriched.motifs)[[1]],
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)

