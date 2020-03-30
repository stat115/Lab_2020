library(data.table)
library(dplyr)

dt <- fread("atac_pbmc_500_nextgem_fragments.tsv.gz", header = FALSE)
bcs <- fread("atac_pbmc_500_data/barcodes.tsv.gz", header = FALSE)[[1]] %>% as.character()
dt <- dt[dt$V4 %in% bcs,]
write.table(dt, file = "atac_pbmc_500_data/fragments.tsv", sep ="\t", quote =FALSE, row.names = FALSE, col.names = FALSE)
