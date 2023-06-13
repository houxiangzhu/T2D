library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
library(motifmatchr)
library(ggseqlogo)
library(chromVAR)

## Motif analysis with Signac
### loading and viewing data
setwd("/Users/biocore/Documents/Projects/scATACseq/CRI-BIO-842")
t2d <- readRDS("Full29_ATAC_integration_majorityvoting_060222.rds")
#t2d_2 <- readRDS("Full29_ATAC_integration_majorityvoting_060222.rds")

# p1 <- DimPlot(t2d, 
#               group.by = "Majority_RNA20",
#               label = TRUE, 
#               pt.size = 0.1,
#               label.size = 4,
#               repel = T)
# 
# pdf("majority_voting_20_dimplot.pdf", width = 10, height = 8)
# p1
# dev.off()

### Removing peaks that are on scaffolds

main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)

keep.peaks <- as.logical(seqnames(granges(t2d)) %in% main.chroms)

t2d <- t2d[keep.peaks, ]

### Adding motif information to the Seurat object

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
t2d <- AddMotifs(
  object = t2d,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)


## Computing motif activities by chromVAR

t2d <- RunChromVAR(
  object = t2d,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

DefaultAssay(t2d) <- 'chromvar'

save(t2d, file = "chromVAR.Rdata")





