# The procedure mostly follows Seurat tutorial: https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html
library(Seurat)
library(Signac)

#setwd
setwd("/Users/biocore/Documents/Projects/scATACseq/CRI-BIO-842/RNA_Velocity")

# read in t2d obj
seurat_obj <- readRDS("../t2d_final.rds")

# subset for the same cells in the jointly filtered anndata object
filtered.cells = read.table("filtered_cells.txt", header = F)
filtered.cells = filtered.cells$V1

seurat_obj_subset = subset(seurat_obj, cells = filtered.cells)

# preprocess RNA
DefaultAssay(seurat_obj_subset) <- "RNA"
seurat_obj_subset <- NormalizeData(seurat_obj_subset)
seurat_obj_subset <- FindVariableFeatures(seurat_obj_subset)
seurat_obj_subset <- ScaleData(seurat_obj_subset, do.scale = F) # not scaled for consistency with scVelo (optionally, use SCTransform)
seurat_obj_subset <- RunPCA(seurat_obj_subset, verbose = FALSE)
seurat_obj_subset <- RunUMAP(seurat_obj_subset, dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_') # optional

# preprocess ATAC
DefaultAssay(seurat_obj_subset) <- "peaks"
seurat_obj_subset <- RunTFIDF(seurat_obj_subset)
seurat_obj_subset <- FindTopFeatures(seurat_obj_subset, min.cutoff = 'q0')
seurat_obj_subset <- RunSVD(seurat_obj_subset)
seurat_obj_subset <- RunUMAP(seurat_obj_subset, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_") # optional

# find weighted nearest neighbors
seurat_obj_subset <- FindMultiModalNeighbors(seurat_obj_subset, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50), k.nn = 50)
seurat_obj_subset <- RunUMAP(seurat_obj_subset, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_") # optional

# extract neighborhood graph
nn_idx <- seurat_obj_subset@neighbors$weighted.nn@nn.idx
nn_dist <- seurat_obj_subset@neighbors$weighted.nn@nn.dist
nn_cells <- seurat_obj_subset@neighbors$weighted.nn@cell.names

# save neighborhood graph
write.table(nn_idx, "nn_idx.txt", sep = ',', row.names = F, col.names = F, quote = F)
write.table(nn_dist, "nn_dist.txt", sep = ',', row.names = F, col.names = F, quote = F)
write.table(nn_cells, "nn_cells.txt", sep = ',', row.names = F, col.names = F, quote = F)