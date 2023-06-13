# READ saved t2d object
t2d <- readRDS("t2d_final.rds")

Idents(t2d) <- "celltype_final"

# Subset Seurat object based on identity class
t2d_subset = subset(x = t2d, idents = "Beta6")

# Set identity classes to sample_group_col_name
Idents(t2d_subset) <- "sample_group"

# Downsample the number of cells per identity class
#set.seed(8)
t2d_subset = subset(x = t2d_subset, downsample = 300)

# split fragmetn file to LN, OW, OB, T2D by cell identities

SplitFragments(
  t2d_subset,
  outdir = "fragment_files"
)