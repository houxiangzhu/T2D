# READ saved t2d object
t2d <- readRDS("t2d_final.rds")

#### change the gene annotation to hg19
## Add gene annotation from a gtf file
gtf <- rtracklayer::import('gencode.v19.annotation.gtf')

gtf$tx_id = gtf$transcript_id

gtf$gene_biotype = gtf$gene_type

seqlevelsStyle(gtf) <- "UCSC"

# set gene annotations
Annotation(t2d) <- gtf

# get gene annotation information
#Annotation(t2d)
#### done add hg19 annotation


Idents(t2d) = "celltype_final"

t2d_subset = subset(t2d, idents = c("Beta1", "Beta2", "Beta3",
                                    "Beta4", "Beta5", "Beta6"))

#########
# create a new celltype_group column
t2d_subset@meta.data$celltype_group = paste(t2d_subset@meta.data$celltype_final, t2d_subset@meta.data$sample_group, sep = "_")

# set levels
pheno = c("LN", "OW", "OB", "T2D")
celltypes = sort(unique(t2d_subset$celltype_final))

ordered.level = c()

for (c in celltypes) {
  for (p in pheno){
    c_p = paste(c,p,sep = "_")
    ordered.level = append(ordered.level, c_p)
  }
  
}

# Set identity classes to celltype_col_name
Idents(object = t2d_subset) <- "celltype_group"

# sort levels
Idents(t2d_subset) <- factor(x = Idents(t2d_subset), levels = ordered.level)

#############
# add links for hic data 2

hi.c.df = read.csv("./HIC2/data/results/GSM5677610_Beta_1.2000_pt1.tsv", header = F, skip = 1, sep = "\t")

hi.c.df.filtered = hi.c.df %>%
  dplyr::filter(V1 == V4) %>%
  dplyr::filter(V7 <= 0.1) %>%
  mutate(chr = V1,
         start = ifelse((V2+V3)/2 <= (V5+V6)/2, round((V2+V3)/2), round((V5+V6)/2)),
         end = ifelse((V2+V3)/2 <= (V5+V6)/2, round((V5+V6)/2), round((V2+V3)/2)),
         score = V8) %>%
  dplyr::select(chr,start,end, score)

link.granges = makeGRangesFromDataFrame(hi.c.df.filtered, keep.extra.columns=T)


Links(t2d_subset) <- link.granges

# Plot tracking plot
set.seed(2)

genes = c("PFKFB2", "ATP2A3", "CHL1", "GLRA1", "RASGRF1", 
          "SLC2A2", "WSCD2", "P2RY1", "SCARB1", "ANKRD27")

for (gene in genes){
  
  p = CoveragePlot(
    object = t2d_subset,
    region = gene,
    annotation = T,
    peaks = F,
    extend.upstream = 500000,
    extend.downstream = 500000, 
    features = gene,
    expression.assay = "RNA",
    expression.slot = "data",
    links = TRUE
  )
  
  svglite(paste0(gene,"_HiC_2.svg"), width = 20, height = 10)
  print(p)
  dev.off()
  
  
  
}