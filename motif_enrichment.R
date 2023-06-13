#"OB"  "T2D" "LN"  "OW"

seurat_obj <- readRDS("t2d_final.rds")

DefaultAssay(seurat_obj) <- 'chromvar'

celltypes_vector =  c(
  "Activated stellate", 
  "Quiescent stellate",
  "Beta1",
  "Beta2",
  "Beta3",
  "Beta4",
  "Beta5",
  "Beta6"
)

celltype_col_name = "celltype_final" 
sample_group_col_name = "sample_group"

logfc.threshold = 0.15
p.value = 0.05

ident1_char = "OB"
ident2_char = "LN"

for (celltype in celltypes_vector){
  
  # Set identity classes to celltype_col_name
  Idents(object = seurat_obj) <- celltype_col_name
  
  # Subset Seurat object based on identity class
  seurat_obj_subset = subset(x = seurat_obj, idents = celltype)
  
  # Set identity classes to sample_group_col_name
  Idents(object = seurat_obj_subset) <- sample_group_col_name
  
  differential.activity <- FindMarkers(
    object = seurat_obj_subset,
    ident.1 = ident1_char,
    ident.2 = ident2_char,
    only.pos = F,
    mean.fxn = rowMeans,
    fc.name = "avg_diff",
    logfc.threshold = 0
  )
  
  differential.activity$motif.id = rownames(differential.activity)
  motif.id.name.map = read.csv("Motif.ID.Name.Map.csv")
  differential.activity = merge(differential.activity, motif.id.name.map, 
                                by = "motif.id", all.x = T)
  differential.activity = arrange(differential.activity, p_val)
  
  differential.activity$keyvals = "grey30"
  differential.activity$keyvals[differential.activity$p_val <= p.value & differential.activity$avg_diff >= logfc.threshold] = "#ff7f0e"
  differential.activity$keyvals[differential.activity$p_val <= p.value & differential.activity$avg_diff <= -logfc.threshold] = "#2ca02c"
  
  differential.activity = arrange(differential.activity, keyvals)
  
  keyvals = differential.activity$keyvals
  names(keyvals)[keyvals == 'grey30'] <- 'NS'
  names(keyvals)[keyvals == '#ff7f0e'] <- 'UP'
  names(keyvals)[keyvals == '#2ca02c'] <- 'DOWN'
  
  #keyvals = factor(keyvals, levels = c("grey30","#ff7f0e","#2ca02c")) # ADD here
  
  p = volcano.plot(d.mtx = differential.activity,
                   lab = NA,
                   #selectLab = c(differential.activity.pos$motif.name[1:15], differential.activity.neg$motif.name[1:15]),
                   x = "avg_diff",
                   y = "p_val",
                   p.value = p.value,
                   logfc.threshold = logfc.threshold,
                   celltype = celltype,
                   colCustom = keyvals
  )
  
  
  svglite(paste0(celltype,"_",paste0(ident1_char,collapse = "_"),"vs",paste0(ident2_char,collapse = "_"),"_Significant_Motifs (P<=",p.value,"&LogFC>=",logfc.threshold,").svg"), width = 10, height = 10)
  plot(p)
  dev.off()
  
}


# volcano.plot function
volcano.plot <- function(d.mtx = differential.activity,
                         lab = NA,
                         #selectLab = c(differential.activity.pos$motif.name[1:15], differential.activity.neg$motif.name[1:15]),
                         x = "avg_diff",
                         y = "p_val",
                         p.value = p.value,
                         logfc.threshold = logfc.threshold,
                         celltype = celltype,
                         colCustom = keyvals){
  
  p = EnhancedVolcano(d.mtx,
                      lab = NA,
                      #selectLab = selectLab,
                      x = x,
                      y = y,
                      colCustom = colCustom,
                      #legendLabels=c('NS','Down','UP'), # ADD here
                      pCutoff = p.value,
                      FCcutoff = logfc.threshold,
                      pointSize = 5,
                      labSize = 6,
                      #drawConnectors = TRUE,
                      #endsConnectors = "last",
                      boxedLabels = F,
                      max.overlaps = Inf,
                      #title = paste0(celltype,"_",ident1_char,"vs",ident2_char,"_Significant_Motifs (P<=",p.value,"&LogFC>=",logfc.threshold,")"),
                      subtitle = NULL,
                      caption = NULL,
                      titleLabSize = 20,
                      legendLabSize = 18,
                      title = NULL
  )+
    theme(plot.title=element_text(hjust=0.5))
  
  return(p)
}