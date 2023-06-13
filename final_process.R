load(file = "chromVAR.Rdata")

### Add the final cell types label column
t2d@meta.data$celltype_final = t2d@meta.data$Majority_RNA10

rest = !(t2d@meta.data$celltype_final %in% c("Ductal", "Delta", "Acinar"))

t2d@meta.data$celltype_final[rest] = t2d@meta.data$Yan_celltype[rest]

save(t2d, file = "t2d_final_withNoneGroup.Rdata")

# Set identity classes to an existing column in meta data
Idents(object = t2d) <- "celltype_final"

# Subset Seurat object based on identity class
t2d = subset(x = t2d, idents = c("None"), invert = TRUE)

# t2d_acinar = subset(x = t2d, idents = c("Acinar"))

# add the groups info
ln = c("HI_17", "HI_19", "HI_8", "HI_11", "HI_5", "HI_21", "HI_25")
ow = c("HI_7","HI_31","HI_9","HI_18","HI_20")
ob = c("HI_15","HI_3","HI_6","HI_22","HI_1","HI_26","HI_32")
T2D = c("HI_16","HI_4","HI_23","HI_24","HI_12","HI_14","HI_2","HI_10","HI_27","HI_28")

t2d@meta.data$sample_group = t2d@meta.data$orig.ident

ln_in = (t2d@meta.data$sample_group %in% ln)
ow_in = (t2d@meta.data$sample_group %in% ow)
ob_in = (t2d@meta.data$sample_group %in% ob)
T2D_in = (t2d@meta.data$sample_group %in% T2D)

t2d@meta.data$sample_group[ln_in] = "LN"
t2d@meta.data$sample_group[ow_in] = "OW"
t2d@meta.data$sample_group[ob_in] = "OB"
t2d@meta.data$sample_group[T2D_in] = "T2D"

# save the final object
#save(t2d, file = "t2d_final.Rdata")
saveRDS(t2d, file = "t2d_final.rds")

p <- DimPlot(t2d, 
             group.by = "celltype_final",
             label = TRUE, 
             pt.size = 0.1,
             label.size = 4,
             repel = T)

p1 <- DimPlot(t2d, 
              group.by = "sample_group",
              label = TRUE, 
              pt.size = 0.1,
              label.size = 4,
              repel = T)

pdf("t2d_final_dimplot.pdf", width = 20, height = 8)
p+p1
dev.off()