#ATAC-seq integration analysis code
#tmp Ren
#06-07-2022

#the code shows how to analyze the data saved at sftp://gardner.cri.uchicago.edu/group/bioinformatics/Projects/CRI-BIO-842-HMS-Kulkarni-yli-zyren/ATAC_seq_integration_sharing_with_Houxiang/Full29_ATAC_integration_majorityvoting_060222.rds
library(Seurat)
library(Signac)
library(harmony)
library(distances)

#filtered_combined_051022.rds files are in sftp://gardner.cri.uchicago.edu/group/bioinformatics/shared/projectReports/CRI-BIO-842-owiewieisdsafaogiadfaofdosaa/ATAC_seq_processed_rds_data_29sample_060222.zip
for (i in c(1:32)){
  tmp<-readRDS(paste0("HI_",i,"_filtered_combined_051022.rds"))
  combined <- merge(
    x = combined,
    y = tmp
  )
  rm(tmp)
}
meta<-combined@meta.data
combined <- merge(
  x = HI_1,
  y = list(HI_2,HI_3,HI_4,HI_5,HI_6,HI_7,HI_8,HI_9,HI_10,HI_11,HI_12,HI_13,HI_14,HI_15,HI_16,HI_17,HI_18,HI_19,HI_20,HI_21,HI_22,HI_23,HI_24,HI_25,HI_26,HI_27,HI_28,HI_31,HI_32),
  add.cell.ids = c("HI_1","HI_2","HI_3","HI_4","HI_5","HI_6","HI_7","HI_8","HI_9","HI_10","HI_11","HI_12","HI_13","HI_14","HI_15","HI_16","HI_17","HI_18","HI_19","HI_20","HI_21","HI_22","HI_23","HI_24","HI_25","HI_26","HI_27","HI_28","HI_31","HI_32")
)



combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 50)
combined <- RunSVD(combined, assay = "peaks")

table(meta$orig.ident)
combined <- RunHarmony(
  object = combined,
  group.by.vars = 'orig.ident',
  reduction = 'lsi',
  assay.use = 'peaks',
  max.iter.harmony = 20,
  project.dim = FALSE
)

combined <- FindNeighbors(object = combined, reduction="harmony", dims = 2:30)
combined <- FindClusters(object = combined, verbose = FALSE, algorithm = 2, resolution = 0.6)

combined <- RunUMAP(object = combined, reduction = 'harmony', dims = 2:15)
DimPlot(combined, group.by = 'orig.ident', pt.size = 0.1) + ggplot2::ggtitle("Harmony integration\ncolor by sample HI1 to 32 (no 13,29,30)")

######################################
#
#label the cell type using Yan's label
######################################
Yan_RNA<-readRDS(file="~/full29samp_rmSamp13_29_30_integration_results.rds")
#data for Yan_RNA is at location:sftp://gardner.cri.uchicago.edu/gpfs/data/biocore-analysis/CRI-BIO-842-HMS-Kulkarni-yli-zyren/integration_29samples_rm13_29_30
Yan_RNA_meta<-Yan_RNA@meta.data
meta_ATAC<-combined@meta.data
table(Yan_RNA_meta$orig.ident)
Yan_RNA_meta$orig.ident[Yan_RNA_meta$orig.ident=="HI10"]<-"HI_10"
Yan_RNA_meta$orig.ident[Yan_RNA_meta$orig.ident=="HI5"]<-"HI_5"
Yan_RNA_meta$orig.ident[Yan_RNA_meta$orig.ident=="HI8"]<-"HI_8"
Yan_RNA_meta$orig.ident[Yan_RNA_meta$orig.ident=="HI9"]<-"HI_9"
tmp<-combined@meta.data
tmp$RowName2<-paste0(tmp$orig.ident,"_",substr(row.names(tmp),0,18))
Yan_RNA_meta<-Yan_RNA_meta[paste0(Yan_RNA_meta$orig.ident,"_",substr(row.names(Yan_RNA_meta),0,18)) %in% tmp$RowName2, ]
tmp<-tmp[match(paste0(Yan_RNA_meta$orig.ident,"_",substr(row.names(Yan_RNA_meta),0,18)),tmp$RowName2),]

#make sure the number matches.
sum(paste0(Yan_RNA_meta$orig.ident,"_",substr(row.names(Yan_RNA_meta),0,18))==tmp$RowName2)


#ask Yan for this txt file.
Yan_cluster_name <- read.delim("~/finalizedCellClusters_rmOutlier_6betaSubCluster_4tmp.txt")

Yan_cluster_name<-Yan_cluster_name[match(tmp$Yan_celllabel,Yan_cluster_name$X),]
sum(is.na(tmp$Yan_celllabel==Yan_cluster_name$X))
tmp$Yan_celllabel[is.na(tmp$Yan_celllabel==Yan_cluster_name$X)][1]
tmp$Yan_celllabel[2]==Yan_cluster_name$X[2]
head(tmp$Yan_celllabel==Yan_cluster_name$X)
Yan_cluster_name$cellClusterNames[is.na(tmp$Yan_celllabel==Yan_cluster_name$X)][1]

sum(table(Yan_cluster_name$cellClusterNames))
sum(is.na(Yan_cluster_name$cellClusterNames))
Yan_cluster_name$cellClusterNames[is.na(Yan_cluster_name$cellClusterNames)]<-"None"
tmp$Yan_celltype<-Yan_cluster_name$cellClusterNames

tmp2<-tmp[!(tmp$Yan_celllabel %in% Yan_cluster_name$X),]

combined@meta.data<-tmp


cells<-row.names(tmp)[tmp$Yan_celltype!="None"]
DimPlot(combined,cells = cells,pt.size = 0.1,label = T, group.by = 'Yan_celltype') + ggplot2::ggtitle("29 Sample ATAC Harmony Integration\ncolor by Updated RNA Cell Type")


DimPlot(combined,pt.size = 0.1,label = T, group.by = 'Yan_celltype') + ggplot2::ggtitle("29 Sample ATAC Harmony Integration\ncolor by Updated RNA Cell Type")

tmp<-combined@meta.data

tmp$Type[tmp$orig.ident %in% c("HI_5","HI_8","HI_11","HI_17","HI_19","HI_21","HI_25")]<-"Lean-ND"
tmp$Type[tmp$orig.ident %in% c("HI_1","HI_3","HI_6","HI_13","HI_15","HI_22","HI_26")]<-"Obese-ND"
tmp$Type[tmp$orig.ident %in% c("HI_7","HI_9","HI_18","HI_20")]<-"Overweight-ND"
tmp$Type[tmp$orig.ident %in% c("HI_2","HI_4","HI_10","HI_12","HI_14","HI_16","HI_23","HI_24","HI_27","HI_28")]<-"T2D"
tmp$Type<-factor(tmp$Type,levels=c("Lean-ND","Overweight-ND","Obese-ND","T2D"))
combined@meta.data<-tmp


p1 <- DimPlot(combined, group.by = 'orig.ident', pt.size = 0.1) + ggplot2::ggtitle("Harmony integration")

p2 <- DimPlot(combined,pt.size = 0.1,label = T, group.by = 'Type') + ggplot2::ggtitle("Harmony integration")

ggarrange(p1,p2)


##################################################
#calculate the sample using 10 nearest neigbhor voting
###################################################
UMAP<-combined@reductions$umap@cell.embeddings
dim(UMAP)
head(row.names(UMAP))
colnames(UMAP)
meta<-combined@meta.data
meta<-merge(meta,UMAP,by="row.names")
head(row.names(meta))

# Euclidean distances
my_distances1 <- distances::distances(my_data_points)
# Euclidean distances in only one dimension
UMAP<-as.data.frame(UMAP)
my_distances <-distances::distances(meta,
                                    dist_variables = c("UMAP_1","UMAP_2"))
min(my_distances[,1])
major_RNA_type<-as.data.frame(matrix(NA,nrow=172030,ncol=2))
major_RNA_type10<-major_RNA_type

for (i in (1:172030)){
  tmp<-NULL
  tmp2<-NULL
  order<-NULL
  tmp<-my_distances[,i]
  order<-names(tmp[order(tmp)[1:11]])[2:51]
  tmp2<-table(meta$Yan_celltype[as.numeric(order)])
  major_RNA_type[i,1]<-names(tmp2[tmp2==max(tmp2)])[1]
  major_RNA_type[i,2]<-max(tmp2)[1]
  
}



meta$Majority_RNA10<-major_RNA_type$V1
combined@meta.data<-meta

cells<-row.names(meta)[meta$Majority_RNA10!="None"]

DimPlot(combined,pt.size = 0.1,raster=TRUE,label = T, group.by = 'Majority_RNA10') + ggplot2::ggtitle("29 Sample ATAC Harmony Integration\ncolor by Majority Voting 10 neighbors RNA Cell Type")

DimPlot(combined,pt.size = 0.1,label = T, group.by = 'peaks_snn_res.0.6') + ggplot2::ggtitle("Harmony integration\ncolor by ATAC Res.0.3")

cells<-row.names(meta)[meta$Majority_RNA10!="None"]
DimPlot(combined,cells=cells,pt.size = 0.1,raster=TRUE,label = T, group.by = 'Majority_RNA10') + ggplot2::ggtitle("29 Sample ATAC Harmony Integration\ncolor by Majority Voting 10 neighbors RNA Cell Type")

table<-table(meta$Majority_RNA10,meta$Yan_celltype)

heatmap(table,Rowv = NULL, Colv =NULL ,display_numbers =F,scale="column",cluster_rows=F,cluster_cols=F,legend=F,xlab = "label from Majority 10 voting", ylab = "RNA-seq label")


row.names(meta)<-meta$Row.names

DimPlot(combined,pt.size = 0.1,ncol=6,raster=TRUE,label = F,group.by="orig.ident",split.by = "Majority_RNA10") + ggplot2::ggtitle("29 Sample ATAC Harmony Integration\ncolor by Majority Voting 10 neighbors RNA Cell Type")

saveRDS(combined,file="group/bioinformatics/Projects/CRI-BIO-842-HMS-Kulkarni-yli-zyren/ATAC_seq_integration_sharing_with_Houxiang/Full29_ATAC_integration_majorityvoting_060222.rds")
