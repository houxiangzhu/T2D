# READ saved t2d poject
t2d <- readRDS("t2d_final.rds")

DefaultAssay(t2d) <- 'chromvar'

motif.id.name.map = read.csv("Motif.ID.Name.Map.csv")

motif.id.name.map = motif.id.name.map[match(rownames(t2d@assays$chromvar@data),
                                            motif.id.name.map$motif.id),]

rownames(t2d@assays$chromvar@data) = motif.id.name.map$motif.name

Idents(object = t2d) <- "celltype_final"

t2d = subset(x = t2d, idents = c("Beta1", "Beta2", "Beta3", "Beta4", "Beta5", "Beta6"))

all.markers <- FindAllMarkers(object = t2d, 
                              only.pos = T,
                              mean.fxn = rowMeans,
                              fc.name = "avg_diff",
                              logfc.threshold = 0
)

rownames(all.markers) = all.markers$X
all.markers$X = NULL

all.markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_diff) -> top

top$cluster = as.character(top$cluster)

top %>%
  arrange(cluster, desc(avg_diff)) -> top

write.csv(top, file = "beta.pos.markers.top20.csv")

# scale the data
t2d_scaled = ScaleData(object = t2d,
                       features = rownames(t2d@assays$chromvar@data)
)
Idents(t2d_scaled) <- factor(x = Idents(t2d_scaled), levels = sort(unique(t2d_scaled$celltype_final)))


### captures count of cells for the ident with the fewest
maxcells  <- min(table(Idents(t2d_scaled)))

pdf("betas.only.markers.heatmap.downsampled.top20.pdf", width = 30, height = 10)

DoHeatmap(subset(t2d_scaled, downsample = maxcells), 
          features = top$gene,
          angle = 45) +
  labs(fill = "Z-score") +
  theme(text=element_text(size=30))

dev.off()