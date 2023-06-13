library(dplyr)

genes = read.csv("CI_genes.txt", header = F)
genes = genes$V1

tss = read.csv("gencode.v34.annotation.genes.tss_flanking2000.csv")

tss.selected = tss %>%
  filter(gene_name %in% genes)


tss.bed = tss.selected %>%
  mutate(start = start -1) %>%
  dplyr::select(seqnames, start, end)

write.table(tss.bed, file = "CI_genes_tss.bed", col.names = F, row.names = F, quote = F, sep = "\t")

# generate bed file for the six CI genes only

#library(dplyr)

genes = read.csv("CI_genes_six.txt", header = F)
genes = genes$V1

tss = read.csv("gencode.v34.annotation.genes.tss_flanking2000.csv")

tss.selected = tss %>%
  filter(gene_name %in% genes)


tss.bed = tss.selected %>%
  mutate(start = start -1) %>%
  dplyr::select(seqnames, start, end)

write.table(tss.bed, file = "CI_genes_six_tss.bed", col.names = F, row.names = F, quote = F, sep = "\t")