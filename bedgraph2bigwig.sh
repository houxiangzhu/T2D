# Convert the merged sorted bedgraph file to bigwig file using ucsc's bedgraphtobigwig command

conda activate bedgraphtobigwig

bedGraphToBigWig LN.merged.sorted.bed hg38.chrom.sizes LN.bw

bedGraphToBigWig OW.merged.sorted.bed hg38.chrom.sizes OW.bw

bedGraphToBigWig OB.merged.sorted.bed hg38.chrom.sizes OB.bw

bedGraphToBigWig T2D.merged.sorted.bed hg38.chrom.sizes T2D.bw