
# Merge each fragment bed file using bedtools merge to generate one merged bedgraph file
bedtools merge -i LN.bed -c 5 -o sum > LN.merged.bed

bedtools merge -i OB.bed -c 5 -o sum > OB.merged.bed

bedtools merge -i OW.bed -c 5 -o sum > OW.merged.bed

bedtools merge -i T2D.bed -c 5 -o sum > T2D.merged.bed


# Sort the merged bed file

bedtools sort -i LN.merged.bed > LN.merged.sorted.bed

bedtools sort -i OW.merged.bed > OW.merged.sorted.bed

bedtools sort -i OB.merged.bed > OB.merged.sorted.bed

bedtools sort -i T2D.merged.bed > T2D.merged.sorted.bed

