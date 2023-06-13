cd /Users/biocore/Documents/Projects/scATACseq/CRI-BIO-842/fragment_files

conda activate deeptools

computeMatrix reference-point --referencePoint center -b 10000 -a 10000 -R CI_genes_tss.bed -S LN.bw OW.bw OB.bw T2D.bw --binSize 200 -o All_TSS_Binsum.gz --averageTypeBins sum


plotHeatmap -m All_TSS_Binsum.gz -out All_TSS_Binsum.10kb.pdf --refPointLabel TSS --regionsLabel "Beta6 CI genes" --outFileSortedRegions Heatmap1sortedRegions.Binsum.bed --colorList 'white,#ffe0b1,#ff7403,#f70604' --missingDataColor 1

plotProfile -m All_TSS_Binsum.gz -out All_TSS_Binsum.10kb.profile.pdf --perGroup --colors "#1B9E77" "#D95F02" "#7570B3" "#E7298A" --refPointLabel TSS --regionsLabel ""
