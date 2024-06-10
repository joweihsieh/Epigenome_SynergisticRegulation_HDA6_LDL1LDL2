computeMatrix scale-regions -p 20 -S ChIP_H3Ac_WT_1.bw ChIP_H3Ac_WT_2.bw ChIP_H3Ac_hda6_1.bw ChIP_H3Ac_hda6_2.bw ChIP_H3Ac_ldl12_1.bw ChIP_H3Ac_ldl12_2.bw ChIP_H3Ac_hda6ldl12_1.bw ChIP_H3Ac_hda6ldl12_2.bw -R tair10_gene2.bed -p 20 -bs 100 -a 3000 -b 3000 -m 6000 -out Gene_H3Ac_3kf.matrix.gz
plotHeatmap -m Gene_H3Ac_3kf.matrix.gz -out Gene_H3Ac_3kf.png --heatmapWidth 10 --colorMap Reds 

computeMatrix scale-regions -p 20 -S ChIP_H3K4me2_WT_1.bw ChIP_H3K4me2_WT_2.bw ChIP_H3K4me2_hda6_1.bw ChIP_H3K4me2_hda6_2.bw ChIP_H3K4me2_ldl12_1.bw ChIP_H3K4me2_ldl12_2.bw ChIP_H3K4me2_hda6ldl12_1.bw ChIP_H3K4me2_hda6ldl12_2.bw -R tair10_gene2.bed -p 20 -bs 100 -a 3000 -b 3000 -m 6000 -out Gene_H3K4me2_3kf.matrix.gz
plotHeatmap -m Gene_H3K4me2_3kf.matrix.gz -out Gene_H3K4me2_3kf.png --heatmapWidth 10 --colorMap Reds


computeMatrix scale-regions -p 20 -S ATAC_WT_150bp_1.bw ATAC_WT_150bp_2.bw ATAC_hda6_150bp_1.bw ATAC_hda6_150bp_2.bw ATAC_ldl12_150bp_1.bw ATAC_ldl12_150bp_2.bw ATAC_hda6ldl12_150bp_1.bw ATAC_hda6ldl12_150bp_2.bw -R tair10_gene2.bed -p 20 -bs 100 -a 3000 -b 3000 -m 6000 -out Gene_ATAC_3kf.matrix.gz
plotHeatmap -m Gene_ATAC_3kf.matrix.gz -out Gene_ATAC_3kf.png --heatmapWidth 10 --colorMap Reds

