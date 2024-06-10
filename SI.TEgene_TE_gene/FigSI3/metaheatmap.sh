computeMatrix scale-regions -p 20 -S ChIP_H3Ac_WT_1.bw ChIP_H3Ac_WT_2.bw ChIP_H3Ac_hda6_1.bw ChIP_H3Ac_hda6_2.bw ChIP_H3Ac_ldl12_1.bw ChIP_H3Ac_ldl12_2.bw ChIP_H3Ac_hda6ldl12_1.bw ChIP_H3Ac_hda6ldl12_2.bw -R tair10_TE_gene_expressed.bed -p 20 -bs 100 -a 500 -b 500 -m 2000 -out TE_gene_expressed_H3Ac.matrix.gz
plotHeatmap -m TE_gene_expressed_H3Ac.matrix.gz -out TE_gene_expressed_H3Ac.png --heatmapWidth 10 --colorMap Reds 

computeMatrix scale-regions -p 20 -S ChIP_H3K4me2_WT_1.bw ChIP_H3K4me2_WT_2.bw ChIP_H3K4me2_hda6_1.bw ChIP_H3K4me2_hda6_2.bw ChIP_H3K4me2_ldl12_1.bw ChIP_H3K4me2_ldl12_2.bw ChIP_H3K4me2_hda6ldl12_1.bw ChIP_H3K4me2_hda6ldl12_2.bw -R tair10_TE_gene_expressed.bed -p 20 -bs 100 -a 500 -b 500 -m 2000 -out TE_gene_expressed_H3K4me2.matrix.gz
plotHeatmap -m TE_gene_expressed_H3K4me2.matrix.gz -out TE_gene_expressed_H3K4me2.png --heatmapWidth 10 --colorMap Reds



computeMatrix scale-regions -p 20 -S ChIP_H3Ac_WT_1.bw ChIP_H3Ac_WT_2.bw ChIP_H3Ac_hda6_1.bw ChIP_H3Ac_hda6_2.bw ChIP_H3Ac_ldl12_1.bw ChIP_H3Ac_ldl12_2.bw ChIP_H3Ac_hda6ldl12_1.bw ChIP_H3Ac_hda6ldl12_2.bw -R tair10_TE_element_expressed.bed -p 20 -bs 100 -a 500 -b 500 -m 2000 -out TE_element_expressed_H3Ac.matrix.gz
plotHeatmap -m TE_element_expressed_H3Ac.matrix.gz -out TE_element_expressed_H3Ac.png --heatmapWidth 10 --colorMap Reds 

computeMatrix scale-regions -p 20 -S ChIP_H3K4me2_WT_1.bw ChIP_H3K4me2_WT_2.bw ChIP_H3K4me2_hda6_1.bw ChIP_H3K4me2_hda6_2.bw ChIP_H3K4me2_ldl12_1.bw ChIP_H3K4me2_ldl12_2.bw ChIP_H3K4me2_hda6ldl12_1.bw ChIP_H3K4me2_hda6ldl12_2.bw -R tair10_TE_element_expressed.bed -p 20 -bs 100 -a 500 -b 500 -m 2000 -out TE_element_expressed_H3K4me2.matrix.gz
plotHeatmap -m TE_element_expressed_H3K4me2.matrix.gz -out TE_element_expressed_H3K4me2.png --heatmapWidth 10 --colorMap Reds
