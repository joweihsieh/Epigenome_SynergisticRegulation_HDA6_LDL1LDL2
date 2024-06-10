computeMatrix scale-regions -S ad-kc4.perfect.24nt.siRNA.merged.bw ad-kc5.perfect.24nt.siRNA.merged.bw ad-kc6.perfect.24nt.siRNA.merged.bw -R TE_CHH_WT_cluster1.bed TE_CHH_WT_cluster2.bed TE_CHH_WT_cluster3.bed -p 20 -bs 100 -a 2000 -b 2000 -m 2000 -out TE_si24_WT.matrix.gz --missingDataAsZero

computeMatrix scale-regions -S ad-ka45.perfect.24nt.siRNA.merged.bw ad-ka46.perfect.24nt.siRNA.merged.bw ad-ka7.perfect.24nt.siRNA.merged.bw -R TE_CHH_hda6_cluster1.bed TE_CHH_hda6_cluster2.bed TE_CHH_hda6_cluster3.bed -p 20 -bs 100 -a 2000 -b 2000 -m 2000 -out TE_si24_hda6.matrix.gz --missingDataAsZero

computeMatrix scale-regions -S ad-k1245.perfect.24nt.siRNA.merged.bw ad-k1246.perfect.24nt.siRNA.merged.bw ad-k127.perfect.24nt.siRNA.merged.bw -R TE_CHH_ldl12_cluster1.bed TE_CHH_ldl12_cluster2.bed TE_CHH_ldl12_cluster3.bed -p 20 -bs 100 -a 2000 -b 2000 -m 2000 -out TE_si24_ldl12.matrix.gz --missingDataAsZero

computeMatrix scale-regions -S ad-k12a4.perfect.24nt.siRNA.merged.bw ad-k12a5.perfect.24nt.siRNA.merged.bw ad-k12a6.perfect.24nt.siRNA.merged.bw -R TE_CHH_hda6ldl12_cluster1.bed TE_CHH_hda6ldl12_cluster2.bed TE_CHH_hda6ldl12_cluster3.bed -p 20 -bs 100 -a 2000 -b 2000 -m 2000 -out TE_si24_hda6ldl12.matrix.gz --missingDataAsZero

perl setMax.pl TE_si24_WT.matrix.gz|gzip > TE_si24_WTl.matrix.gz
perl setMax.pl TE_si24_hda6.matrix.gz|gzip > TE_si24_hda6l.matrix.gz
perl setMax.pl TE_si24_ldl12.matrix.gz|gzip > TE_si24_ldl12l.matrix.gz
perl setMax.pl TE_si24_hda6ldl12.matrix.gz|gzip > TE_si24_hda6ldl12l.matrix.gz

plotHeatmap -m TE_si24_WTl.matrix.gz -out TE_si24_WT.png --heatmapWidth 10 --sortRegions no --zMax 4
plotHeatmap -m TE_si24_hda6l.matrix.gz -out TE_si24_hda6.png --heatmapWidth 10 --sortRegions no --zMax 4
plotHeatmap -m TE_si24_ldl12l.matrix.gz -out TE_si24_ldl12.png --heatmapWidth 10 --sortRegions no --zMax 4
plotHeatmap -m TE_si24_hda6ldl12l.matrix.gz -out TE_si24_hda6ldl12.png --heatmapWidth 10 --sortRegions no --zMax 4

