#computeMatrix scale-regions -S At-C1.bam.CHH.bw At-C2.bam.CHH.bw -R tair10_TE.bed -p 20 -bs 100 -a 2000 -b 2000 -m 2000 -out TE_CHH_WT.matrix.gz
#computeMatrix scale-regions -S At-a1.bam.CHH.bw At-a2.bam.CHH.bw -R tair10_TE.bed -p 20 -bs 100 -a 2000 -b 2000 -m 2000 -out TE_CHH_hda6.matrix.gz
#computeMatrix scale-regions -S At-121.bam.CHH.bw At-122.bam.CHH.bw -R tair10_TE.bed -p 20 -bs 100 -a 2000 -b 2000 -m 2000 -out TE_CHH_ldl12.matrix.gz
#computeMatrix scale-regions -S At-a121.bam.CHH.bw At-a122.bam.CHH.bw -R tair10_TE.bed -p 20 -bs 100 -a 2000 -b 2000 -m 2000 -out TE_CHH_hda6ldl12.matrix.gz

#perl removeNaMatrix.pl TE_CHH_WT.matrix.gz |gzip > TE_CHH_WTs.matrix.gz
#perl removeNaMatrix.pl TE_CHH_hda6.matrix.gz |gzip > TE_CHH_hda6s.matrix.gz
#perl removeNaMatrix.pl TE_CHH_ldl12.matrix.gz |gzip > TE_CHH_ldl12s.matrix.gz
#perl removeNaMatrix.pl TE_CHH_hda6ldl12.matrix.gz |gzip > TE_CHH_hda6ldl12s.matrix.gz

perl removeNaMatrixIntersect.pl TE_CHH_WT.matrix.gz |gzip > TE_CHH_WTi.matrix.gz
perl removeNaMatrixIntersect.pl TE_CHH_hda6.matrix.gz |gzip > TE_CHH_hda6i.matrix.gz
perl removeNaMatrixIntersect.pl TE_CHH_ldl12.matrix.gz |gzip > TE_CHH_ldl12i.matrix.gz
perl removeNaMatrixIntersect.pl TE_CHH_hda6ldl12.matrix.gz |gzip > TE_CHH_hda6ldl12i.matrix.gz


plotHeatmap -m TE_CHH_WTi.matrix.gz -out TE_CHH_WTi.png --kmeans 3 --heatmapWidth 10 --outFileSortedRegions TE_CHH_WTi_sort.txt --zMax 0.18 --colorMap Reds --yMax 0.12 --yMin 0.0
plotHeatmap -m TE_CHH_hda6i.matrix.gz -out TE_CHH_hda6i.png --kmeans 3 --heatmapWidth 10 --outFileSortedRegions TE_CHH_hda6i_sort.txt --zMax 0.18 --colorMap Reds --yMax 0.12 --yMin 0.0
plotHeatmap -m TE_CHH_ldl12i.matrix.gz -out TE_CHH_ldl12i.png --kmeans 3 --heatmapWidth 10 --outFileSortedRegions TE_CHH_ldl12i_sort.txt --zMax 0.18 --colorMap Reds --yMax 0.12 --yMin 0.0
plotHeatmap -m TE_CHH_hda6ldl12i.matrix.gz -out TE_CHH_hda6ldl12i.png --kmeans 3 --heatmapWidth 10 --outFileSortedRegions TE_CHH_hda6ldl12i_sort.txt --zMax 0.18 --colorMap Reds --yMax 0.12 --yMin 0.0

#computeMatrix scale-regions -S At-C1.bam.CHH.bw At-C2.bam.CHH.bw At-a1.bam.CHH.bw At-a2.bam.CHH.bw At-121.bam.CHH.bw At-122.bam.CHH.bw At-a121.bam.CHH.bw At-a122.bam.CHH.bw -R tair10_TE.bed -p 20 -bs 100 -a 2000 -b 2000 -m 2000 -out TE_CHH_all.matrix.gz
#plotHeatmap -m TE_CHH_all.matrix.gz -out TE_CHH_all.png --kmeans 3 --heatmapWidth 10

plotHeatmap -m TE_CHH_WTi.matrix.gz -out TE_CHH_WTi_oc.png --kmeans 3 --heatmapWidth 10 --outFileSortedRegions TE_CHH_WTi_sort.txt --zMax 0.18 --yMax 0.12 --yMin 0.0
plotHeatmap -m TE_CHH_hda6i.matrix.gz -out TE_CHH_hda6i_oc.png --kmeans 3 --heatmapWidth 10 --outFileSortedRegions TE_CHH_hda6i_sort.txt --zMax 0.18 --yMax 0.12 --yMin 0.0
plotHeatmap -m TE_CHH_ldl12i.matrix.gz -out TE_CHH_ldl12i_oc.png --kmeans 3 --heatmapWidth 10 --outFileSortedRegions TE_CHH_ldl12i_sort.txt --zMax 0.18 --yMax 0.12 --yMin 0.0
plotHeatmap -m TE_CHH_hda6ldl12i.matrix.gz -out TE_CHH_hda6ldl12i_oc.png --kmeans 3 --heatmapWidth 10 --outFileSortedRegions TE_CHH_hda6ldl12i_sort.txt --zMax 0.18 --yMax 0.12 --yMin 0.0




