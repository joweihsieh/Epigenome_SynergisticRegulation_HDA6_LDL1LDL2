## CMT OR DRM loci
#https://www.biorxiv.org/content/biorxiv/early/2021/07/31/2021.07.31.454434.full.pdf
#/work1/home/joweihsieh/genome/Arabidopsis/TAIR10
CMT_TE.bed
# 4486
DRM_TE.bed
# 3039
CMT_and_DRM_TE.bed
# 80


bedtools intersect -a tair10_TE.bed -b DRM_TE.bed -wa -wb >  tair10_TEgene_DRM.txt
#37
bedtools intersect -a tair10_TE.bed -b CMT_TE.bed -wa -wb >  tair10_TEgene_CMT.txt
#3857
bedtools intersect -a tair10_TE.bed -b CMT_and_DRM_TE.bed -wa -wb >  tair10_TEgene_CMT_and_DRM.txt
#4



