# TE family: 
# Enrichment code: /work1/home/joweihsieh/20190731_HDA6/Analysis/enrichment/1.classification.R
# cluster: /work1/home/joweihsieh/20190731_HDA6/Analysis/CMT2_RdDM/upTE_union_value_binding_H3Ac_H3K4me2_WGBS_type_original_101_001_cluster_20220804.txt

###########################################
#/work1/home/yenmr/genome/Arabidopsis_thaliana/Ensembl/TAIR10/Annotation/Genes/Arabidosis_TE_families.txt
TE_family = read.table("Arabidosis_TE_families.txt",header=T)
TE_family[is.na(TE_family$Transposon_Super_Family),'Transposon_Super_Family'] = 'Unassigned'

#/work1/home/joweihsieh/20190731_HDA6/Analysis/CMT2_RdDM/TE_CHH_WTi_sort.txt
data = read.table("TE_CHH_WTi_sort.txt",sep="\t",header=F)


WT_cluster_1 =  data[data$V13=="cluster_1",]
WT_cluster_2 =  data[data$V13=="cluster_2",]
WT_cluster_3 =  data[data$V13=="cluster_3",]

TE_family_WT_cluster_1 = TE_family[TE_family$gene_id%in%WT_cluster_1$V4,]
TE_family_WT_cluster_1$type = "WT_cluster_1"

TE_family_WT_cluster_2 = TE_family[TE_family$gene_id%in%WT_cluster_2$V4,]
TE_family_WT_cluster_2$type = "WT_cluster_2"

TE_family_WT_cluster_3 = TE_family[TE_family$gene_id%in%WT_cluster_3$V4,]
TE_family_WT_cluster_3$type = "WT_cluster_3"




TE_family_mydata = rbind(TE_family_WT_cluster_1,TE_family_WT_cluster_2,TE_family_WT_cluster_3)
matrix = as.data.frame.matrix(table(TE_family_mydata$Transposon_Super_Family,TE_family_mydata$type))


matrix$submafily_sum = apply(matrix, 1, function(x) sum(x))

total = sum(matrix$submafily_sum)


matrix$WT_cluster_1_enriched = log((matrix[1]/sum(matrix[1]))/(matrix$submafily_sum/total),2)
matrix$WT_cluster_2_enriched = log((matrix[2]/sum(matrix[2]))/(matrix$submafily_sum/total),2)

matrix$WT_cluster_3_enriched = log((matrix[3]/sum(matrix[3]))/(matrix$submafily_sum/total),2)

mydata = data.frame(matrix)
write.table(matrix, "cluster_TE_subfamily_enrichment.txt", sep="\t", quote=F)