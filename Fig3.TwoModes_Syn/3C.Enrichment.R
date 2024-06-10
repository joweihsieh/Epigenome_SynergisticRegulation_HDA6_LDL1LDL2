#/work1/home/yenmr/genome/Arabidopsis_thaliana/Ensembl/TAIR10/Annotation/Genes/Arabidosis_TE_families.txt
TE_family = read.table("Arabidosis_TE_families.txt",header=T)
TE_family[is.na(TE_family$Transposon_Super_Family),'Transposon_Super_Family'] = 'Unassigned'

#/work1/home/joweihsieh/20190731_HDA6/Analysis/heatmap/upTE_union_value_binding_H3Ac_H3K4me2_WGBS_type_original_101_001_20220802_ATACseq.txt
data=read.table("upTE_union_value_binding_H3Ac_H3K4me2_WGBS_type_original_101_001_20220802_ATACseq.txt",sep="\t",header=T)

#class

Syn_101=data[data$types=="1_0_1",] # HDA6-predominant
Syn_001=data[data$types=="0_0_1",] # double-locked
Anta=data[data$types=="strong_weak_anta",]


#### expressed TEs
#~/20190731_HDA6/Analysis/heatmap/Exprssed_TE.txt
Data3=read.table("Exprssed_TE.txt",header=T)



TE_family_Syn_101=TE_family[TE_family$gene_id%in%Syn_101$Row.names,]
TE_family_Syn_101$type="Syn_101"
TE_family_Syn_001=TE_family[TE_family$gene_id%in%Syn_001$Row.names,]
TE_family_Syn_001$type="Syn_001"
TE_family_Anta=TE_family[TE_family$gene_id%in%Anta$Row.names,]
TE_family_Anta$type="Anta"




TE_family_expressed_TE=TE_family[TE_family$gene_id%in%Data3$Row.names,]
TE_family_expressed_TE$type="expressed_TE"


TE_family_mydata=rbind(TE_family_Syn_101,TE_family_Syn_001,TE_family_Anta, TE_family_expressed_TE)
matrix = as.data.frame.matrix(table(TE_family_mydata$Transposon_Super_Family,TE_family_mydata$type))


number=matrix$expressed_TE
#matrix$number=number

total=sum(matrix$expressed_TE)

#Enrich = function(matrix) { log2(apply(matrix,2, function(x) {x /sum(x)})/(rowSums(matrix)/sum(matrix)))}
#Enrich = function(matrix) { log2(apply(matrix,2, function(x) {x /sum(x)})/number/total)}


matrix$Syn_101_enrich=log((matrix[1]/sum(matrix[1]))/(matrix$expressed_TE/total),2)
matrix$Syn_001_enrich=log((matrix[2]/sum(matrix[2]))/(matrix$expressed_TE/total),2)
matrix$Anta_enrich=log((matrix[3]/sum(matrix[3]))/(matrix$expressed_TE/total),2)



write.table(matrix,"enrichment_anta_syn_001_101_20230607.txt",sep="\t",quote=F,row.names=F,col.names=F)

