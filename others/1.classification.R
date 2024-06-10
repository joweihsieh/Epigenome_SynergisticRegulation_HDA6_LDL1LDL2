#/work1/home/yenmr/genome/Arabidopsis_thaliana/Ensembl/TAIR10/Annotation/Genes/Arabidosis_TE_families.txt
TE_family = read.table("Arabidosis_TE_families.txt", header = T)
TE_family[is.na(TE_family$Transposon_Super_Family),'Transposon_Super_Family'] = 'Unassigned'

#/work1/home/yenmr/project/HDAC6/RNA/raw/DETE_L12a.txt
DETE_L12a=read.table("DETE_L12a.txt", header = T, row.names = 'Row.names')
colnames(DETE_L12a)[22] = 'DETE_L12a'

TE_family2 = merge(TE_family,DETE_L12a[22], by.x = 'gene_id', by.y = 'row.names')

#/work1/home/joweihsieh/20190731_HDA6/Analysis/metaplot/Synergistic_DETE.bed
syn = read.table("Synergistic_DETE.bed",header=F)
colnames(syn)[4] = 'gene_id'


TE_family2_non_syn=TE_family2[!TE_family2$gene_id%in%syn$gene_id,]
TE_family2_syn=TE_family2[TE_family2$gene_id%in%syn$gene_id,]


matrix_L12a = as.data.frame.matrix(table(TE_family2$Transposon_Super_Family,TE_family2$DETE_L12a))
matrix_non_syn = as.data.frame.matrix(table(TE_family2_non_syn$Transposon_Super_Family,TE_family2_non_syn$DETE_L12a))
matrix_syn = as.data.frame.matrix(table(TE_family2_syn$Transposon_Super_Family,TE_family2_syn$DETE_L12a))

Enrich = function(matrix) { log2(apply(matrix,2, function(x) {x /sum(x)})/(rowSums(matrix)/sum(matrix)))}

colnames(matrix_L12a) = c('Down_L12a','No_L12a','Up_L12a')
colnames(matrix_non_syn) = c('Down_non_syn','No_non_syn','Up_non_syn')
colnames(matrix_syn) = c('Down_syn','No_syn','Up_syn')

Enrich_L12a = Enrich(matrix_L12a)
Enrich_syn = Enrich(matrix_syn)
Enrich_non_syn = Enrich(matrix_non_syn)

Enrich_L12a_1 = Enrich_L12a[c(2:5,7,9:12),]
Enrich_syn_1 = Enrich_syn[c(2:5,7,9:12),]
Enrich_non_syn_1 = Enrich_non_syn[c(2:5,7,9:12),]

Up_count = cbind(matrix_L12a['Up_L12a'], matrix_syn['Up_syn'],matrix_non_syn['Up_non_syn'])
Down_count = cbind(matrix_L12a['Down_L12a'], matrix_syn['Down_syn'],matrix_non_syn['Down_non_syn'])


Up_count = Up_count[c(2:5,7,9:12),]
Down_count = Down_count[c(2:5,7,9:12),]

Up_count_percentage = apply(Up_count, 2, function(x){x*100/sum(x, na.rm = T)})
Down_count_percentage = apply(Down_count, 2, function(x){x*100/sum(x ,na.rm = T)})







