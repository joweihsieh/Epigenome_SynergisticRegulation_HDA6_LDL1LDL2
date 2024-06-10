library(edgeR)


a = read.table("Count_TE.txt")
a = a[,c(1:6,10:15)]

#ChIP_hda6_1.bam ChIP_hda6_2.bam ChIP_ldl1_1.bam ChIP_ldl1_2.bam ChIP_hda6_input.bam ChIP_ldl1_input.bam -bed tair10_TE.bed > Count_TE.txt
colnames(a) = c("Chr","Str","End","Gene_ID","Size","Direction", "hda6_1", "hda6_2", "ldl1_1", "ldl1_2", "input_1", "input_2")
row.names(a) = a$Gene_ID


group <- c(rep("hda6", 2) , rep("ldl1", 2), rep("input", 2)) 

y <- DGEList(a[7:12], group=group)
#keep = rowSums(cpm(y[,1:4]) > 0) >= 1 & rowSums(cpm(y[,5:7]) > 0) >= 1
#y = y[keep,]
y <- calcNormFactors(y)
y <- estimateCommonDisp(y, verbose=TRUE)
y <- estimateTagwiseDisp(y)


c = data.frame(cpm(y))
write.table(c, file = "CPM_TE.txt", quote=F, sep="\t")

et_hda6 <- exactTest(y, pair=c("input","hda6"))
gene_hda6  = topTags(et_hda6,n = 40000000)$table
#Ggene_hda6$DHMG = ifelse(gene_hda6$FDR < 0.05 & gene_hda6$logFC >= 0.5849625, 'Yes','No')
gene_hda6$DHMG = ifelse(gene_hda6$FDR < 0.05 & gene_hda6$logFC >= 0, 'Yes','No')
hda6 = merge(c,gene_hda6,by.x="row.names",by.y = "row.names")
write.table(hda6,file = "DHMG_TE_hda6.txt",quote=F,sep='\t',row.names=F)

et_ldl1 <- exactTest(y, pair=c("input","ldl1"))
gene_ldl1  = topTags(et_ldl1,n = 40000000)$table
#gene_ldl1$DHMG = ifelse(gene_ldl1$FDR < 0.05 & gene_ldl1$logFC >= 1, 'Yes','No')
gene_ldl1$DHMG = ifelse(gene_ldl1$FDR < 0.05 & gene_ldl1$logFC >= 0, 'Yes','No')
ldl1 = merge(c,gene_ldl1,by.x="row.names",by.y = "row.names")
write.table(ldl1,file = "DHMG_TE_ldl1.txt",quote=F,sep='\t',row.names=F)


hda6_bind = hda6[hda6$DHMG =='Yes',]
ldl1_bind = ldl1[ldl1$DHMG =='Yes',]
dim(hda6_bind[hda6_bind$Row.names %in% ldl1_bind$Row.names,])

