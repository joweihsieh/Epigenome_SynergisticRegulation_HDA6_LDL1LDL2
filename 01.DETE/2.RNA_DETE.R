library(edgeR)

count = read.table("tair10_TE_count.txt",row.names="V4")
count$V5 = count$V3 - count$V2
count = count[,c(1:5,9:17)]
colnames(count)= c("Chromosome","Start","End","lengths","Direction",'C_6','C_7','C_7_2','L12_6','L12_7','L12a_6','L12a_7','a_6','a_7')

group <- c(rep("C", 3), rep("L12", 2), rep("L12a",2), rep("a",2))
y <- DGEList(count[6:14], group=group, genes = count['lengths'])

y <- calcNormFactors(y)
y <- estimateCommonDisp(y, verbose=TRUE)
y <- estimateTagwiseDisp(y)
keep <- rowSums(y$pseudo.counts[,1:9] > 10) > 1
y = y[keep,]


c = data.frame(rpkm(y))
c$C_mean = apply(c[1:3],1,mean)
c$L12_mean = apply(c[4:5],1,mean)
c$L12a_mean = apply(c[6:7],1,mean)
c$a_mean = apply(c[8:9],1,mean)


c2 = merge(count[1:5],c,by ='row.names')
write.table(c2,file = "rpkm_table_split.txt", row.names=F,quote=F, sep="\t")



et_ctrl <- exactTest(y, pair=c("C","L12"))
gene_ctrl  = topTags(et_ctrl,n = 40000000)$table
gene_ctrl$DETE = ifelse(gene_ctrl$FDR < 0.05 & gene_ctrl$logFC >= 1, 'Up',ifelse(gene_ctrl$FDR < 0.05 & gene_ctrl$logFC <= -1,'Down','No'))
ctrl = merge(c2,gene_ctrl[,c(2,4,5,6)],by.x="Row.names",by.y = "row.names")
write.table(ctrl,file = "DETE_L12.txt",quote=F,sep='\t',row.names=F)

et_ctrl <- exactTest(y, pair=c("C","L12a"))
gene_ctrl  = topTags(et_ctrl,n = 40000000)$table
gene_ctrl$DETE = ifelse(gene_ctrl$FDR < 0.05 & gene_ctrl$logFC >= 1, 'Up',ifelse(gene_ctrl$FDR < 0.05 & gene_ctrl$logFC <= -1,'Down','No'))
ctrl = merge(c2,gene_ctrl[,c(2,4,5,6)],by.x="Row.names",by.y = "row.names")
write.table(ctrl,file = "DETE_L12a.txt",quote=F,sep='\t',row.names=F)

et_ctrl <- exactTest(y, pair=c("C","a"))
gene_ctrl  = topTags(et_ctrl,n = 40000000)$table
gene_ctrl$DETE = ifelse(gene_ctrl$FDR < 0.05 & gene_ctrl$logFC >= 1, 'Up',ifelse(gene_ctrl$FDR < 0.05 & gene_ctrl$logFC <= -1,'Down','No'))
ctrl = merge(c2,gene_ctrl[,c(2,4,5,6)],by.x="Row.names",by.y = "row.names")
write.table(ctrl,file = "DETE_a.txt",quote=F,sep='\t',row.names=F)



