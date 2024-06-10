library(edgeR)

c_6 = read.table("C_6.txt")
c_7 = read.table("C_7.txt")
c_7_2 = read.table("C_7_2.txt")
L12_6 = read.table("L12_6.txt")
L12_7 = read.table("L12_7.txt")
L12a_6 = read.table("L12a_6.txt")
L12a_7 = read.table("L12a_7.txt")
a_6 = read.table("a_6.txt")
a_7 = read.table("a_7.txt")

TE = read.table("tair10_TE.bed")
gene = read.table("tair10_gene2.bed")

gene_size = read.table("gene_size.txt",header=F,row.names = 'V1')
colnames(gene_size) = "lengths"

count = cbind.data.frame(row.names= c_6$V1, c_6 = c_6$V2, c_7=c_7$V2, c_7_2=c_7_2$V2, L12_6=L12_6$V2, L12_7=L12_7$V2, L12a_6=L12a_6$V2, L12a_7=L12a_7$V2, a_6=a_6$V2, a_7=a_7$V2)
x = count[row.names(count) %in% row.names(gene_size),]						 

#x = count[1:(nrow(count) -5),]
#write.table(x,file='combined_count2.txt',quote=F, sep="\t")
#ngene_size = gene_size[order(match(row.names(x),ngene_size$V1)),]
m = match(row.names(x),row.names(gene_size))
ngene_size = gene_size[m,,drop=F]
#group <- c("col","axe","jmj28","jmj28axe")
group <- c(rep("c", 3) , rep("L12", 2), rep("L12a",2),rep("a",2))
y <- DGEList(x, group=group, genes = ngene_size)
keep = rowSums(cpm(y) > 0) >= 1

y = y[keep,]
#y$samples$lib.size <- colSums(y$counts)
y <- calcNormFactors(y)
#y$common.dispersion = 0.2
#y$common.dispersion = 0.1157891
y <- estimateCommonDisp(y, verbose=TRUE)
y <- estimateTagwiseDisp(y)
keep <- rowSums(y$pseudo.counts > 10) > 1
y = y[keep,]


c = data.frame(rpkm(y))
#write.table(c, file = "rpkm_2.txt", quote=F, sep="\t")


et_ctrl <- exactTest(y, pair=c("c","a"))
gene_ctrl  = topTags(et_ctrl,n = 40000000)$table
gene_ctrl$DEG = ifelse(gene_ctrl$FDR < 0.05 & gene_ctrl$logFC >= 1, 'Up',ifelse(gene_ctrl$FDR < 0.05 & gene_ctrl$logFC <= -1,'Down','No'))
gene_ctrl = gene_ctrl[row.names(gene_ctrl) %in% gene$V4,]
table(gene_ctrl$DEG)
write.table(gene_ctrl,file = "DEG_hda6.txt",quote=F,sep='\t',row.names=T)

et_ctrl <- exactTest(y, pair=c("c","L12"))
gene_ctrl  = topTags(et_ctrl,n = 40000000)$table
gene_ctrl$DEG = ifelse(gene_ctrl$FDR < 0.05 & gene_ctrl$logFC >= 1, 'Up',ifelse(gene_ctrl$FDR < 0.05 & gene_ctrl$logFC <= -1,'Down','No'))
gene_ctrl = gene_ctrl[row.names(gene_ctrl) %in% gene$V4,]
table(gene_ctrl$DEG)


write.table(gene_ctrl,file = "DEG_ldl12.txt",quote=F,sep='\t',row.names=T)

et_ctrl <- exactTest(y, pair=c("c","L12a"))
gene_ctrl  = topTags(et_ctrl,n = 40000000)$table
gene_ctrl$DEG = ifelse(gene_ctrl$FDR < 0.05 & gene_ctrl$logFC >= 1, 'Up',ifelse(gene_ctrl$FDR < 0.05 & gene_ctrl$logFC <= -1,'Down','No'))
gene_ctrl = gene_ctrl[row.names(gene_ctrl) %in% gene$V4,]
table(gene_ctrl$DEG)

write.table(gene_ctrl,file = "DEG_hda6ldl12.txt",quote=F,sep='\t',row.names=T)

