fpkm_1 <- read.table("fpkm_GeneNTE.txt")
fpkm_2 <- read.table("fpkm_TEgene.txt",header=T, row.names="Row.names")
fpkm_3 <- read.table("fpkm_TE.txt",header=T, row.names="Row.names")
fpkm_3$WT <- apply(fpkm_3[,6:8],1,mean)
fpkm_3$hda6 <- apply(fpkm_3[,13:14],1,mean)
fpkm_3$ldl12 <- apply(fpkm_3[,9:10],1,mean)
fpkm_3$hda6ldl12 <- apply(fpkm_3[,11:12],1,mean)

TEgene <- read.table("tair10_TE.bed")
TEelement <- read.table("tair10_Transposon_element_sorted.bed")
gene <- read.table("tair10_gene2.bed")
expressedTEgene <- read.table("/work4/home/yenmr/project/HDAC6/RNA/raw/DETE_a.txt",header=T,  row.names="Row.names")
expressedGene <- read.table("~/project/HDAC6/RNA/raw/DEG_hda6.txt",header=T)
fpkm_gene <- fpkm_1[row.names(fpkm_1) %in% row.names(expressedGene),]
fpkm_gene$WT_mean <- apply(fpkm_gene[,1:3],1,mean)
fpkm_gene$hda6_mean <- apply(fpkm_gene[,8:9],1,mean)
fpkm_gene$ldl12_mean <- apply(fpkm_gene[,4:5],1,mean)
fpkm_gene$hda6ldl12_mean <- apply(fpkm_gene[,6:7],1,mean)

fpkm_gene <- log2(fpkm_gene[,10:13] + 1)
names(fpkm_gene) <- c('WT','hda6','ldl12','hda6ldl12')
pdf(file = "Exp_boxplot_gene.pdf", width=6, height=6)
boxplot(fpkm_gene$WT, fpkm_gene$hda6, fpkm_gene$ldl12, fpkm_gene$hda6ldl12, main="Expression of genes", names=names(fpkm_gene),ylab = "Log2(FPKM + 1)",ylim = c(0,18))
dev.off()
apply(fpkm_gene,2,mean)
apply(fpkm_gene,2,median)
t.test(fpkm_gene$WT, fpkm_gene$hda6,paired=T)$p.value
t.test(fpkm_gene$WT, fpkm_gene$ldl12,paired=T)$p.value
t.test(fpkm_gene$WT, fpkm_gene$hda6ldl12,paired=T)$p.value

fpkm_TEgene <- log2(expressedTEgene[,c(15,18,16,17)] + 1)
names(fpkm_TEgene) <- c('WT','hda6','ldl12','hda6ldl12')
pdf(file = "Exp_boxplot_TEgene.pdf", width=6, height=6)
boxplot(fpkm_TEgene$WT, fpkm_TEgene$hda6, fpkm_TEgene$ldl12, fpkm_TEgene$hda6ldl12, main="Expression of TE genes", names=names(fpkm_TEgene),ylab = "Log2(FPKM + 1)",ylim = c(0,18))
dev.off()
apply(fpkm_TEgene,2,mean)
apply(fpkm_TEgene,2,median)
t.test(fpkm_TEgene$WT, fpkm_TEgene$hda6,paired=T)$p.value
t.test(fpkm_TEgene$WT, fpkm_TEgene$ldl12,paired=T)$p.value
t.test(fpkm_TEgene$WT, fpkm_TEgene$hda6ldl12,paired=T)$p.value



fpkm_TEelement <- log2(fpkm_3[row.names(fpkm_3) %in% TEelement$V4,c(15,18,16,17)]+ 1)
names(fpkm_TEelement) <- c('WT','hda6','ldl12','hda6ldl12')
pdf(file = "Exp_boxplot_TEelement.pdf", width=6, height=6)
boxplot(fpkm_TEelement$WT, fpkm_TEelement$hda6, fpkm_TEelement$ldl12, fpkm_TEelement$hda6ldl12, main="Expression of TE elements", names=names(fpkm_TEelement),ylab = "Log2(FPKM + 1)",ylim = c(0,18))
dev.off()
apply(fpkm_TEelement,2,mean)
apply(fpkm_TEelement,2,median)
t.test(fpkm_TEelement$WT, fpkm_TEelement$hda6,paired=T)$p.value
t.test(fpkm_TEelement$WT, fpkm_TEelement$ldl12,paired=T)$p.value
t.test(fpkm_TEelement$WT, fpkm_TEelement$hda6ldl12,paired=T)$p.value



