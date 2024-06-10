library(pheatmap)
library(RColorBrewer)


DETE_L12 = read.table("DETE_L12.txt",header=T, row.names = 'Row.names')
DETE_L12a = read.table("DETE_L12a.txt",header=T, row.names = 'Row.names')
DETE_a = read.table("DETE_a.txt",header=T, row.names = 'Row.names')


colnames(DETE_L12)[22] = 'DETE_L12'
colnames(DETE_L12a)[22] = 'DETE_L12a'
colnames(DETE_a)[22] = 'DETE_a'

logRPKM = log2(DETE_L12[6:14] + 1)

pca <- prcomp(t(logRPKM))
summary(pca)
vars <- (pca$sdev)^2
pc1_percent = vars[1] / sum(vars)
pc2_percent = vars[2] / sum(vars)
value = pca$x

col_WT <- 'black'
col_hda6 <- '#4F8FF0'
col_ldl12 <- '#ED7232'
col_hda6ldl12 <- '#7EAB55'

pdf(file = 'pca_hda6_RNA.pdf',height= 4.24, width = 3.75)
#plot(value[,'PC1'],value[,'PC2'], pch = c(rep(1,3),rep(2,2),rep(3,2),rep(4,2)),xlab = paste0('PC1 (',as.integer(pc1_percent * 100),'%)'),ylab = paste0('PC2 (',as.integer(pc2_percent * 100),'%)'),xlim = c(-80,100),ylim = c(-50,50))
plot(value[,'PC1'],value[,'PC2'], pch = 19,xlab = paste0('PC1 (',as.integer(pc1_percent * 100),'%)'),ylab = paste0('PC2 (',as.integer(pc2_percent * 100),'%)'),xlim = c(-80,100),ylim = c(-50,50), col = c(rep(col_WT,3),rep(col_hda6,2),rep(col_ldl12,2),rep(col_hda6ldl12,2)))
#legend("topleft", pch = c(1,4,2,3), c('WT','hda6','ldl1/2','hda6/ldl1/2'),bty='n')
dev.off()

pdf(file = 'pca_hda6_RNA_legend.pdf',height= 4, width = 3.75)
plot(NA,xlim = c(-80,100),ylim = c(-100,80),axes=F,xlab = '',ylab = '')
legend("topleft", pch = 19, col = c(col_WT, col_hda6, col_ldl12, col_hda6ldl12), c('WT','hda6','ldl1/2','hda6/ldl1/2'),bty='n')
dev.off()

pdf(file = 'pca_hda6_RNA_zoom.pdf',height= 3, width = 3)
plot(value[,'PC1'],value[,'PC2'], pch = 19, col = c(rep(col_WT,3),rep(col_hda6,2),rep(col_ldl12,2),rep(col_hda6ldl12,2)),xlab = '',ylab = '',xlim = c(-74,-65),ylim = c(-1,4))
dev.off()


#############
# TE subfamilies assignment
#############

TE_family = read.table("Arabidosis_TE_families.txt",header=T)
TE_family[is.na(TE_family$Transposon_Super_Family),'Transposon_Super_Family'] = 'Unassigned'

TE_family2 = merge(TE_family,DETE_a[22],by.x = 'gene_id',by.y = 'row.names')
TE_family3 = merge(TE_family2,DETE_L12[22],by.x = 'gene_id',by.y = 'row.names')
TE_family4 = merge(TE_family3,DETE_L12a[22],by.x = 'gene_id',by.y = 'row.names')
Exprssed_TE = merge(DETE_L12[1:18],TE_family4[c(1,5,6,7,4)],by.x = 'row.names',by.y = 'gene_id')
write.table(Exprssed_TE,file = 'Exprssed_TE.txt',row.names=F, quote=F, sep='\t')


#############
# Enrichment of TE subfamilies
#############

matrix_L12 = as.data.frame.matrix(table(TE_family4$Transposon_Super_Family,TE_family4$DETE_L12))
matrix_L12a = as.data.frame.matrix(table(TE_family4$Transposon_Super_Family,TE_family4$DETE_L12a))
matrix_a = as.data.frame.matrix(table(TE_family4$Transposon_Super_Family,TE_family4$DETE_a))
Enrich = function(matrix) { log2(apply(matrix,2, function(x) {x /sum(x)})/(rowSums(matrix)/sum(matrix)))}
colnames(matrix_L12) = c('Down_L12','No_L12','Up_L12')
colnames(matrix_L12a) = c('Down_L12a','No_L12a','Up_L12a')
colnames(matrix_a) = c('Down_a','No_a','Up_a')
Enrich_L12 = Enrich(matrix_L12)
Enrich_L12a = Enrich(matrix_L12a)
Enrich_a = Enrich(matrix_a)



Up_count = cbind(matrix_a['Up_a'],matrix_L12['Up_L12'],matrix_L12a['Up_L12a'])
Down_count = cbind(matrix_a['Down_a'],matrix_L12['Down_L12'],matrix_L12a['Down_L12a'])
Up_count = Up_count[c(2:5,7,9:12),]
Down_count = Down_count[c(2:5,7,9:12),]

Up_count_percentage=apply(Up_count, 2, function(x){x*100/sum(x,na.rm=T)})
Down_count_percentage=apply(Down_count, 2, function(x){x*100/sum(x,na.rm=T)})


coul = brewer.pal(9, "Set1") 

png('Up_count_percentage.png',height= 2000, width = 2000, res=400)
barplot(as.matrix(Up_count_percentage[nrow(Up_count_percentage):1,]),beside = F,col= rev(coul),ylab = 'Percentage')
dev.off()

png('Up_count.png',height= 2000, width = 2000, res=400)
barplot(as.matrix(Up_count[nrow(Up_count):1,]),beside = F,col= rev(coul),ylab = 'Count')
dev.off()

png('Down_count_percentage.png',height= 2000, width = 2000, res=400)
barplot(as.matrix(Down_count_percentage[nrow(Down_count_percentage):1,]),beside = F,col= rev(coul),ylab = 'Percentage')
dev.off()

png('Down_count.png',height= 2000, width = 2000, res=400)
barplot(as.matrix(Down_count[nrow(Down_count):1,]),beside = F,col= rev(coul),ylab = 'Count')
dev.off()

png('Legend.png',height= 1000, width = 700, res=400)
par(mar = c(0,0,0,0))
plot(NA,ylim =c(0,10),xlim = c(0,10),axes=F,ylab = '',xlab = '')
legend("topright", legend = row.names(Up_count_percentage), fill = coul,bty = "n")
dev.off()

Up_enrich = cbind(Enrich_a[,'Up_a',drop=F],Enrich_L12[,'Up_L12',drop=F],Enrich_L12a[,'Up_L12a',drop=F])
Up_enrich = Up_enrich[c(2:5,7,9:12),]
Down_enrich = cbind(Enrich_a[,'Down_a',drop=F],Enrich_L12[,'Down_L12',drop=F],Enrich_L12a[,'Down_L12a',drop=F])
Down_enrich = Down_enrich[c(2:5,7,9:12),]

png('Down_enrich.png',height= 1800, width = 2000, res=400)
barplot(Down_enrich,beside=T,ylim = c(-3,3),col= coul,ylab = 'Log2 fold enrichment')
abline(v = c(10.5,20.5),lty = 2)
dev.off()

png('Up_enrich.png',height= 1800, width = 2000, res=400)
barplot(Up_enrich,beside=T,ylim = c(-3,3),col= coul,ylab = 'Log2 fold enrichment')
abline(v = c(10.5,20.5),lty = 2)
dev.off()

