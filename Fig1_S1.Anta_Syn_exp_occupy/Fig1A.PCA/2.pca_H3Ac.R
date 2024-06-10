wt1 <- read.table("bw_ChIP/ChIP_H3Ac_WT_1.tab")
wt2 <- read.table("bw_ChIP/ChIP_H3Ac_WT_2.tab")
hda1 <- read.table("bw_ChIP/ChIP_H3Ac_hda6_1.tab")
hda2 <- read.table("bw_ChIP/ChIP_H3Ac_hda6_2.tab")
ldl1 <- read.table("bw_ChIP/ChIP_H3Ac_ldl12_1.tab")
ldl2 <- read.table("bw_ChIP/ChIP_H3Ac_ldl12_2.tab")
dm1 <- read.table("bw_ChIP/ChIP_H3Ac_hda6ldl12_1.tab")
dm2 <- read.table("bw_ChIP/ChIP_H3Ac_hda6ldl12_2.tab")
outfile <- 'pca_hda6_H3Ac_TE.pdf'

df <- data.frame(row.names = wt1$V1, wt1=wt1$V6, wt2=wt2$V6, hda1=hda1$V6, hda2=hda2$V6, ldl1=ldl1$V6, ldl2=ldl2$V6, dm1=dm1$V6, dm2=dm2$V6)
df_log <- log2(df+1)

pca <- prcomp(t(df_log))
summary(pca)
vars <- (pca$sdev)^2
pc1_percent = vars[1] / sum(vars)
pc2_percent = vars[2] / sum(vars)
value = pca$x

col_WT <- 'black'
col_hda6 <- '#4F8FF0'
col_ldl12 <- '#ED7232'
col_hda6ldl12 <- '#7EAB55'

maxi <- max(abs(value[,1:2])) + 5
pdf(file = outfile, height= 4.25, width = 3.75)
plot(value[,'PC1'], value[,'PC2'], pch = 19, xlab = paste0('PC1 (',as.integer(pc1_percent * 100),'%)'), ylab = paste0('PC2 (',as.integer(pc2_percent * 100),'%)'), xlim = c(-maxi, maxi), ylim = c(-maxi, maxi), col = c(rep(col_WT,2),rep(col_hda6,2), rep(col_ldl12,2), rep(col_hda6ldl12,2)), main = "H3Ac TE")
#legend("topleft", pch = c(1,4,2,3), c('WT','hda6','ldl1/2','hda6/ldl1/2'),bty='n')
dev.off()

