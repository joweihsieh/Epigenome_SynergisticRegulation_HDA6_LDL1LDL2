
plotPCA <- function(name) {

	wt1 <- read.table(paste0("bw_WGBS/At-C1.bam.",name,".tab"))
	wt2 <- read.table(paste0("bw_WGBS/At-C2.bam.",name,".tab"))
	hda1 <- read.table(paste0("bw_WGBS/At-a1.bam.",name,".tab"))
	hda2 <- read.table(paste0("bw_WGBS/At-a2.bam.",name,".tab"))
	ldl1 <- read.table(paste0("bw_WGBS/At-121.bam.",name,".tab"))
	ldl2 <- read.table(paste0("bw_WGBS/At-122.bam.",name,".tab"))
	dm1 <- read.table(paste0("bw_WGBS/At-a121.bam.",name,".tab"))
	dm2 <- read.table(paste0("bw_WGBS/At-a122.bam.",name,".tab"))

	outfile <- paste0("pca_hda6_",name,"_TE.pdf")

	df <- data.frame(row.names = wt1$V1, wt1=wt1$V6, wt2=wt2$V6, hda1=hda1$V6, hda2=hda2$V6, ldl1=ldl1$V6, ldl2=ldl2$V6, dm1=dm1$V6, dm2=dm2$V6)

	pca <- prcomp(t(df))
	summary(pca)
	vars <- (pca$sdev)^2
	pc1_percent = vars[1] / sum(vars)
	pc2_percent = vars[2] / sum(vars)
	value = pca$x

	col_WT <- 'black'
	col_hda6 <- '#4F8FF0'
	col_ldl12 <- '#ED7232'
	col_hda6ldl12 <- '#7EAB55'

	maxi <- max(abs(value[,1:2])) * 1.1
	pdf(file = outfile, height= 4.25, width = 3.75)
	plot(value[,'PC1'], value[,'PC2'], pch = 19, xlab = paste0('PC1 (',as.integer(pc1_percent * 100),'%)'), ylab = paste0('PC2 (',as.integer(pc2_percent * 100),'%)'), xlim = c(-maxi, maxi), ylim = c(-maxi, maxi), col = c(rep(col_WT,2),rep(col_hda6,2), rep(col_ldl12,2), rep(col_hda6ldl12,2)), main = paste(name, "TE"))
	#legend("topleft", pch = c(1,4,2,3), c('WT','hda6','ldl1/2','hda6/ldl1/2'),bty='n')
	dev.off()
}
col_WT <- 'black'
col_hda6 <- '#4F8FF0'
col_ldl12 <- '#ED7232'
col_hda6ldl12 <- '#7EAB55'
pdf(file = 'pca_hda6_CG_legend.pdf',height= 4, width = 3.75)
plot(NA,xlim = c(-80,100),ylim = c(-100,80),axes=F,xlab = '',ylab = '')
legend("topleft", pch = 19, col = c(col_WT, col_hda6, col_ldl12, col_hda6ldl12), c('WT','hda6','ldl1/2','hda6/ldl1/2'),bty='n')
dev.off()

plotPCA('CG')
plotPCA('CHG')
plotPCA('CHH')


