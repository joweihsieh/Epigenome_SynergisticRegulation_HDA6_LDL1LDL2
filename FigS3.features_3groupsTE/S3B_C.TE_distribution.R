library(ggplot2)

TE_anta <- read.table("tair10_antagonistic_TE.bed")
TE_001 <- read.table("tair10_synergistic_TE_001.bed")
TE_101 <- read.table("tair10_synergistic_TE_101.bed")
TE_gene_distance <- read.table("TE_closed_TSS.txt")
TE_gene_distance$dis1 <- abs(TE_gene_distance$V6 - TE_gene_distance$V2)
TE_gene_distance$dis2 <- abs(TE_gene_distance$V6 - TE_gene_distance$V3)
TE_gene_distance$dis2gene <- apply(TE_gene_distance[,c('dis1','dis2')], 1, min)
set_max = 50000
TE_gene_distance$dis2gene <- ifelse(TE_gene_distance$dis2gene < set_max, TE_gene_distance$dis2gene, set_max) 
peri = read.table("At-C1.bam.CG.peri.bed")

TE_anta$type <- "Anta"
TE_001$type <- '001'
TE_101$type <- '101'
TE_anta$col <- 'red'
TE_001$col <- 'blue'
TE_101$col <- 'darkgreen'

TE <- rbind(TE_anta[,c('V1','V2','V3','V4','type','col')],TE_001[,c('V1','V2','V3','V4','type','col')],TE_101[,c('V1','V2','V3','V4','type','col')])
TE$type <- factor(TE$type, levels = c('Anta','001','101'))
TE$lengths <- TE$V3 - TE$V2
TE$dis2gene <- TE_gene_distance[match(TE$V4, TE_gene_distance$V4), 'dis2gene']
p <- ggplot(TE, aes(x=type, y=lengths)) + geom_violin() + geom_boxplot(width=0.1) + ylab("Bp") + ggtitle("TE lengths distribution") + xlab("") + theme_minimal()
ggsave("TE_Violin_size_antagonistic_synergistic.pdf")

p <- ggplot(TE, aes(x=type, y=lengths)) + geom_boxplot() + ylab("Bp") + ggtitle("TE lengths distribution") + xlab("") + theme_minimal()
ggsave("TE_Boxplot_size_antagonistic_synergistic.pdf")

paste('lengths','Anta','101',t.test(TE[TE$type=='Anta','lengths'], TE[TE$type=='101','lengths'])$p.value)
paste('lengths','Anta','001',t.test(TE[TE$type=='Anta','lengths'], TE[TE$type=='001','lengths'])$p.value)
paste('lengths','001','101',t.test(TE[TE$type=='001','lengths'], TE[TE$type=='101','lengths'])$p.value)
paste('dis2gene','Anta','101',t.test(TE[TE$type=='Anta','dis2gene'], TE[TE$type=='101','dis2gene'])$p.value)
paste('dis2gene','Anta','001',t.test(TE[TE$type=='Anta','dis2gene'], TE[TE$type=='001','dis2gene'])$p.value)
paste('dis2gene','001','101',t.test(TE[TE$type=='001','dis2gene'], TE[TE$type=='101','dis2gene'])$p.value)

paste('wilcox lengths','Anta','101',wilcox.test(TE[TE$type=='Anta','lengths'], TE[TE$type=='101','lengths'])$p.value)
paste('wilcox lengths','Anta','001',wilcox.test(TE[TE$type=='Anta','lengths'], TE[TE$type=='001','lengths'])$p.value)
paste('wilcox lengths','001','101',wilcox.test(TE[TE$type=='001','lengths'], TE[TE$type=='101','lengths'])$p.value)
paste('wilcox dis2gene','Anta','101',wilcox.test(TE[TE$type=='Anta','dis2gene'], TE[TE$type=='101','dis2gene'])$p.value)
paste('wilcox dis2gene','Anta','001',wilcox.test(TE[TE$type=='Anta','dis2gene'], TE[TE$type=='001','dis2gene'])$p.value)
paste('wilcox dis2gene','001','101',wilcox.test(TE[TE$type=='001','dis2gene'], TE[TE$type=='101','dis2gene'])$p.value)

p <- ggplot(TE, aes(x=type, y=dis2gene)) + geom_violin() + ylab("Bp") + geom_boxplot(width=0.1) + ggtitle("TE to gene distance (max 50 kb)") + xlab("") + theme_minimal()
ggsave("TE_Violin_dis2gene_antagonistic_synergistic.pdf") + theme_minimal()
p <- ggplot(TE, aes(x=type, y=dis2gene)) + geom_boxplot() + ylab("Bp") + ggtitle("TE to gene distance (max 50 kb)") + xlab("") + theme_minimal()
ggsave("TE_Boxplot_dis2gene_antagonistic_synergistic.pdf") + theme_minimal()


genome_size <- c(30427671, 19698289, 23459830, 18585056, 26975502)

pdf("TE_genome_distribution.pdf", width = 6, height=4)
width <- 0.4
scale <- 1000000
#plot(NA, xlim = c(0, 35), ylim = c(0,6),xlab = "", ylab = "", yaxt = "n", axes = FALSE)
plot(NA, type = "n", xlab = "Genome position (Mb)", yaxt = "n", bty = "n", ylab = "Chromosome", xlim = c(0, 35), ylim = c(0,6))
y_base = 5

for (i in 1:nrow(TE)) {
	x <- TE[i, "V2"] / scale
	y <- TE[i, "V1"]
	col <- TE[i, "col"]
	y_add <- ifelse(TE[i, "type"]=='Anta', width, ifelse(TE[i, "type"]=='101', width * 2, width * 3 ))

	lines(x = c(x,x), y = c(y_base + width - y, y_base - y), col = col, lwd = 0.2)
	#lines(x = c(x,x), y = c(y_base + y_add - width - y, y_base + y_add - y), col = col, lwd = 0.5)
}

for (y in seq(1,5)) {
	x <- genome_size[y] / scale
	rect(xleft = 0, ybottom = y_base + width - y , xright = x, ytop = y_base - y, border = "black")
#	rect(xleft = 0, ybottom = y_base + width - y , xright = x, ytop = y_base - y + width * 2, border = "gray50")
#	rect(xleft = 0, ybottom = y_base - y , xright = x, ytop = y_base - y + width * 3, border = "black")
}

for (i in 1:5) {
	  text(x = 0, y = y_base - i + width / 2, labels = i, pos = 2, col = "black")
}
dev.off()



pdf("TE_genome_distribution2.pdf", width = 6, height=4)
width <- 0.2
scale <- 1000000
plot(NA, type = "n", xlab = "Genome position (Mb)", yaxt = "n", bty = "n", ylab = "Chromosome", xlim = c(0, 35), ylim = c(0,6))
y_base = 5

for (i in 1:nrow(peri)) {
	x1 <- peri[i, "V2"] / scale
	x2 <- peri[i, "V3"] / scale
    y <- peri[i, "V1"]
	
	rect(xleft = x1, ybottom = y_base - y , xright = x2, ytop = y_base - y + width * 3, col = "gray80", border = NA)

}


for (i in 1:nrow(TE)) {
	x <- TE[i, "V2"] / scale
    y <- TE[i, "V1"]
	col <- TE[i, "col"]
	y_add <- ifelse(TE[i, "type"]=='Anta', width, ifelse(TE[i, "type"]=='101', width * 2, width * 3 ))

	#lines(x = c(x,x), y = c(y_base + width - y, y_base - y), col = col, lwd = 0.2)
	lines(x = c(x,x), y = c(y_base + y_add - width - y, y_base + y_add - y), col = col, lwd = 0.5)
}

for (y in seq(1,5)) {
	x <- genome_size[y] / scale
#   rect(xleft = 0, ybottom = y_base + width - y , xright = x, ytop = y_base - y, border = "black")
    rect(xleft = 0, ybottom = y_base + width - y , xright = x, ytop = y_base - y + width * 2, border = "gray50")
	rect(xleft = 0, ybottom = y_base - y , xright = x, ytop = y_base - y + width * 3, border = "black")
}

for (i in 1:5) {
	text(x = 0, y = y_base - i + width / 2, labels = i, pos = 2, col = "black")
}
dev.off()
