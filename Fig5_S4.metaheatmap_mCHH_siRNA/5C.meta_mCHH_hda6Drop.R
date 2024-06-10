require(zoo)
WT = read.table(gzfile("TE_CHH_WTi.matrix.gz"),skip = 1)
hda6 = read.table(gzfile("TE_CHH_hda6i.matrix.gz"),skip = 1)
#hda6ldl12 = read.table(gzfile("TE_CHH_hda6ldl12l.matrix.gz"),skip = 1)
#ldl12 = read.table(gzfile("TE_CHH_ldl12l.matrix.gz"),skip = 1)
WTa = rbind(WT[c(4,7:66)],setNames(WT[c(4,67:126)],names(WT[c(4,7:66)])))
hda6a = rbind(hda6[c(4,7:66)],setNames(hda6[c(4,67:126)],names(hda6[c(4,7:66)])))
#ldl12a = rbind(ldl12[c(4,7:66)],setNames(ldl12[c(4,67:126)],names(ldl12[c(4,7:66)])),setNames(ldl12[c(4,127:186)],names(ldl12[c(4,7:66)])))
#hda6ldl12a = rbind(hda6ldl12[c(4,7:66)],setNames(hda6ldl12[c(4,67:126)],names(hda6ldl12[c(4,7:66)])),setNames(hda6ldl12[c(4,127:186)],names(hda6ldl12[c(4,7:66)])))

#cluster_WT = read.table("TE_CHH_WTi_sort.txt")
#cluster_hda6 = read.table("TE_CHH_hda6i_sort.txt")
#cluster_hda6ldl12 = read.table("TE_CHH_hda6ldl12i_sort.txt")
#cluster_ldl12 = read.table("TE_CHH_ldl12i_sort.txt")
cluster_wt_hda6 = read.table("TE_wt_hda6.txt",header=T)
cluster_wtCH_hda6CM = cluster_wt_hda6[cluster_wt_hda6$deepTools_group_wt=='cluster_1' & cluster_wt_hda6$deepTools_group_hda6=='cluster_2',]
cluster_wtCH_hda6CH = cluster_wt_hda6[cluster_wt_hda6$deepTools_group_wt=='cluster_1' & cluster_wt_hda6$deepTools_group_hda6=='cluster_1',]

WT_cluster_wtCH_hda6CM = apply(WTa[WTa$V4 %in% cluster_wtCH_hda6CM$name,2:61],2,mean,na.rm=T)
WT_cluster_wtCH_hda6CH = apply(WTa[WTa$V4 %in% cluster_wtCH_hda6CH$name,2:61],2,mean,na.rm=T)
hda6_cluster_wtCH_hda6CM = apply(hda6a[hda6a$V4 %in% cluster_wtCH_hda6CM$name,2:61],2,mean,na.rm=T)
hda6_cluster_wtCH_hda6CH = apply(hda6a[hda6a$V4 %in% cluster_wtCH_hda6CH$name,2:61],2,mean,na.rm=T)


col_WT = "black"
#col_hda6 = "#CA232A"
col_hda6 = '#4472C4'
#col_ldl12 = "#66B1FC"
#col_hda6ldl12 =  "#6F2EFF"

#col_cluster_1 = "#555B6E"
#col_cluster_2 = "#E5989B" 
#col_cluster_3 = "#FFD6BA"
#col_cluster_CHCM = 'black'
#col_cluster_CHCH = '#4472C4'


png(filename="Meta_TE_clusterCHCM_CHH.png",width=2000,height=2000,res=400)
par(mar=c(5,5,2,2))
plot(NA,ylab= expression(paste('CHH methylation (%)')),xlab=NA,col='white',cex =.6,pch=20,type='p',cex.main = 1,frame = FALSE,ylim=c(0,15),xlim=c(0,60),xaxt="n",axes=TRUE,cex.lab=1.3,main="CH->CM")
ind <- c("-2 Kb","","","2 Kb")
xx <- c(0,20,40,60)
axis(1,at=xx,labels=ind)
abline(v=20,lty=2,col='grey70')
abline(v=40,lty=2,col='grey70')
lines(rollapply(zoo(WT_cluster_wtCH_hda6CM *100), width = 3, by = 1, FUN = mean, align = "center"), col = col_WT, lwd = 5)
lines(rollapply(zoo(hda6_cluster_wtCH_hda6CM *100), width = 3, by = 1, FUN = mean, align = "center"), col = col_hda6, lwd = 5)
#lines(rollapply(zoo(WT_cluster_3), width = 3, by = 1, FUN = mean, align = "center"), col = col_cluster_3, lwd = 3)
#legend("topright",c("","swi3b","hda6"),text.font=c(1,3,3), bty='n',lty=1,col = c(col_col,col_b,col_sil1),cex=0.8,lwd = 3)
dev.off()


png(filename="Meta_TE_clusterCHCH_CHH.png",width=2000,height=2000,res=400)
par(mar=c(5,5,2,2))
plot(NA,ylab= expression(paste('CHH methylation (%)')),xlab=NA,col='white',cex =.6,pch=20,type='p',cex.main = 1,frame = FALSE,ylim=c(0,15),xlim=c(0,60),xaxt="n",axes=TRUE,cex.lab=1.3,main="CH->CH")
ind <- c("-2 Kb","","","2 Kb")
xx <- c(0,20,40,60)
axis(1,at=xx,labels=ind)
abline(v=20,lty=2,col='grey70')
abline(v=40,lty=2,col='grey70')
lines(rollapply(zoo(WT_cluster_wtCH_hda6CH *100), width = 3, by = 1, FUN = mean, align = "center"), col = col_WT, lwd = 5)
lines(rollapply(zoo(hda6_cluster_wtCH_hda6CH *100), width = 3, by = 1, FUN = mean, align = "center"), col = col_hda6, lwd = 5)
#lines(rollapply(zoo(hda6_cluster_3), width = 3, by = 1, FUN = mean, align = "center"), col = col_cluster_3, lwd = 3)
#legend("topright",c("","swi3b","hda6"),text.font=c(1,3,3), bty='n',lty=1,col = c(col_col,col_b,col_sil1),cex=0.8,lwd = 3)
dev.off()


#png(filename="Meta_TE_ldl12_smRNA_24nt_CHH_cluster.png",width=2000,height=2000,res=400)
#par(mar=c(5,5,2,2))
#plot(NA,ylab= expression(paste('24-nt SmRNA abundance (RPKM)')),xlab=NA,col='white',cex =.6,pch=20,type='p',cex.main = 1,frame = FALSE,ylim=c(0,1.2),xlim=c(0,60),xaxt="n",axes=TRUE,cex.lab=1.3,main="ldl12")
#ind <- c("-2 Kb","","","2 Kb")
#xx <- c(0,20,40,60)
#axis(1,at=xx,labels=ind)
#abline(v=20,lty=2,col='grey70')
#abline(v=40,lty=2,col='grey70')
#lines(rollapply(zoo(ldl12_cluster_1), width = 3, by = 1, FUN = mean, align = "center"), col = col_cluster_1, lwd = 3)
#lines(rollapply(zoo(ldl12_cluster_2), width = 3, by = 1, FUN = mean, align = "center"), col = col_cluster_2, lwd = 3)
#lines(rollapply(zoo(ldl12_cluster_3), width = 3, by = 1, FUN = mean, align = "center"), col = col_cluster_3, lwd = 3)
#legend("topright",c("","swi3b","hda6"),text.font=c(1,3,3), bty='n',lty=1,col = c(col_col,col_b,col_sil1),cex=0.8,lwd = 3)
#dev.off()


#png(filename="Meta_TE_hda6ldl12_smRNA_24nt_CHH_cluster.png",width=2000,height=2000,res=400)
#par(mar=c(5,5,2,2))
#plot(NA,ylab= expression(paste('24-nt SmRNA abundance (RPKM)')),xlab=NA,col='white',cex =.6,pch=20,type='p',cex.main = 1,frame = FALSE,ylim=c(0,1.2),xlim=c(0,60),xaxt="n",axes=TRUE,cex.lab=1.3,main="hda6ldl12")
#ind <- c("-2 Kb","","","2 Kb")
#xx <- c(0,20,40,60)
#axis(1,at=xx,labels=ind)
#abline(v=20,lty=2,col='grey70')
#abline(v=40,lty=2,col='grey70')
#lines(rollapply(zoo(hda6ldl12_cluster_1), width = 3, by = 1, FUN = mean, align = "center"), col = col_cluster_1, lwd = 3)
#lines(rollapply(zoo(hda6ldl12_cluster_2), width = 3, by = 1, FUN = mean, align = "center"), col = col_cluster_2, lwd = 3)
#lines(rollapply(zoo(hda6ldl12_cluster_3), width = 3, by = 1, FUN = mean, align = "center"), col = col_cluster_3, lwd = 3)
#legend("topright",c("","swi3b","hda6"),text.font=c(1,3,3), bty='n',lty=1,col = c(col_col,col_b,col_sil1),cex=0.8,lwd = 3)
#dev.off()

