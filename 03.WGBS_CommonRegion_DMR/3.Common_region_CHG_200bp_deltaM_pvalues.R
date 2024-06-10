CHG= read.table("UnionSite_Region_methylation.CHG.txt",sep="\t",header=T)

start=7
replicates=2
treatments=5
ncols=seq(start,(ncol(CHG)-1),replicates)

for (i in ncols){
	CHG[[i]]=CHG[[i]]+0.00001
}


samples_interval=seq(replicates,(replicates*(treatments-1)),2)

name=c("hda6_C","ldl12_C","hda6ldl12_C","jmj28_C")

#CHG$p_value=apply(CHG[,7:10],1,function(x){t.test(x[1:2], x[3:4], alternative="two.sided")$p.value})

for (i in (samples_interval)){
	colnames=paste0("p_value_",name[(i/2)])
	CHG[[colnames]]=apply(CHG[,c(start,start+1,start+i,start+i+1)],1,function(x){t.test(x[1:2],x[3:4],alternative="two.sided")$p.value})
}


for (i in ncols){
	CHG[[i]]=CHG[[i]]-0.00001
}


for (i in ncols){
	colnames=paste0(colnames(CHG)[i],colnames(CHG)[i+1])
	CHG[,colnames]=apply(CHG[,c(i,i+1)],1,mean)
}

#CHG=CHG[,1:25]
control_mean_nth=ncol(CHG)-treatments+1


for (i in ((control_mean_nth+1):ncol(CHG))){
	colnames=paste0(colnames(CHG)[i],"_",colnames(CHG)[control_mean_nth])
	CHG[,colnames]=CHG[,i]-CHG[,control_mean_nth]
}

write.table(CHG,"Common_regions_CHG_200bp_deltaM_pvalues.txt",sep="\t",col.names=T,row.names=F,quote=F)

