CG= read.table("UnionSite_Region_methylation.CG.txt",sep="\t",header=T)

start=7
replicates=2
treatments=5
ncols=seq(start,(ncol(CG)-1),replicates)

for (i in ncols){
	CG[[i]]=CG[[i]]+0.00001
}


samples_interval=seq(replicates,(replicates*(treatments-1)),2)

name=c("hda6_C","ldl12_C","hda6ldl12_C","jmj28_C")

#CG$p_value=apply(CG[,7:10],1,function(x){t.test(x[1:2], x[3:4], alternative="two.sided")$p.value})

for (i in (samples_interval)){
	colnames=paste0("p_value_",name[(i/2)])
	CG[[colnames]]=apply(CG[,c(start,start+1,start+i,start+i+1)],1,function(x){t.test(x[1:2],x[3:4],alternative="two.sided")$p.value})
}


for (i in ncols){
	CG[[i]]=CG[[i]]-0.00001
}


for (i in ncols){
	colnames=paste0(colnames(CG)[i],colnames(CG)[i+1])
	CG[,colnames]=apply(CG[,c(i,i+1)],1,mean)
}

#CG=CG[,1:25]
control_mean_nth=ncol(CG)-treatments+1


for (i in ((control_mean_nth+1):ncol(CG))){
	colnames=paste0(colnames(CG)[i],"_",colnames(CG)[control_mean_nth])
	CG[,colnames]=CG[,i]-CG[,control_mean_nth]
}

write.table(CG,"Common_regions_200bp_deltaM_pvalues.txt",sep="\t",col.names=T,row.names=F,quote=F)

