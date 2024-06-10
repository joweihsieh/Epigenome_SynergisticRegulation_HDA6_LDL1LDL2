

CHH= read.table("UnionSite_Region_methylation.CHH.txt",sep="\t",header=T)

start=7
replicates=2
treatments=5
ncols=seq(start,(ncol(CHH)-1),replicates)

for (i in ncols){
	CHH[[i]]=CHH[[i]]+0.00001
}


samples_interval=seq(replicates,(replicates*(treatments-1)),2)

name=c("hda6_C","ldl12_C","hda6ldl12_C","jmj28_C")

#CHH$p_value=apply(CHH[,7:10],1,function(x){t.test(x[1:2], x[3:4], alternative="two.sided")$p.value})

for (i in (samples_interval)){
	colnames=paste0("p_value_",name[(i/2)])
	CHH[[colnames]]=apply(CHH[,c(start,start+1,start+i,start+i+1)],1,function(x){t.test(x[1:2],x[3:4],alternative="two.sided")$p.value})
}


for (i in ncols){
	CHH[[i]]=CHH[[i]]-0.00001
}


for (i in ncols){
	colnames=paste0(colnames(CHH)[i],colnames(CHH)[i+1])
	CHH[,colnames]=apply(CHH[,c(i,i+1)],1,mean)
}

#CHH=CHH[,1:25]
control_mean_nth=ncol(CHH)-treatments+1


for (i in ((control_mean_nth+1):ncol(CHH))){
	colnames=paste0(colnames(CHH)[i],"_",colnames(CHH)[control_mean_nth])
	CHH[,colnames]=CHH[,i]-CHH[,control_mean_nth]
}

write.table(CHH,"Common_regions_CHH_200bp_deltaM_pvalues.txt",sep="\t",col.names=T,row.names=F,quote=F)

