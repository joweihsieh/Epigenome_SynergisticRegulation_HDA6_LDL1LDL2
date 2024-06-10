library(ggplot2)
library(gplots)
library(RColorBrewer)
library(pheatmap)
library(dendextend)


hda6=read.table("DETE_a.txt",header=T)
l12=read.table("DETE_L12.txt",header=T)
hda6_l12=read.table("DETE_L12a.txt",header=T)

a=data.frame(Row.names=hda6[hda6$DETE=="Up","Row.names"])
b=data.frame(Row.names=l12[l12$DETE=="Up","Row.names"])
c=data.frame(Row.names=hda6_l12[hda6_l12$DETE=="Up","Row.names"])


ID_1=merge(a,b,by="Row.names",all=T)
ID=merge(ID_1,c,by="Row.names",all=T)


for (i in 1:nrow(ID)){
	if (ID[i,"Row.names"]%in%a$Row.names){
		ID[i,"hda6"]=1
	} else {ID[i,"hda6"]=0}
}

for (i in 1:nrow(ID)){
	if (ID[i,"Row.names"]%in%b$Row.names){
		ID[i,"ldl1/2"]=1
	} else {ID[i,"ldl1/2"]=0}
}

for (i in 1:nrow(ID)){
	if (ID[i,"Row.names"]%in%c$Row.names){
		ID[i,"hda6/ldl1/2"]=1
	} else {ID[i,"hda6/ldl1/2"]=0}
}



TE=read.table("Exprssed_TE.txt",header=T)

data=merge(ID,TE[,c(1,2,3,4)],by="Row.names")

write.table(data,"upTE_union_ID.txt",col.names=T,row.names=F,sep="\t",quote=F)
write.table(data[,c(5,6,7,1)],"upTE_union_ID.bed",col.names=F,row.names=F,sep="\t",quote=F)



################ tab for union upDETE
#cannot conduct this script under sal
library(magrittr)
FileList=list.files(pattern=".bam.CG.bw")

length(FileList)
list = FileList %>% sapply(function(x)strsplit(x,'-')[[1]][2])%>% sapply(function(x)strsplit(x,'\\.')[[1]][1]) 

datalist=list()
for (ff in seq(1:length(FileList))){
        File1=FileList[ff]
        File2=paste0("WGBS_CG_",list[ff],"_upTE_union.tab")
        code=paste("~/bin/bigWigAverageOverBed",File1,"upTE_union_ID.bed", File2)
        datalist[[ff]]=code
}

dataset=do.call(rbind,datalist)

write.table(dataset,"Run_CG_tab.sh",quote=F,sep=" ",row.names=F,col.names=F)


#####
FileList2=list.files(pattern=".bam.CHG.bw")

length(FileList2)
list2 = FileList2 %>% sapply(function(x)strsplit(x,'-')[[1]][2])%>% sapply(function(x)strsplit(x,'\\.')[[1]][1]) 

datalist2=list()
for (ff in seq(1:length(FileList2))){
        File1=FileList2[ff]
        File2=paste0("WGBS_CHG_",list2[ff],"_upTE_union.tab")
        code=paste("~/bin/bigWigAverageOverBed",File1,"upTE_union_ID.bed", File2)
        datalist2[[ff]]=code
}
dataset2=do.call(rbind,datalist2)

write.table(dataset2,"Run_CHG_tab.sh",quote=F,sep=" ",row.names=F,col.names=F)

#####

FileList3=list.files(pattern=".bam.CHH.bw")

length(FileList3)
list3 = FileList3 %>% sapply(function(x)strsplit(x,'-')[[1]][2])%>% sapply(function(x)strsplit(x,'\\.')[[1]][1]) 

datalist3=list()
for (ff in seq(1:length(FileList3))){
        File1=FileList3[ff]
        File2=paste0("WGBS_CHH_",list3[ff],"_upTE_union.tab")
        code=paste("~/bin/bigWigAverageOverBed",File1,"upTE_union_ID.bed", File2)
        datalist3[[ff]]=code
}
dataset3=do.call(rbind,datalist3)


write.table(dataset3,"Run_CHH_tab.sh",quote=F,sep=" ",row.names=F,col.names=F)





