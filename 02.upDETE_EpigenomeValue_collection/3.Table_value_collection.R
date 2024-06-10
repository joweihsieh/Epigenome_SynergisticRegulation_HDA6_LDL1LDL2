library(ggplot2)
library(gplots)
library(RColorBrewer)
library(pheatmap)
library(dendextend)


#################### 
# union up-DETE in mutants
#################### 

hda6 = read.table("DETE_a.txt",header = T)
l12 = read.table("DETE_L12.txt",header = T)
hda6_l12 = read.table("DETE_L12a.txt",header = T)

a = data.frame(Row.names = hda6[hda6$DETE == "Up","Row.names"])
b = data.frame(Row.names = l12[l12$DETE == "Up","Row.names"])
c = data.frame(Row.names = hda6_l12[hda6_l12$DETE == "Up","Row.names"])


ID_1 = merge(a,b,by = "Row.names",all = T)
ID = merge(ID_1,c,by = "Row.names",all = T)


for (i in 1:nrow(ID)){
        if (ID[i,"Row.names"]%in%a$Row.names){
                ID[i,"hda6"] = 1
        } else {ID[i,"hda6"] = 0}
}

for (i in 1:nrow(ID)){
        if (ID[i,"Row.names"]%in%b$Row.names){
                ID[i,"ldl1/2"] = 1
        } else {ID[i,"ldl1/2"] = 0}
}

for (i in 1:nrow(ID)){
        if (ID[i,"Row.names"]%in%c$Row.names){
                ID[i,"hda6/ldl1/2"] = 1
        } else {ID[i,"hda6/ldl1/2"] = 0}
}



TE = read.table("Exprssed_TE.txt",header = T)
data = merge(ID,TE[,c(1,2,3,4)],by = "Row.names")


#################### 
# ADD binding (but we do not use this at the end)
#################### 

hda6_b = read.table("ChIP_hda6_peaks.narrowPeak.position",header = F)
ldl1_b = read.table("ChIP_ldl1_peaks.narrowPeak.position",header = F)



for (i in 1:nrow(ID)){
        if (ID[i,"Row.names"]%in%hda6_b$V14){
                ID[i,"hda6_b"] = "Binding"
        } else {ID[i,"hda6_b"] = "noBinding"}
}



for (i in 1:nrow(ID)){
        if (ID[i,"Row.names"]%in%ldl1_b$V14){
                ID[i,"ldl1_b"] = "Binding"
        } else {ID[i,"ldl1_b"] = "noBinding"}
}

for (i in 1:nrow(ID)){
        if (ID[i,"hda6"] == 1 & ID[i,"ldl1/2"] == 0 & ID[i,"hda6/ldl1/2"] == 0){
                ID[i,"type"] = "hda6_specific"
        } else if (ID[i,"hda6"] == 0 & ID[i,"ldl1/2"] == 1 & ID[i,"hda6/ldl1/2"] == 0){
                ID[i,"type"] = "ldl1_2_specific"
        } else if (ID[i,"hda6"] == 0 & ID[i,"ldl1/2"] == 0 & ID[i,"hda6/ldl1/2"] == 1){
                ID[i,"type"] = "hda6_ldl1_2_specific"
        } else if (ID[i,"hda6"] == 1 & ID[i,"ldl1/2"] == 0 & ID[i,"hda6/ldl1/2"] == 1){
                ID[i,"type"] = "hda6_triple_common"
        } else if (ID[i,"hda6"] == 1 & ID[i,"ldl1/2"] == 1 & ID[i,"hda6/ldl1/2"] == 1){
                ID[i,"type"] = "common"
        } else {ID[i,"type"] = "None"}
}

#################### 
# ADD histone modification
#################### 

H3Ac_hda6_1 = read.table("ChIP_H3Ac_hda6_1_upTE_union.tab",sep = "\t",col.names = c("Row.names","a","b","c","d","H3Ac_hda6_1"))
H3Ac_hda6_2 = read.table("ChIP_H3Ac_hda6_2_upTE_union.tab",sep = "\t",col.names = c("Row.names","a","b","c","d","H3Ac_hda6_2"))

H3Ac_ldl12_1 = read.table("ChIP_H3Ac_ldl12_1_upTE_union.tab",sep = "\t",col.names = c("Row.names","a","b","c","d","H3Ac_ldl12_1"))
H3Ac_ldl12_2 = read.table("ChIP_H3Ac_ldl12_2_upTE_union.tab",sep = "\t",col.names = c("Row.names","a","b","c","d","H3Ac_ldl12_2"))

H3Ac_hda6ldl12_1 = read.table("ChIP_H3Ac_hda6ldl12_1_upTE_union.tab",sep = "\t",col.names = c("Row.names","a","b","c","d","H3Ac_hda6ldl12_1"))
H3Ac_hda6ldl12_2 = read.table("ChIP_H3Ac_hda6ldl12_2_upTE_union.tab",sep = "\t",col.names = c("Row.names","a","b","c","d","H3Ac_hda6ldl12_2"))


H3Ac_WT_1 = read.table("ChIP_H3Ac_WT_1_upTE_union.tab",sep = "\t",col.names = c("Row.names","a","b","c","d","H3Ac_WT_1"))
H3Ac_WT_2 = read.table("ChIP_H3Ac_WT_2_upTE_union.tab",sep = "\t",col.names = c("Row.names","a","b","c","d","H3Ac_WT_2"))

H3K4me2_hda6_1 = read.table("ChIP_H3K4me2_hda6_1_upTE_union.tab",sep = "\t",col.names = c("Row.names","a","b","c","d","H3K4me2_hda6_1"))
H3K4me2_hda6_2 = read.table("ChIP_H3K4me2_hda6_2_upTE_union.tab",sep = "\t",col.names = c("Row.names","a","b","c","d","H3K4me2_hda6_2"))

H3K4me2_ldl12_1 = read.table("ChIP_H3K4me2_ldl12_1_upTE_union.tab",sep = "\t",col.names = c("Row.names","a","b","c","d","H3K4me2_ldl12_1"))
H3K4me2_ldl12_2 = read.table("ChIP_H3K4me2_ldl12_2_upTE_union.tab",sep = "\t",col.names = c("Row.names","a","b","c","d","H3K4me2_ldl12_2"))

H3K4me2_hda6ldl12_1 = read.table("ChIP_H3K4me2_hda6ldl12_1_upTE_union.tab",sep = "\t",col.names = c("Row.names","a","b","c","d","H3K4me2_hda6ldl12_1"))
H3K4me2_hda6ldl12_2 = read.table("ChIP_H3K4me2_hda6ldl12_2_upTE_union.tab",sep = "\t",col.names = c("Row.names","a","b","c","d","H3K4me2_hda6ldl12_2"))


H3K4me2_WT_1 = read.table("ChIP_H3K4me2_WT_1_upTE_union.tab",sep = "\t",col.names = c("Row.names","a","b","c","d","H3K4me2_WT_1"))
H3K4me2_WT_2 = read.table("ChIP_H3K4me2_WT_2_upTE_union.tab",sep = "\t",col.names = c("Row.names","a","b","c","d","H3K4me2_WT_2"))

a = merge(ID,H3Ac_hda6_1[,c(1,6)],by = "Row.names")
b = merge(a,H3Ac_hda6_2[,c(1,6)],by = "Row.names")
c = merge(b,H3Ac_ldl12_1[,c(1,6)],by = "Row.names")
d = merge(c,H3Ac_ldl12_2[,c(1,6)],by = "Row.names")
e = merge(d,H3Ac_hda6ldl12_1[,c(1,6)],by = "Row.names")
f = merge(e,H3Ac_hda6ldl12_2[,c(1,6)],by = "Row.names")
g = merge(f,H3Ac_WT_1[,c(1,6)],by = "Row.names")
h = merge(g,H3Ac_WT_2[,c(1,6)],by = "Row.names")


a = merge(h,H3K4me2_hda6_1[,c(1,6)],by = "Row.names")
b = merge(a,H3K4me2_hda6_2[,c(1,6)],by = "Row.names")
c = merge(b,H3K4me2_ldl12_1[,c(1,6)],by = "Row.names")
d = merge(c,H3K4me2_ldl12_2[,c(1,6)],by = "Row.names")
e = merge(d,H3K4me2_hda6ldl12_1[,c(1,6)],by = "Row.names")
f = merge(e,H3K4me2_hda6ldl12_2[,c(1,6)],by = "Row.names")
g = merge(f,H3K4me2_WT_1[,c(1,6)],by = "Row.names")
h = merge(g,H3K4me2_WT_2[,c(1,6)],by = "Row.names")

h$H3Ac_hda6_mean = (h$H3Ac_hda6_1+h$H3Ac_hda6_2)/2
h$H3Ac_ldl12_mean = (h$H3Ac_ldl12_1+h$H3Ac_ldl12_2)/2
h$H3Ac_hda6ldl12_mean = (h$H3Ac_hda6ldl12_1+h$H3Ac_hda6ldl12_2)/2
h$H3Ac_WT_mean = (h$H3Ac_WT_1+h$H3Ac_WT_2)/2

h$H3K4me2_hda6_mean = (h$H3K4me2_hda6_1+h$H3K4me2_hda6_2)/2
h$H3K4me2_ldl12_mean = (h$H3K4me2_ldl12_1+h$H3K4me2_ldl12_2)/2
h$H3K4me2_hda6ldl12_mean = (h$H3K4me2_hda6ldl12_1+h$H3K4me2_hda6ldl12_2)/2
h$H3K4me2_WT_mean = (h$H3K4me2_WT_1+h$H3K4me2_WT_2)/2


h$H3Ac_fold_hda6 = log((h$H3Ac_hda6_mean+0.001)/(h$H3Ac_WT_mean+0.001),2)
h$H3Ac_fold_ldl12 = log((h$H3Ac_ldl12_mean+0.001)/(h$H3Ac_WT_mean+0.001),2)
h$H3Ac_fold_hda6ldl12 = log((h$H3Ac_hda6ldl12_mean+0.001)/(h$H3Ac_WT_mean+0.001),2)

h$H3K4me2_fold_hda6 = log((h$H3K4me2_hda6_mean+0.001)/(h$H3K4me2_WT_mean+0.001),2)
h$H3K4me2_fold_ldl12 = log((h$H3K4me2_ldl12_mean+0.001)/(h$H3K4me2_WT_mean+0.001),2)
h$H3K4me2_fold_hda6ldl12 = log((h$H3K4me2_hda6ldl12_mean+0.001)/(h$H3K4me2_WT_mean+0.001),2)


#original value
data_value =  h[,c(1:7,24:31)]

#################### 
# ADD WGBS
#################### 

library(magrittr)

#list file
FileList  =  list.files(patter = "WGBS_C")

CG_list = FileList[c(5,6,9,10,1,2,7,8)]
CHG_list = FileList[c(15,16,19,20,11,12,17,18)]
CHH_list = FileList[c(25,26,29,30,21,22,27,28)]


list1  =  CG_list %>% sapply(function(x)strsplit(x,'_')[[1]][2])
list2  =  CG_list %>% sapply(function(x)strsplit(x,'_')[[1]][3])
list3 = paste0(list1,"_",list2)

list4  =  CHG_list %>% sapply(function(x)strsplit(x,'_')[[1]][2])
list5  =  CHG_list %>% sapply(function(x)strsplit(x,'_')[[1]][3])
list6 = paste0(list4,"_",list5)

list7  =  CHH_list %>% sapply(function(x)strsplit(x,'_')[[1]][2])
list8  =  CHH_list %>% sapply(function(x)strsplit(x,'_')[[1]][3])
list9 = paste0(list7,"_",list8)

list_all = c(CG_list,CHG_list,CHH_list)
colnames_all = c(list3,list6,list9)

# read and combine files
datalist = list()
for(ff in 1:length(list_all)){
  File  =  list_all[ff]
  TFdata  =  read.table(File,sep = "\t",header = F)
  datalist[[ff]] = TFdata[,c(1,6)]
}

dataset = do.call(cbind,datalist)
meth_number = seq(4,ncol(dataset),2)
meth_table = dataset[,c(1:2,meth_number)]
colnames(meth_table) = c("Row.names",colnames_all)

number = length(meth_table)
# 

for (i in seq(2,ncol(meth_table),2)){
  names = paste0(colnames(meth_table)[i],"_",colnames(meth_table)[(i+1)],"_mean")
  meth_table[,names] = ((meth_table[,i]+meth_table[,(i+1)])/2)*100
}


meth_Value_table = meth_table[,c(1,26:37)]
colnames(meth_Value_table) = c("Row.names","CG_col_mean","CG_hda6_mean","CG_ldl12_mean","CG_hda6ldl12_mean","CHG_col_mean","CHG_hda6_mean","CHG_ldl12_mean","CHG_hda6ldl12_mean","CHH_col_mean","CHH_hda6_mean","CHH_ldl12_mean","CHH_hda6ldl12_mean")

##### delta
start = number + 1
treatments = 4

for (i in seq(start,ncol(meth_table),treatments)){
  name2 = paste0("delta_",colnames(meth_table)[i+1],"_",colnames(meth_table)[(i)])
  meth_table[,name2] = (meth_table[,i+1]-meth_table[,(i)])
  name3 = paste0("delta_",colnames(meth_table)[i+2],"_",colnames(meth_table)[(i)])
  meth_table[,name3] = (meth_table[,i+2]-meth_table[,(i)])
  name4 = paste0("delta_",colnames(meth_table)[i+3],"_",colnames(meth_table)[(i)])
  meth_table[,name4] = (meth_table[,i+3]-meth_table[,(i)])
}

colnames(meth_table)[(start+(treatments*3)):ncol(meth_table)] = c("delta_CG_hda6_col","delta_CG_ldl12_col","delta_CG_hda6ldl12_col","delta_CHG_hda6_col","delta_CHG_ldl12_col","delta_CHG_hda6ldl12_col","delta_CHH_hda6_col","delta_CHH_ldl12_col","delta_CHH_hda6ldl12_col")
data2 = meth_table[,c(1,(start+(treatments*3)):ncol(meth_table))]
mydata = merge(data,data2,by="Row.names")
#write.table(mydata,"upTE_union_binding_H3Ac_H3K4me2_WGBS.txt",col.names=T,row.names=F,sep="\t",quote=F)

#################### 
# ADD RNA-seq - I
#################### 

RNA = cbind(hda6[,c(1,16:19,20:23)],l12[,c(20:23)],hda6_l12[,c(20:23)])
colnames(RNA)[6:length(RNA)] = c("hda6_logFC","hda6_PValue","hda6_FDR","hda6_DETE","l12_logFC","l12_PValue","l12_FDR","l12_DETE","hda6l12_logFC","hda6l12_PValue","hda6l12_FDR","hda6l12_DETE")
mydata2 = merge(mydata, RNA, by = "Row.names")
mydata2$hda6l12_hda6 = mydata2$hda6l12_logFC - mydata2$hda6_logFC
write.table(mydata2, "upTE_union_value_binding_H3Ac_H3K4me2_WGBS.txt", col.names = T, row.names = F, sep = "\t", quote = F)


#################### 
# ADD RNA-seq - II
#################### 
RNA = hda6[,c(1,16:19)]

mydata = merge(data_value,meth_Value_table,by = "Row.names")
mydata2 = merge(mydata,RNA,by = "Row.names")


table = read.table("upTE_union_value_binding_H3Ac_H3K4me2_WGBS_type.txt",sep = "\t",header = T)
mydata3 = merge(mydata2,table[,c("Row.names","Chr","Start","End","length","strand","Syn_anta")],by = "Row.names")

write.table(mydata3,"upTE_union_value_binding_H3Ac_H3K4me2_WGBS_type_original.txt",sep = "\t",quote = F,row.names = F,col.names = T)


#################### 
# Effect calculation
#################### 
DETE = mydata3
DETE_exp_log = log((DETE[,c(28:31)]+1),2)
DETE_exp_log$effect2 = (DETE_exp_log[,"L12a_mean"]-DETE_exp_log[,"C_mean"])-((DETE_exp_log[,"a_mean"]-DETE_exp_log[,"C_mean"])+(DETE_exp_log[,"L12_mean"]-DETE_exp_log[,"C_mean"]))
DETE_1 = cbind(DETE,DETE_exp_log$effect2)

syn = DETE_1[DETE_1[, "DETE_exp_log$effect2"]> 0,]
anta = DETE_1[DETE_1[, "DETE_exp_log$effect2"]< 0,]

strong_syn = DETE_1[DETE_1[, "DETE_exp_log$effect2"]>= 1,]
strong_anta = DETE_1[DETE_1[, "DETE_exp_log$effect2"]<= (-1),]


for (i in 1:nrow(DETE_1)){
  if (DETE_1[i,"Row.names"]%in%syn$Row.names){
    DETE_1[i,"Syn_anta"]="syn"
  } else {DETE_1[i,"Syn_anta"]="anta"}
}


for (i in 1:nrow(DETE_1)){
  if (DETE_1[i,"Row.names"]%in%strong_syn$Row.names){
    DETE_1[i,"Syn_anta"]="strong_syn"
  } else if (DETE_1[i,"Row.names"]%in%strong_anta$Row.names){
    DETE_1[i,"Syn_anta"]="strong_anta"
  }
}

write.table(DETE_1,"upTE_union_value_binding_H3Ac_H3K4me2_WGBS_type_original_20220713.txt",col.names=T,row.names=F,sep="\t",quote=F)

#################### 
# Classify TEs into antagonistic, HDA6-predominant (1_0_1), and double_locked (0_0_1)
#################### 
new = read.table("upTE_union_value_binding_H3Ac_H3K4me2_WGBS_type_original_20220713.txt",sep="\t",header=T)
others = read.table("upTE_union_value_binding_H3Ac_H3K4me2_WGBS_type_effects_H3K9me2_smRNA_ATAC.txt",sep="\t",header=T)


Data = merge(new, others[,c(1, 8:22, 27:38)], by="Row.names")
Data2 = Data

Data2[,c(8:15)] = log(Data2[,c(8:15)]+1,2)
Data2[,c(28:31)] = log(Data2[,c(28:31)]+1,2)


type_1_0_1 = Data2 [Data2[,"Syn_anta"]=="syn" | Data2[,"Syn_anta"]=="strong_syn"& Data2[,"hda6"]== "1" & Data2[,"ldl1.2"]== "0" & Data2[,"hda6.ldl1.2"]== "1" ,]
type_1_0_1$types = "1_0_1"

type_0_0_1 = Data2 [Data2[,"hda6"]== "0" & Data2[,"ldl1.2"]== "0" & Data2[,"hda6.ldl1.2"]== "1" ,]
type_0_0_1$types = "0_0_1"


type_weak_strong_anta = Data2 [Data2[,"Syn_anta"]=="anta" | Data2[,"Syn_anta"]=="strong_anta",]
type_weak_strong_anta$types = "strong_weak_anta"


Data3 = rbind(type_1_0_1,type_0_0_1,type_weak_strong_anta)


write.table(Data3,"upTE_union_value_binding_H3Ac_H3K4me2_WGBS_type_original_101_001_20220802.txt",col.names=T,row.names=F,sep="\t",quote=F)


#################### 
# ADD cluster info from Fig5
#################### 


WT = read.table("TE_CHH_WTi_sort.txt",sep="\t",header=F)
colnames(WT)[c(4,13)]=c("Row.names","cluster_WT")


hda6 = read.table("TE_CHH_hda6i_sort.txt",sep="\t",header=F)
colnames(hda6)[c(4,13)]=c("Row.names","cluster_hda6")


ldl12 = read.table("TE_CHH_ldl12i_sort.txt",sep="\t",header=F)
colnames(ldl12)[c(4,13)]=c("Row.names","cluster_ldl12")


hda6ldl12 = read.table("TE_CHH_hda6ldl12i_sort.txt",sep="\t",header=F)
colnames(hda6ldl12)[c(4,13)]=c("Row.names","cluster_hda6ldl12")



Data4 = merge(Data3, WT[,c(4,13)], by ="Row.names")
Data5 = merge(Data4, hda6[,c(4,13)], by ="Row.names")
Data6 = merge(Data5, ldl12[,c(4,13)], by ="Row.names")
Data7 = merge(Data6, hda6ldl12[,c(4,13)], by ="Row.names")

table(Data7[,c("types","cluster_WT")])
table(Data7[,c("cluster_WT","cluster_hda6")])
table(Data7[,c("cluster_WT","cluster_ldl12")])
table(Data7[,c("cluster_WT","cluster_hda6ldl12")])


Data7_001= Data7[Data7$types=="0_0_1",]
Data7_101= Data7[Data7$types=="1_0_1",]
Data7_anta= Data7[Data7$types=="strong_weak_anta",]


table(Data7_001[,c("cluster_WT","cluster_hda6")])
table(Data7_001[,c("cluster_WT","cluster_ldl12")])
table(Data7_001[,c("cluster_WT","cluster_hda6ldl12")])


table(Data7_101[,c("cluster_WT","cluster_hda6")])
table(Data7_101[,c("cluster_WT","cluster_ldl12")])
table(Data7_101[,c("cluster_WT","cluster_hda6ldl12")])


table(Data7_anta[,c("cluster_WT","cluster_hda6")])
table(Data7_anta[,c("cluster_WT","cluster_ldl12")])
table(Data7_anta[,c("cluster_WT","cluster_hda6ldl12")])

write.table(Data7,"upTE_union_value_binding_H3Ac_H3K4me2_WGBS_type_original_101_001_cluster_20220804.txt",col.names=T,row.names=F,sep="\t",quote=F)



