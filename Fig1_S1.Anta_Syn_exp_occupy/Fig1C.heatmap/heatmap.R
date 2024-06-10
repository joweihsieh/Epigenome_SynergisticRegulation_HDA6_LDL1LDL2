library(ggplot2)
library(grid)
library(gridExtra)
library(zoo)
library(ComplexHeatmap)
require(circlize)
library(gplots)
library(RColorBrewer)


DETE = read.table("upTE_union_value_binding_H3Ac_H3K4me2_WGBS_type_original.txt",sep="\t",header=T)
DETE_exp = DETE[,c("C_mean","L12_mean", "a_mean", "L12a_mean")]
DETE_exp_2 = DETE_exp[rev(order(DETE_exp[,"a_mean"])),]


Zoo_rank_data0=zoo(DETE_exp_2$C_mean)
Zoo_rank_data0_1=rollapply(Zoo_rank_data0, 10, mean)

Zoo_rank_data1=zoo(DETE_exp_2$L12_mean)
Zoo_rank_data1_1=rollapply(Zoo_rank_data1, 10, mean)


Zoo_rank_data2=zoo(DETE_exp_2$a_mean)
Zoo_rank_data2_1=rollapply(Zoo_rank_data2, 10, mean)

Zoo_rank_data3=zoo(DETE_exp_2$L12a_mean)
Zoo_rank_data3_1=rollapply(Zoo_rank_data3, 10, mean)


row = seq(1:length(Zoo_rank_data0_1))
data2 = data.frame(row, Zoo_rank_data0_1, Zoo_rank_data1_1,Zoo_rank_data2_1,Zoo_rank_data3_1)
colnames(data2)=c("row","WT","LDL12","HDA6","HDA6LDL12")

data2_filtered = data2
data2_filtered[data2_filtered < 0] = NA


data2_log = log(data2_filtered[,c(2:5)] + 1, 2)
effect2 = (data2_log[,"HDA6LDL12"]-data2_log[,"WT"])-((data2_log[,"HDA6"]-data2_log[,"WT"])+(data2_log[,"LDL12"]-data2_log[,"WT"]))


Mydata3_1 = data.matrix(data2_log)
rownames(Mydata3_1) = rep('', nrow(data2_log))

#save rstudio pdf 10*15
mycol = colorpanel(n=40,low="white",mid="#F6CF7F",high="darkred")
png("heatmap_upDETE_col_mutants_barchart_net_effect_sortby_hda6_with_line_noeffect.png", width=1200, height=800, res=300) 
heatmap.2(as.matrix(Mydata3_1[,c(4,3,2,1)]),key.title=NA,labRow=NA,col=mycol,trace="column",tracecol = "#797979",sepwidth=c(0.02,0.02,0.02),cex.lab=.4,mar=c(5,2),dendrogram='none', Rowv=FALSE, Colv=FALSE,colsep=c(1:8))
dev.off()



data2_log$effect = effect2
mycol2 = colorpanel(n=40,low="white",mid="white",high="white")
png("heatmap_upDETE_col_mutants_barchart_net_effect_sortby_hda6_with_effect.png", width=1200, height=800, res=300) 
heatmap.2(as.matrix(data2_log[,c(5,5)]),key.title=NA,labRow=NA,col=mycol2,trace="column",tracecol = "#797979",sepwidth=c(0.02,0.02,0.02),cex.lab=.4,mar=c(5,2),dendrogram='none', Rowv=FALSE, Colv=FALSE,colsep=c(1:8))
dev.off()


### manually combine these two plots




