library(ggplot2)
library(grid)
library(gridExtra)
library(zoo)

data = read.table("upTE_union_value_binding_H3Ac_H3K4me2_WGBS.txt",sep="\t",header=T)

data$effect = (data$hda6l12_logFC)-(data$hda6_logFC+data$l12_logFC)

data$logFC_H3Ac_effect = data$H3Ac_fold_hda6ldl12-data$H3Ac_fold_hda6-data$H3Ac_fold_ldl12
data$logFC_H3K4me2_effect = data$H3K4me2_fold_hda6ldl12-data$H3K4me2_fold_hda6-data$H3K4me2_fold_ldl12
data$CG_effect = data$delta_CG_hda6ldl12_col-data$delta_CG_hda6_col-data$delta_CG_ldl12_col
data$CHG_effect = data$delta_CHG_hda6ldl12_col-data$delta_CHG_hda6_col-data$delta_CHG_ldl12_col
data$CHH_effect = data$delta_CHH_hda6ldl12_col-data$delta_CHH_hda6_col-data$delta_CHH_ldl12_col


rank_data = data[order(data$effect),]

anta = length(rank_data[rank_data$effect<0,1])

#expression
#-9~11
Zoo_rank_data0 = zoo(rank_data$effect)
Zoo_rank_data0_1 = rollapply(Zoo_rank_data0,50, mean)

#H3K4me2
#-2~3
Zoo_rank_data = zoo(rank_data$logFC_H3K4me2_effect)
Zoo_rank_data_1 = rollapply(Zoo_rank_data, 50, mean)

#H3Ac
#-3~
Zoo_rank_data2 = zoo(rank_data$logFC_H3Ac_effect)
Zoo_rank_data2_1 = rollapply(Zoo_rank_data2, 50, mean)


#CG
#-43~78
Zoo_rank_data3 = zoo(rank_data$CG_effect)
Zoo_rank_data3_1 = rollapply(Zoo_rank_data3, 50, mean)



#CHG
Zoo_rank_data4 = zoo(rank_data$CHG_effect)
Zoo_rank_data4_1 = rollapply(Zoo_rank_data4, 50, mean)


#CHH
Zoo_rank_data5 = zoo(rank_data$CHH_effect)
Zoo_rank_data5_1 = rollapply(Zoo_rank_data5, 50, mean)



row = seq(1:length(Zoo_rank_data5_1))
data2 = data.frame(row,Zoo_rank_data0_1, Zoo_rank_data_1,Zoo_rank_data2_1,Zoo_rank_data3_1,Zoo_rank_data4_1,Zoo_rank_data5_1)
colnames(data2) = c("row","expression","H3K4me2","H3Ac","CG","CHG","CHH")


write.table(data2,"ranking_expression_change_plot_mC_H3K4me2_H3Ac_change_smooth_value_net_new.txt",sep="\t",quote=F,row.names=F,col.names=T)

###### expression

png("ranking_expression_change_smooth_net.png", width=1200, height=800, res=300) # start export
ggplot(data2,aes(x=row))+
        geom_hline(yintercept = 0,color = "black", linetype="dashed",size=0.5)+
        geom_vline(xintercept = anta,color = "grey", size=0.5)+
        #geom_line(aes(y=H3K4me2),color="#F17E08")+
        #geom_line(aes(y=H3Ac),color="#7E0002")+
        geom_line(aes(y=expression*1/3),color="black")+
        xlab("Antagonism  -> Synergism")+
        #ylab("histone modification changes")+
        #scale_y_continuous(name = paste("Histone modifications FC"),limits=c(-2.5,0.5),sec.axis = sec_axis(~.*9/2,name = "Expression FC"))+
        scale_y_continuous(name = paste("Histone modifications logFC"),sec.axis = sec_axis(~.*3/1,name = "expression logFC"))+  
        theme_classic()
dev.off()



###### expression and histone changes
png("ranking_expression_change_plot_H3K4me2_H3Ac_change_smooth_net.png", width=1200, height=800, res=300) # start export
ggplot(data2,aes(x=row))+
        geom_hline(yintercept = 0,color = "black", linetype="dashed",size=0.5)+
        geom_vline(xintercept = anta,color = "grey", size=0.5)+
        geom_line(aes(y=H3K4me2),color="#F17E08")+
        geom_line(aes(y=H3Ac),color="#7E0002")+
        geom_line(aes(y=expression*1/3),color="black")+
        xlab("Antagonism  -> Synergism")+
        ylab("histone modification changes")+
        #scale_y_continuous(name = paste("Histone modifications FC"),limits=c(-2.5,0.5),sec.axis = sec_axis(~.*9/2,name = "Expression FC"))+
        scale_y_continuous(name = paste("Histone modifications logFC"),sec.axis = sec_axis(~.*3/1,name = "expression logFC"))+  
        theme_classic()
dev.off()


###### expression and methylation changes
png("ranking_expression_change_plot_mC_change_smooth_net.png", width=1200, height=800, res=300) # start export
ggplot(data2,aes(x=row))+
        geom_hline(yintercept = 0,color = "black", linetype="dashed",size=0.5)+
        geom_vline(xintercept = anta,color = "grey", size=0.5)+
        geom_line(aes(y=CG),color="#9ACC98")+
        geom_line(aes(y=CHG),color="#87CEFA")+
        geom_line(aes(y=CHH),color="#FC6665")+
        geom_line(aes(y=expression*3),color="black")+
        xlab("Antagonism  -> Synergism")+
        ylab("histone modification changes")+
        #scale_y_continuous(name = paste("Histone modifications FC"),limits=c(-2.5,0.5),sec.axis = sec_axis(~.*9/2,name = "Expression FC"))+
        scale_y_continuous(name = paste("delta methylation levels"),sec.axis = sec_axis(~.*1/3,name = "expression logFC"))+
        theme_classic()
dev.off()

