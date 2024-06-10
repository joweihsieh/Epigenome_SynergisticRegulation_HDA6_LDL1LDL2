library(ggplot2)
library(grid)
library(gridExtra)
library(zoo)

setwd("/Users/joweihsieh/Dropbox/Prof.Chen/project/HDA6/20220516_all_data_ready/results/ranking_exp/20230714_信賴區間")


data=read.table("upTE_union_value_binding_H3Ac_H3K4me2_WGBS.txt",sep="\t",header=T)

data$effect=(data$hda6l12_logFC)-(data$hda6_logFC+data$l12_logFC)

data$logFC_H3Ac_effect=data$H3Ac_fold_hda6ldl12-data$H3Ac_fold_hda6-data$H3Ac_fold_ldl12
data$logFC_H3K4me2_effect=data$H3K4me2_fold_hda6ldl12-data$H3K4me2_fold_hda6-data$H3K4me2_fold_ldl12
data$CG_effect=data$delta_CG_hda6ldl12_col-data$delta_CG_hda6_col-data$delta_CG_ldl12_col
data$CHG_effect=data$delta_CHG_hda6ldl12_col-data$delta_CHG_hda6_col-data$delta_CHG_ldl12_col
data$CHH_effect=data$delta_CHH_hda6ldl12_col-data$delta_CHH_hda6_col-data$delta_CHH_ldl12_col


rank_data=data[order(data$effect),]

# expression

Zoo_rank_data0=zoo(rank_data$effect)
Zoo_rank_data0_1=rollapply(Zoo_rank_data0,50, mean)


table_for_CI = matrix(0, nrow=51, ncol=length(Zoo_rank_data0_1))
table_for_CI[1,] = Zoo_rank_data0_1


for (i in 1:ncol(table_for_CI)){
	# 1:50, 2:51
	table_for_CI[2:51,i] = Zoo_rank_data0[(0+i):(49+i)]
}

table_for_CI_t = data.frame(t(table_for_CI))
colnames(table_for_CI_t)[1]=c("mean")
table_for_CI_t$row=seq(1:nrow(table_for_CI_t))

write.table(table_for_CI_t,"ranking_expression_change_smooth_value_net_new_originl.txt",sep="\t",quote=F,row.names=F,col.names=T)


#############
#H3K4me2
#-2~3
Zoo_rank_data=zoo(rank_data$logFC_H3K4me2_effect)
Zoo_rank_data_1=rollapply(Zoo_rank_data, 50, mean)

#H3Ac
#-3~
Zoo_rank_data2=zoo(rank_data$logFC_H3Ac_effect)
Zoo_rank_data2_1=rollapply(Zoo_rank_data2, 50, mean)


#CG
#-43~78
Zoo_rank_data3=zoo(rank_data$CG_effect)
Zoo_rank_data3_1=rollapply(Zoo_rank_data3, 50, mean)



#CHG
Zoo_rank_data4=zoo(rank_data$CHG_effect)
Zoo_rank_data4_1=rollapply(Zoo_rank_data4, 50, mean)


#CHH
Zoo_rank_data5=zoo(rank_data$CHH_effect)
Zoo_rank_data5_1=rollapply(Zoo_rank_data5, 50, mean)



row=seq(1:length(Zoo_rank_data5_1))
data2 <- data.frame(row,Zoo_rank_data0_1, Zoo_rank_data_1,Zoo_rank_data2_1,Zoo_rank_data3_1,Zoo_rank_data4_1,Zoo_rank_data5_1)
colnames(data2)=c("row","expression","H3K4me2","H3Ac","CG","CHG","CHH")


write.table(data2,"ranking_expression_change_plot_mC_H3K4me2_H3Ac_change_smooth_value_net_new.txt",sep="\t",quote=F,row.names=F,col.names=T)

############### https://www.roelpeters.be/add-confidence-interval-line-ggplot/


library(magrittr)
library(dplyr)
library(ggplot2)


data3=read.table("ranking_expression_change_smooth_value_net_new_originl.txt",sep="\t",header=T)
data3$sd = apply(data3[2:51], 1 , function(x) sd(x))
data3$mean_high = data3$mean+data3$sd
data3$mean_low = data3$mean-data3$sd


data4=read.table("ranking_expression_change_plot_mC_H3K4me2_H3Ac_change_smooth_value_net_new.txt",sep="\t",header=T)

combined_data = cbind(data3,data4[2:ncol(data4)])
##### pearson


pearson_corr_H3K4me2 <- cor(combined_data$expression, combined_data$H3K4me2, method = "pearson")
pearson_corr_H3Ac <- cor(combined_data$expression, combined_data$H3Ac, method = "pearson")
pearson_corr_CG <- cor(combined_data_anta$expression, combined_data_anta$CG, method = "pearson")
pearson_corr_CHG <- cor(combined_data_syn$expression, combined_data_syn$CHG, method = "pearson")
pearson_corr_CHH <- cor(combined_data_syn$expression, combined_data_syn$CHH, method = "pearson")


##### spearman
combined_data$rank_Exp <- c(1:dim(combined_data)[1])
combined_data[order(combined_data$H3K4me2),"rank_H3K4me2"] <- c(1:dim(combined_data)[1])
combined_data[order(combined_data$H3Ac),"rank_H3Ac"] <- c(1:dim(combined_data)[1])
combined_data[order(combined_data$CG),"rank_CG"] <- c(1:dim(combined_data)[1])
combined_data[order(combined_data$CHG),"rank_CHG"] <- c(1:dim(combined_data)[1])
combined_data[order(combined_data$CHH),"rank_CHH"] <- c(1:dim(combined_data)[1])

combined_data_anta <- combined_data[combined_data$expression < 0,] 
combined_data_syn <- combined_data[combined_data$expression > 0,] 

spearman_corr_H3K4me2 <- cor(combined_data$rank_Exp, combined_data$rank_H3K4me2, method = "spearman")
spearman_corr_H3Ac <- cor(combined_data$rank_Exp, combined_data$rank_H3Ac, method = "spearman")

spearman_corr_CG <- cor(combined_data_anta$rank_Exp, combined_data_anta$rank_CG, method = "spearman")
spearman_corr_CHG <- cor(combined_data_syn$rank_Exp, combined_data_syn$rank_CHG, method = "spearman")
spearman_corr_CHH <- cor(combined_data_syn$rank_Exp, combined_data_syn$rank_CHH, method = "spearman")


#####

combined_data$TRUE_NOISE = ifelse(combined_data$mean_high*combined_data$mean_low < 0,"Noise","TRUE")

Noise_length = length(combined_data[combined_data$TRUE_NOISE =="Noise",1])
# 29

Noise_number = which(combined_data$TRUE_NOISE =="Noise") 
Noise_start = Noise_number[1]
Noise_end = Noise_number[length(Noise_number)]



png("ranking_expression_change_smooth_net_CI.png", width=1200, height=800, res=300) # start export
ggplot(data3,aes(x=row))+
        geom_hline(yintercept = 0,color = "black", linetype="dashed",size=0.2)+

        geom_vline(xintercept = Noise_start,color = "blue", linetype="dashed", size=0.2)+
        geom_vline(xintercept = Noise_end,color = "blue", linetype="dashed", size=0.2)+

        geom_line(aes(y=mean*1/3),color="black",size=0.5)+
  		geom_ribbon(aes(ymin = mean_low*1/3, ymax = mean_high*1/3), alpha = 0.2)+

        xlab("Antagonism  -> Synergism")+
        scale_y_continuous(name = paste("Histone modifications logFC"),sec.axis = sec_axis(~.*3/1,name = "expression logFC"))+  
        theme_classic()
dev.off()


png("ranking_expression_change_H3Ac_H3k4me2_smooth_net_CI.png", width=1200, height=800, res=300) # start export
ggplot(combined_data,aes(x=row))+
        geom_hline(yintercept = 0,color = "black", linetype="dashed",size=0.2)+

        geom_line(aes(y=H3K4me2),color="#F17E08",size=0.5)+
        geom_line(aes(y=H3Ac),color="#7E0002",size=0.5)+

        geom_vline(xintercept = Noise_start,color = "blue", linetype="dashed", size=0.2)+
        geom_vline(xintercept = Noise_end,color = "blue", linetype="dashed", size=0.2)+

        geom_line(aes(y=mean*1/3),color="black",size=0.5)+
  		geom_ribbon(aes(ymin = mean_low*1/3, ymax = mean_high*1/3), alpha = 0.2)+

        xlab("Antagonism  -> Synergism")+
        scale_y_continuous(name = paste("Histone modifications logFC"),sec.axis = sec_axis(~.*3/1,name = "expression logFC"))+  
        theme_classic()
dev.off()

############

Noise_start_index_exp = combined_data[c(Noise_start:Noise_end),"expression"][1]
#-0.155477287

Noise_end_index_exp = combined_data[c(Noise_start:Noise_end),"expression"][Noise_length]
#0.192212498

data=read.table("upTE_union_value_binding_H3Ac_H3K4me2_WGBS.txt",sep="\t",header=T)

data$effect=(data$hda6l12_logFC)-(data$hda6_logFC+data$l12_logFC)

data$logFC_H3Ac_effect=data$H3Ac_fold_hda6ldl12-data$H3Ac_fold_hda6-data$H3Ac_fold_ldl12
data$logFC_H3K4me2_effect=data$H3K4me2_fold_hda6ldl12-data$H3K4me2_fold_hda6-data$H3K4me2_fold_ldl12
data$CG_effect=data$delta_CG_hda6ldl12_col-data$delta_CG_hda6_col-data$delta_CG_ldl12_col
data$CHG_effect=data$delta_CHG_hda6ldl12_col-data$delta_CHG_hda6_col-data$delta_CHG_ldl12_col
data$CHH_effect=data$delta_CHH_hda6ldl12_col-data$delta_CHH_hda6_col-data$delta_CHH_ldl12_col


rank_data=data[order(data$effect),]


rank_data$Cut_Keep = ifelse((rank_data$effect < Noise_end_index_exp & rank_data$effect > Noise_start_index_exp ),"Cutout","Keep")


Cutout_number = which(rank_data$Cut_Keep =="Cutout") 


################################# 2SD

library(magrittr)
library(dplyr)
library(ggplot2)


data3=read.table("ranking_expression_change_smooth_value_net_new_originl.txt",sep="\t",header=T)
data3$sd = apply(data3[2:51], 1 , function(x) sd(x))
data3$mean_high = data3$mean+data3$sd*2
data3$mean_low = data3$mean-data3$sd*2


data4=read.table("ranking_expression_change_plot_mC_H3K4me2_H3Ac_change_smooth_value_net_new.txt",sep="\t",header=T)

combined_data = cbind(data3,data4[2:ncol(data4)])


combined_data$TRUE_NOISE = ifelse(combined_data$mean_high*combined_data$mean_low < 0,"Noise","TRUE")

Noise_length = length(combined_data[combined_data$TRUE_NOISE =="Noise",1])
# 57

Noise_number = which(combined_data$TRUE_NOISE =="Noise") 
Noise_start = Noise_number[1]
Noise_end = Noise_number[length(Noise_number)]



png("ranking_expression_change_smooth_net_CI_2SD.png", width=1200, height=800, res=300) # start export
ggplot(data3,aes(x=row))+
        geom_hline(yintercept = 0,color = "black", linetype="dashed",size=0.2)+

        geom_vline(xintercept = Noise_start,color = "blue", linetype="dashed", size=0.2)+
        geom_vline(xintercept = Noise_end,color = "blue", linetype="dashed", size=0.2)+

        geom_line(aes(y=mean*1/3),color="black",size=0.5)+
  		geom_ribbon(aes(ymin = mean_low*1/3, ymax = mean_high*1/3), alpha = 0.2)+

        xlab("Antagonism  -> Synergism")+
        scale_y_continuous(name = paste("Histone modifications logFC"),sec.axis = sec_axis(~.*3/1,name = "expression logFC"))+  
        theme_classic()
dev.off()


############

Noise_start_index_exp = combined_data[c(Noise_start:Noise_end),"expression"][1]
#-0.2733835

Noise_end_index_exp = combined_data[c(Noise_start:Noise_end),"expression"][Noise_length]
#0.4001392

data=read.table("upTE_union_value_binding_H3Ac_H3K4me2_WGBS.txt",sep="\t",header=T)

data$effect=(data$hda6l12_logFC)-(data$hda6_logFC+data$l12_logFC)

data$logFC_H3Ac_effect=data$H3Ac_fold_hda6ldl12-data$H3Ac_fold_hda6-data$H3Ac_fold_ldl12
data$logFC_H3K4me2_effect=data$H3K4me2_fold_hda6ldl12-data$H3K4me2_fold_hda6-data$H3K4me2_fold_ldl12
data$CG_effect=data$delta_CG_hda6ldl12_col-data$delta_CG_hda6_col-data$delta_CG_ldl12_col
data$CHG_effect=data$delta_CHG_hda6ldl12_col-data$delta_CHG_hda6_col-data$delta_CHG_ldl12_col
data$CHH_effect=data$delta_CHH_hda6ldl12_col-data$delta_CHH_hda6_col-data$delta_CHH_ldl12_col


rank_data=data[order(data$effect),]


rank_data$Cut_Keep = ifelse((rank_data$effect < Noise_end_index_exp & rank_data$effect > Noise_start_index_exp ),"Cutout","Keep")


Cutout_number = which(rank_data$Cut_Keep =="Cutout") 

# 220~271



################################# 3SD

library(magrittr)
library(dplyr)
library(ggplot2)


data3=read.table("ranking_expression_change_smooth_value_net_new_originl.txt",sep="\t",header=T)
data3$sd = apply(data3[2:51], 1 , function(x) sd(x))
data3$mean_high = data3$mean+data3$sd*3
data3$mean_low = data3$mean-data3$sd*3


data4=read.table("ranking_expression_change_plot_mC_H3K4me2_H3Ac_change_smooth_value_net_new.txt",sep="\t",header=T)

combined_data = cbind(data3,data4[2:ncol(data4)])


combined_data$TRUE_NOISE = ifelse(combined_data$mean_high*combined_data$mean_low < 0,"Noise","TRUE")

Noise_length = length(combined_data[combined_data$TRUE_NOISE =="Noise",1])
# 82

Noise_number = which(combined_data$TRUE_NOISE =="Noise") 
Noise_start = Noise_number[1]
Noise_end = Noise_number[length(Noise_number)]



png("ranking_expression_change_smooth_net_CI_3SD.png", width=1200, height=800, res=300) # start export
ggplot(data3,aes(x=row))+
        geom_hline(yintercept = 0,color = "black", linetype="dashed",size=0.2)+

        geom_vline(xintercept = 180,color = "blue", linetype="dashed", size=0.2)+
        geom_vline(xintercept = Noise_end,color = "blue", linetype="dashed", size=0.2)+

        geom_line(aes(y=mean*1/3),color="black",size=0.5)+
  		geom_ribbon(aes(ymin = mean_low*1/3, ymax = mean_high*1/3), alpha = 0.2)+

        xlab("Antagonism  -> Synergism")+
        scale_y_continuous(name = paste("Histone modifications logFC"),sec.axis = sec_axis(~.*3/1,name = "expression logFC"))+  
        theme_classic()
dev.off()




############

Noise_start_index_exp = combined_data[c(180:Noise_end),"expression"][1]
# -0.3853322

Noise_end_index_exp = combined_data[c(180:Noise_end),"expression"][(Noise_length-2)]
# 0.5615747

data=read.table("upTE_union_value_binding_H3Ac_H3K4me2_WGBS.txt",sep="\t",header=T)

data$effect=(data$hda6l12_logFC)-(data$hda6_logFC+data$l12_logFC)

data$logFC_H3Ac_effect=data$H3Ac_fold_hda6ldl12-data$H3Ac_fold_hda6-data$H3Ac_fold_ldl12
data$logFC_H3K4me2_effect=data$H3K4me2_fold_hda6ldl12-data$H3K4me2_fold_hda6-data$H3K4me2_fold_ldl12
data$CG_effect=data$delta_CG_hda6ldl12_col-data$delta_CG_hda6_col-data$delta_CG_ldl12_col
data$CHG_effect=data$delta_CHG_hda6ldl12_col-data$delta_CHG_hda6_col-data$delta_CHG_ldl12_col
data$CHH_effect=data$delta_CHH_hda6ldl12_col-data$delta_CHH_hda6_col-data$delta_CHH_ldl12_col


rank_data=data[order(data$effect),]


rank_data$Cut_Keep = ifelse((rank_data$effect < Noise_end_index_exp & rank_data$effect > Noise_start_index_exp ),"Cutout","Keep")


Cutout_number = which(rank_data$Cut_Keep =="Cutout") 








################################# plot together with histone marks
png("ranking_expression_change_smooth_net_CI.png", width=1200, height=800, res=300) # start export
ggplot(combined_data,aes(x=row),size=1)+
        geom_hline(yintercept = 0,color = "black", linetype="dashed",size=0.5)+
        geom_line(aes(y=H3K4me2),color="#F17E08",size=1)+
        geom_line(aes(y=H3Ac),color="#7E0002",size=1)+
        geom_line(aes(y=mean*1/4),color="black",size=1)+

  		geom_ribbon(aes(ymin = mean_low*1/4, ymax = mean_high*1/3), alpha = 0.1)+

        xlab("Antagonism  -> Synergism")+
        scale_y_continuous(name = paste("Histone modifications logFC"),sec.axis = sec_axis(~.*4/1,name = "expression logFC"))+  
        theme_classic()
dev.off()




png("ranking_expression_change_plot_H3K4me2_H3Ac_change_smooth_net.png", width=1200, height=800, res=300) # start export
ggplot(data2,aes(x=row))+
        geom_hline(yintercept = 0,color = "black", linetype="dashed",size=0.5)+
        geom_vline(xintercept = anta,color = "grey", size=0.5)+
        geom_line(aes(y=expression*1/3),color="black")+
        xlab("Antagonism  -> Synergism")+
        ylab("histone modification changes")+
        #scale_y_continuous(name = paste("Histone modifications FC"),limits=c(-2.5,0.5),sec.axis = sec_axis(~.*9/2,name = "Expression FC"))+
        scale_y_continuous(name = paste("Histone modifications logFC"),sec.axis = sec_axis(~.*3/1,name = "expression logFC"))+  
        theme_classic()
dev.off()


