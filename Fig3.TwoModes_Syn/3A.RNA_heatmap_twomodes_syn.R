# RNA-seq heatmap for double_locked and HDA6_predomimant

library(ComplexHeatmap)
require(circlize)
library(gplots)
library(RColorBrewer)

DETE <- read.table("upTE_union_value_binding_H3Ac_H3K4me2_WGBS_type_original_101_001_cluster_20220804.txt", sep = "\t", header = T)

DETE_syn <- DETE[DETE$types != "strong_weak_anta",]
DETE_syn_exp <- DETE_syn[,c("Row.names","C_mean","L12_mean", "a_mean", "L12a_mean","DETE_exp_log.effect2","types")]

DETE_syn_exp_HDA6_predomimant <- DETE_syn_exp[DETE_syn_exp$types == "1_0_1",]
DETE_syn_exp_double_locked <- DETE_syn_exp[DETE_syn_exp$types == "0_0_1",]

col_runif <- colorRamp2(c(0,5,10), c("white","#F6CF7F","darkred"))



ha_e_1 <- Heatmap(as.matrix(DETE_syn_exp_HDA6_predomimant[,c(2:5)]), name = "HDA6_predomimant", col = col_runif, column_dend_height = unit(1, "cm"),
                 row_title="",cluster_columns = F, column_names_side = "bottom", cluster_rows = T, show_row_names = FALSE)


ha_e_2 <- Heatmap(as.matrix(DETE_syn_exp_double_locked[,c(2:5)]), name = "double_locked", col = col_runif, column_dend_height = unit(1, "cm"),
                 row_title="",cluster_columns = F, column_names_side = "bottom", cluster_rows = T, show_row_names = FALSE)



png("heatmap_upDETE_col_mutants_syn_HDA6_predomimant.png", width = 1200, height = 800, res = 300) 

ha_e_1

dev.off()



png("heatmap_upDETE_col_mutants_syn_double_locked.png", width = 1200, height = 800, res = 300) 

ha_e_2

dev.off()



