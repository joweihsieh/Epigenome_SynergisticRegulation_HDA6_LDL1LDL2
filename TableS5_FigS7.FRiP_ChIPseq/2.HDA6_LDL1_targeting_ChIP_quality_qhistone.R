#cp /work1/home/hsieh/9_ProteinChIPs/PeakcallFRip.log /work1/home/joweihsieh/20190731_HDA6/ChIP-seq/20231207_data_quality/qhistone
#cp /work1/home/hsieh/9_ProteinChIPs/PeakcallFRip-hal.log /work1/home/joweihsieh/20190731_HDA6/ChIP-seq/20231207_data_quality/qhistone
#setwd("/Users/joweihsieh/Dropbox/Prof.Chen/project/HDA6/20220516_all_data_ready/results/quality/qhistone")

Data1 <- read.table("PeakcallFRip.log")
Data2 <- read.table("PeakcallFRip-hal.log")


Data1_2 <- rbind(Data1, Data2)
Data1_2 <- Data1_2[!duplicated(Data1_2$V1), ]
#1270
Data1_2_order <- Data1_2[rev(order(Data1_2$V4)),]
write.table(Data1_2_order,"FRiP_qhistone_TFChIP.txt",sep="\t",quote=F,row.names=F,col.names=T)

FRiP_ours <- read.table("FRiP.txt")
for (i in 1:nrow(FRiP_ours)){
#for (i in 1) {
  Data1_2_order$Large_Small <- ifelse(Data1_2_order$V4 > (FRiP_ours[i,1]/100), "Large", "Small")
  FRiP_ours[i,2] <- (table(Data1_2_order$Large_Small)[1]+1)/(nrow(Data1_2_order)+1)

}
colnames(FRiP_ours) <- c("our_FRiP","top_ranking")
write.table(FRiP_ours,"FRiP_our_qhistone_TFChIP.txt",sep="\t",quote=F,row.names=F,col.names=T)



# plotting
ggplot(Data1_2_order, aes(x = V4*100)) +
    geom_density(fill = "grey", color = "black") +
    geom_vline(xintercept = 0.145*100, color = "darkblue", linetype = "dashed", size = 0.5) +
    geom_vline(xintercept = 0.265*100, color = "darkblue", linetype = "dashed", size = 0.5) +
    geom_vline(xintercept = 0.281*100, color = "darkred", linetype = "dashed", size = 0.5) +
    geom_vline(xintercept = 0.1*100, color = "darkred", linetype = "dashed", size = 0.5) +
    labs(title = "Density Plot of FRiP in available ChIP-seq", x = "FRiP score (%)", y = "Density") +
    theme_minimal()
ggsave("DensityPlot_online_protein_ChIP_FRiP_labels_ours.png", width = 6, height = 4, dpi = 300)


ggplot(Data1_2_order, aes(x = V4 * 100)) +
    geom_histogram(binwidth = 1, fill = "white", color = "black") +
    #geom_vline(xintercept = c(0.145 * 100, 0.265 * 100, 0.281 * 100, 0.10 * 100),
    #           color = c("black", "black", "black", "black"),
    #           linetype = "dashed", size = 0.5) +
    labs(title = "Histogram of FRiP in available ChIP-seq", x = "FRiP score (%)", y = "Counts") +
    theme_minimal()
ggsave("Histogram_online_protein_ChIP_FRiP_labels_ours.png", width = 6, height = 4, dpi = 300)



