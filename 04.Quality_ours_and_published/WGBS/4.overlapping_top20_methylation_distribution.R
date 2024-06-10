

#ln -s /work4/home/joweihsieh/20190731_HDA6/WGBS/DMR/200bp/Common_regions_C*_200bp_deltaM_pvalues.txt /work4/home/joweihsieh/20240424_At_meth_ATAC_pub/meth_raw/8.QC


WGBS <- read.table("Common_regions_CG_200bp_deltaM_pvalues.txt", header=T)
WGBS_CHG <- read.table("Common_regions_CHG_200bp_deltaM_pvalues.txt", header=T)
WGBS_CHH <- read.table("Common_regions_CHH_200bp_deltaM_pvalues.txt", header=T)


########################
# our data - density plot for delta methylation between replicates
########################

generate_densityplot_and_tests <- function(variable1, variable2, code, filename) {

  ks_test <- ks.test(variable1, variable2)
  print(ks_test)

  difference <- variable1 - variable2
  
  density_plot <- ggplot() +
    geom_density(aes(x = difference), fill = "skyblue", alpha = 0.5) +
    labs(title = paste0("Density Plot of Differences between ", code, " in Rep1 and Rep2"),
         x = paste0("delta methylation %"),
         y = "Density") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15),
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      legend.text = element_text(size = 15),
      legend.title = element_text(size = 15)
    )
  
  ggsave(filename, plot = density_plot, width = 8, height = 8, dpi = 300)
}


generate_densityplot_and_tests(WGBS$C1, WGBS$C2, "WT", "dist_WT.png")
generate_densityplot_and_tests(WGBS$a1, WGBS$a2, "hda6", "dist_hda6.png")
generate_densityplot_and_tests(WGBS$X121, WGBS$X122, "ldl1/2", "dist_ldl12.png")
generate_densityplot_and_tests(WGBS$a121, WGBS$a122, "hda6/ldl1/2", "dist_hda6ldl12.png")

generate_densityplot_and_tests(WGBS_CHG$C1, WGBS_CHG$C2, "WT", "dist_WT.png")
generate_densityplot_and_tests(WGBS_CHG$a1, WGBS_CHG$a2, "hda6", "dist_hda6.png")
generate_densityplot_and_tests(WGBS_CHG$X121, WGBS_CHG$X122, "ldl1/2", "dist_ldl12.png")
generate_densityplot_and_tests(WGBS_CHG$a121, WGBS_CHG$a122, "hda6/ldl1/2", "dist_hda6ldl12.png")

generate_densityplot_and_tests(WGBS_CHH$C1, WGBS_CHH$C2, "WT", "dist_WT.png")
generate_densityplot_and_tests(WGBS_CHH$a1, WGBS_CHH$a2, "hda6", "dist_hda6.png")
generate_densityplot_and_tests(WGBS_CHH$X121, WGBS_CHH$X122, "ldl1/2", "dist_ldl12.png")
generate_densityplot_and_tests(WGBS_CHH$a121, WGBS_CHH$a122, "hda6/ldl1/2", "dist_hda6ldl12.png")

WGBS$ID = paste0(WGBS$Chr, "_", WGBS$RegionStart, "_", WGBS$RegionEnd)
WGBS_CHG$ID = paste0(WGBS_CHG$Chr, "_", WGBS_CHG$RegionStart, "_", WGBS_CHG$RegionEnd)
WGBS_CHH$ID = paste0(WGBS_CHH$Chr, "_", WGBS_CHH$RegionStart, "_", WGBS_CHH$RegionEnd)


########################
# our data - overlapping
########################

#top20_num <- 83974

calculate_overlap <- function(data, variable1, variable2) {
  data1 <- data[order(data[[variable1]]),]
  data2 <- data[order(data[[variable2]]),]
 
  threshold_variable1 <- quantile(data1[[variable1]], probs = 0.8)
  threshold_variable2 <- quantile(data2[[variable2]], probs = 0.8)
  
  ids_in_threshold_variable1 <- data1[data1[[variable1]] >= threshold_variable1[[1]], "ID"]
  ids_in_threshold_variable2 <- data2[data2[[variable2]] >= threshold_variable2[[1]], "ID"]
  top20_num <- sum(data1[[variable1]] >= threshold_variable1[[1]])

  overlap_count <- length(intersect(ids_in_threshold_variable1, ids_in_threshold_variable2))
  
  return(list(overlap_count = overlap_count, top20_num = top20_num))

}

data <- WGBS
variable1 <- "C1"
variable2 <- "C2"

overlap <- calculate_overlap(WGBS, "C1", "C2")
print(paste("Number of overlapping IDs in top 20% of C1 and C2:", round(overlap[[1]][1]/overlap[[2]][1]*100,2)))

overlap2 <- calculate_overlap(WGBS, "a1", "a2")
print(paste("Number of overlapping IDs in top 20% of a1 and a2:", round(overlap2[[1]][1]/overlap2[[2]][1]*100,2)))

overlap3 <- calculate_overlap(WGBS, "X121", "X122")
print(paste("Number of overlapping IDs in top 20% of a1 and a2:", round(overlap3[[1]][1]/overlap3[[2]][1]*100,2)))

overlap4 <- calculate_overlap(WGBS, "a121", "a122")
print(paste("Number of overlapping IDs in top 20% of a1 and a2:", round(overlap4[[1]][1]/overlap4[[2]][1]*100,2)))

#####
overlap <- calculate_overlap(WGBS_CHG, "C1", "C2")
print(paste("Number of overlapping IDs in top 20% of C1 and C2:", round(overlap[[1]][1]/overlap[[2]][1]*100,2)))

overlap2 <- calculate_overlap(WGBS_CHG, "a1", "a2")
print(paste("Number of overlapping IDs in top 20% of a1 and a2:", round(overlap2[[1]][1]/overlap2[[2]][1]*100,2)))

overlap3 <- calculate_overlap(WGBS_CHG, "X121", "X122")
print(paste("Number of overlapping IDs in top 20% of a1 and a2:", round(overlap3[[1]][1]/overlap3[[2]][1]*100,2)))

overlap4 <- calculate_overlap(WGBS_CHG, "a121", "a122")
print(paste("Number of overlapping IDs in top 20% of a1 and a2:", round(overlap4[[1]][1]/overlap4[[2]][1]*100,2)))

#####
overlap <- calculate_overlap(WGBS_CHH, "C1", "C2")
print(paste("Number of overlapping IDs in top 20% of C1 and C2:", round(overlap[[1]][1]/overlap[[2]][1]*100,2)))

overlap2 <- calculate_overlap(WGBS_CHH, "a1", "a2")
print(paste("Number of overlapping IDs in top 20% of a1 and a2:", round(overlap2[[1]][1]/overlap2[[2]][1]*100,2)))

overlap3 <- calculate_overlap(WGBS_CHH, "X121", "X122")
print(paste("Number of overlapping IDs in top 20% of a1 and a2:", round(overlap3[[1]][1]/overlap3[[2]][1]*100,2)))

overlap4 <- calculate_overlap(WGBS_CHH, "a121", "a122")
print(paste("Number of overlapping IDs in top 20% of a1 and a2:", round(overlap4[[1]][1]/overlap4[[2]][1]*100,2)))

########################
# published data - overlapping
########################
cd /work4/home/joweihsieh/20240424_At_meth_ATAC_pub/meth_raw/8.QC
ln -s /work4/home/yenmr/project/HDAC6/WGBS_public/raw/Region_g*/ ./


#Sample <- "Region_g1"
#context <- "CG"

calculate_overlap_pub <- function(Sample, context) {
  data <- read.table(paste0("./", Sample,"/region_", context, ".txt"), header = TRUE)
  #  data <- read.table("region_CG.txt", header=T)


  data$ID <- paste0(data$Chr, "_", data$Start, "_", data$End)
  
  data1 <- data[order(data[[5]]), ]
  data2 <- data[order(data[[6]]), ]
  
  threshold_variable1 <- quantile(data1[[5]], probs = 0.8)
  threshold_variable2 <- quantile(data2[[6]], probs = 0.8)
  
  ids_in_threshold_variable1 <- data1[data1[[5]] >= threshold_variable1[[1]], "ID"]
  ids_in_threshold_variable2 <- data2[data2[[6]] >= threshold_variable2[[1]], "ID"]
  
  overlap_count <- length(intersect(ids_in_threshold_variable1, ids_in_threshold_variable2))
  
  top20_num <- sum(data1[[5]] >= threshold_variable1)
  
  return(list(overlap_count = overlap_count, top20_num = top20_num))
}

overlap <- calculate_overlap_pub("Region_g1","CG")
print(paste("Number of overlapping IDs in top 20%:", round(overlap[[1]][1]/overlap[[2]][1]*100,2)))

overlap <- calculate_overlap_pub("Region_g2","CG")
print(paste("Number of overlapping IDs in top 20%:", round(overlap[[1]][1]/overlap[[2]][1]*100,2)))

overlap <- calculate_overlap_pub("Region_g3","CG")
print(paste("Number of overlapping IDs in top 20%:", round(overlap[[1]][1]/overlap[[2]][1]*100,2)))

overlap <- calculate_overlap_pub("Region_g4","CG")
print(paste("Number of overlapping IDs in top 20%:", round(overlap[[1]][1]/overlap[[2]][1]*100,2)))

overlap <- calculate_overlap_pub("Region_g5","CG")
print(paste("Number of overlapping IDs in top 20%:", round(overlap[[1]][1]/overlap[[2]][1]*100,2)))

overlap <- calculate_overlap_pub("Region_g6","CG")
print(paste("Number of overlapping IDs in top 20%:", round(overlap[[1]][1]/overlap[[2]][1]*100,2)))


overlap <- calculate_overlap_pub("Region_g7","CG")
print(paste("Number of overlapping IDs in top 20%:", round(overlap[[1]][1]/overlap[[2]][1]*100,2)))

overlap <- calculate_overlap_pub("Region_g8","CG")
print(paste("Number of overlapping IDs in top 20%:", round(overlap[[1]][1]/overlap[[2]][1]*100,2)))

overlap <- calculate_overlap_pub("Region_g9","CG")
print(paste("Number of overlapping IDs in top 20%:", round(overlap[[1]][1]/overlap[[2]][1]*100,2)))


overlap <- calculate_overlap_pub("Region_g10","CG")
print(paste("Number of overlapping IDs in top 20%:", round(overlap[[1]][1]/overlap[[2]][1]*100,2)))


###
overlap <- calculate_overlap_pub("Region_g1","CHG")
print(paste("Number of overlapping IDs in top 20%:", round(overlap[[1]][1]/overlap[[2]][1]*100,2)))

overlap <- calculate_overlap_pub("Region_g2","CHG")
print(paste("Number of overlapping IDs in top 20%:", round(overlap[[1]][1]/overlap[[2]][1]*100,2)))

overlap <- calculate_overlap_pub("Region_g3","CHG")
print(paste("Number of overlapping IDs in top 20%:", round(overlap[[1]][1]/overlap[[2]][1]*100,2)))

overlap <- calculate_overlap_pub("Region_g4","CHG")
print(paste("Number of overlapping IDs in top 20%:", round(overlap[[1]][1]/overlap[[2]][1]*100,2)))

overlap <- calculate_overlap_pub("Region_g5","CHG")
print(paste("Number of overlapping IDs in top 20%:", round(overlap[[1]][1]/overlap[[2]][1]*100,2)))

overlap <- calculate_overlap_pub("Region_g6","CHG")
print(paste("Number of overlapping IDs in top 20%:", round(overlap[[1]][1]/overlap[[2]][1]*100,2)))


overlap <- calculate_overlap_pub("Region_g7","CHG")
print(paste("Number of overlapping IDs in top 20%:", round(overlap[[1]][1]/overlap[[2]][1]*100,2)))

overlap <- calculate_overlap_pub("Region_g8","CHG")
print(paste("Number of overlapping IDs in top 20%:", round(overlap[[1]][1]/overlap[[2]][1]*100,2)))

overlap <- calculate_overlap_pub("Region_g9","CHG")
print(paste("Number of overlapping IDs in top 20%:", round(overlap[[1]][1]/overlap[[2]][1]*100,2)))


overlap <- calculate_overlap_pub("Region_g10","CHG")
print(paste("Number of overlapping IDs in top 20%:", round(overlap[[1]][1]/overlap[[2]][1]*100,2)))


###
overlap <- calculate_overlap_pub("Region_g1","CHH")
print(paste("Number of overlapping IDs in top 20%:", round(overlap[[1]][1]/overlap[[2]][1]*100,2)))

overlap <- calculate_overlap_pub("Region_g2","CHH")
print(paste("Number of overlapping IDs in top 20%:", round(overlap[[1]][1]/overlap[[2]][1]*100,2)))

overlap <- calculate_overlap_pub("Region_g3","CHH")
print(paste("Number of overlapping IDs in top 20%:", round(overlap[[1]][1]/overlap[[2]][1]*100,2)))

overlap <- calculate_overlap_pub("Region_g4","CHH")
print(paste("Number of overlapping IDs in top 20%:", round(overlap[[1]][1]/overlap[[2]][1]*100,2)))

overlap <- calculate_overlap_pub("Region_g5","CHH")
print(paste("Number of overlapping IDs in top 20%:", round(overlap[[1]][1]/overlap[[2]][1]*100,2)))

overlap <- calculate_overlap_pub("Region_g6","CHH")
print(paste("Number of overlapping IDs in top 20%:", round(overlap[[1]][1]/overlap[[2]][1]*100,2)))


overlap <- calculate_overlap_pub("Region_g7","CHH")
print(paste("Number of overlapping IDs in top 20%:", round(overlap[[1]][1]/overlap[[2]][1]*100,2)))

overlap <- calculate_overlap_pub("Region_g8","CHH")
print(paste("Number of overlapping IDs in top 20%:", round(overlap[[1]][1]/overlap[[2]][1]*100,2)))

overlap <- calculate_overlap_pub("Region_g9","CHH")
print(paste("Number of overlapping IDs in top 20%:", round(overlap[[1]][1]/overlap[[2]][1]*100,2)))


overlap <- calculate_overlap_pub("Region_g10","CHH")
print(paste("Number of overlapping IDs in top 20%:", round(overlap[[1]][1]/overlap[[2]][1]*100,2)))














