######################
# get RPKM of published data
######################
#/work4/home/joweihsieh/20190731_HDA6/RNA-seq/published_data/fastq/diff_gene_indi/genes.fpkm_tracking
#nohup /work1/home/yenmr/tools/cufflinks-2.2.1.Linux_x86_64/cuffdiff -p 20 -o diff_gene -L WT,cmt23,cmt2,cmt3,drm12,drm12cmt23,drm12cmt2,drm12cmt3,suvh456 /work1/home/yenmr/genome/Arabidopsis_thaliana/Ensembl/TAIR10/Annotation/Genes/genes.gtf WT_1.bam,WT_2.bam cmt23_1.bam,cmt23_2.bam cmt2_1.bam,cmt2_2.bam cmt3_1.bam,cmt3_2.bam drm12_1.bam,drm12_2.bam drm12cmt23_1.bam,drm12cmt23_2.bam drm12cmt2_1.bam,drm12cmt2_2.bam drm12cmt3_1.bam,drm12cmt3_2.bam suvh456_1.bam,suvh456_2.bam &
nohup /work1/home/yenmr/tools/cufflinks-2.2.1.Linux_x86_64/cuffdiff -p 20 -o diff_gene_indi -L WT_1,WT_2,cmt23_1,cmt23_2,cmt2_1,cmt2_2,cmt3_1,cmt3_2,drm12_1,drm12_2,drm12cmt23_1,drm12cmt23_2,drm12cmt2_1,drm12cmt2_2,drm12cmt3_1,drm12cmt3_2,suvh456_1,suvh456_2 /work1/home/yenmr/genome/Arabidopsis_thaliana/Ensembl/TAIR10/Annotation/Genes/genes.gtf WT_1.bam WT_2.bam cmt23_1.bam cmt23_2.bam cmt2_1.bam cmt2_2.bam cmt3_1.bam cmt3_2.bam drm12_1.bam drm12_2.bam drm12cmt23_1.bam drm12cmt23_2.bam drm12cmt2_1.bam drm12cmt2_2.bam drm12cmt3_1.bam drm12cmt3_2.bam suvh456_1.bam suvh456_2.bam &


######################
# get RPKM of our data
######################
/work1/home/yenmr/tools/cufflinks-2.2.1.Linux_x86_64/cuffdiff -p 20 -o diff_gene -L C,L12,L12a,a /work1/home/yenmr/genome/Arabidopsis_thaliana/Ensembl/TAIR10/Annotation/Genes/genes.gtf C_6.bam,C_7.bam,C_7_2.bam L12_6.bam,L12_7.bam L12a_6.bam,L12a_7.bam a_6.bam,a_7.bam 
cp /work1/home/yenmr/project/HDAC6/RNA/raw/diff_gene/rpkm_hda6.txt /work4/home/joweihsieh/20240424_At_meth_ATAC_pub/RNA_raw


######################
# plot and testing
######################

library(ggplot2)
library(stats)   

# published data
Table <- read.table("genes.fpkm_tracking",sep="\t", header=T)
selected <- colnames(Table)[grepl("FPKM",colnames(Table))]
Table2 <- Table[,c("gene_id",selected)]

# our data
RPKM <- read.table("rpkm_hda6.txt", header=T)


generate_scatterplot <- function(variable1, variable2, code, filename) {
  spearman_correlation <- cor(variable1, variable2, method = "spearman")
  print(paste("Spearman correlation coefficient:", spearman_correlation))
  
  ks_test <- ks.test(variable1, variable2)
  print(ks_test)
  
  #scatterplot <- ggplot(mapping = aes(x = sort(variable1), y = sort(variable2))) + 
  scatterplot <- ggplot(mapping = aes(x = variable1, y = variable2)) + 
    geom_point() +
    geom_abline(aes(slope = 1, intercept = 0), linetype = 2) +
    labs(title = paste0("Spearman:", round(spearman_correlation, digits = 4), "; Kolmogorov-Smirnov test - p=", round(ks_test$p.value, digits = 3)),
         x = paste0("RPKM in Rep1 of ",code),
         y = paste0("RPKM in Rep2 of ",code)) +
    theme_minimal() +
   	theme_bw() +
    theme(
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15),
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      legend.text = element_text(size = 15),
      legend.title = element_text(size = 15)
    ) 	
  
  ggsave(filename, plot = scatterplot, width = 8, height = 8, dpi = 300)
}

# ours
generate_scatterplot(RPKM$C_0, RPKM$C_2, "WT", "scatterplot_WT.png")
generate_scatterplot(RPKM$C_1, RPKM$C_2, "WT", "scatterplot_WT2.png")
generate_scatterplot(RPKM$a_0, RPKM$a_1, "hda6", "scatterplot_hda6.png")
generate_scatterplot(RPKM$L12_0, RPKM$L12_1, "ldl1/2", "scatterplot_ldl12.png")
generate_scatterplot(RPKM$L12a_0, RPKM$L12a_1, "hda6/ldl1/2", "scatterplot_hda6ldl12.png")


# others
generate_scatterplot(Table2[,2], Table2[,3], colnames(Table2)[2], paste0("scatterplot_pub", colnames(Table2)[2], ".png"))
generate_scatterplot(Table2[,4], Table2[,5], "WT", "scatterplot_pubWT.png")

for (i in seq(2, ncol(Table2), 2)) {
  print(generate_scatterplot(Table2[, i], Table2[, (i + 1)], colnames(Table2)[i], paste0("scatterplot_", colnames(Table2)[i], ".png")))
}
