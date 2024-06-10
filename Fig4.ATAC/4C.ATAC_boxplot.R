# Read ATAC-seq data
WT1 <- read.table("ATAC_WT_1_body.tab")
WT2 <- read.table("ATAC_WT_2_body.tab")
hda61 <- read.table("ATAC_hda6_1_body.tab")
hda62 <- read.table("ATAC_hda6_2_body.tab")
hda6ldl121 <- read.table("ATAC_hda6ldl12_1_body.tab")
hda6ldl122 <- read.table("ATAC_hda6ldl12_2_body.tab")
ldl121 <- read.table("ATAC_ldl12_1_body.tab")
ldl122 <- read.table("ATAC_ldl12_2_body.tab")

# Read RNA expression data and calculate means
RNA <- read.table("Exprssed_TE.txt", header = TRUE)
RNA$RNA_WT <- (RNA$C_6 + RNA$C_7 + RNA$C_7_2) / 3
RNA$RNA_hda6 <- (RNA$a_6 + RNA$a_7) / 2
RNA$RNA_ldl12 <- (RNA$L12_6 + RNA$L12_7) / 2
RNA$RNA_hda6ldl12 <- (RNA$L12a_6 + RNA$L12a_7) / 2

# Calculate ATAC-seq means and merge RNA data
df <- data.frame(
  row.names = WT1$V1,
  WT = (WT1$V6 + WT2$V6) / 2,
  hda6 = (hda61$V6 + hda62$V6) / 2,
  ldl12 = (ldl121$V6 + ldl122$V6) / 2,
  hda6ldl12 = (hda6ldl121$V6 + hda6ldl122$V6) / 2
)

# Merge RNA data
df2 <- merge(df, RNA[, c(1, 24:27)], by.x = 'row.names', by.y = 'Row.names')
row.names(df2) <- df2$Row.names
df <- df2[, 2:9]


# Function to create boxplots
create_boxplot <- function(bed_file, output_image, output_data, ylim_range = c(0, 100), ylab_text = 'ATAC-seq abundance (RPKM)', main_text = 'Boxplot', outline = FALSE) {
  selected <- read.table(bed_file)
  df_filtered <- df[row.names(df) %in% selected$V4, ]
  
  png(file = output_image, width = 1500, height = 2000, res = 350)
  boxplot(df_filtered[, c(1, 3, 2, 4)], ylim = ylim_range, ylab = ylab_text, main = main_text, outline = outline, col = NA)
  dev.off()
  
  write.table(df_filtered[, c(1, 3, 2, 4)], file = output_data, sep = '\t', quote = FALSE)

  print(c(colnames(df_filtered)[1], colnames(df_filtered)[2], t.test(df_filtered[,1], df_filtered[,2], paired = TRUE)$p.value))
  print(c(colnames(df_filtered)[1], colnames(df_filtered)[3], t.test(df_filtered[,1], df_filtered[,3], paired = TRUE)$p.value))
  print(c(colnames(df_filtered)[1], colnames(df_filtered)[4], t.test(df_filtered[,1], df_filtered[,4], paired = TRUE)$p.value))
  #print(c(colnames(df_filtered)[2], colnames(df_filtered)[3], t.test(df_filtered[,2], df_filtered[,3], paired = TRUE)$p.value))
  #print(c(colnames(df_filtered)[2], colnames(df_filtered)[4], t.test(df_filtered[,2], df_filtered[,4], paired = TRUE)$p.value))
  #print(c(colnames(df_filtered)[3], colnames(df_filtered)[4], t.test(df_filtered[,3], df_filtered[,4], paired = TRUE)$p.value))

}

# Call the function to generate boxplots
create_boxplot(
  bed_file = "tair10_synergistic_TE_HDA6_predominant.bed",
  output_image = "Boxplot_ATAC_HDA6_predominant.png",
  output_data = "ATAC_value_HDA6_predominant.txt",
  ylim_range = c(0, 100),
  ylab_text = 'ATAC-seq abundance (RPKM)',
  main_text = 'HDA6_predominant',
  outline = FALSE
)

create_boxplot(
  bed_file = "tair10_synergistic_TE_double_locked.bed",
  output_image = "Boxplot_ATAC_double_locked.png",
  output_data = "ATAC_value_double_locked.txt",
  ylim_range = c(0, 100),
  ylab_text = 'ATAC-seq abundance (RPKM)',
  main_text = 'double_locked',
  outline = FALSE
)
