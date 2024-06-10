# Function to run bigWigAverageOverBed for a list of files
run_bigWigAverageOverBed <- function(input_files, bed_file, output_files) {
  if (length(input_files) != length(output_files)) {
    stop("The number of input files must match the number of output files")
  }
  
  for (i in seq_along(input_files)) {
    command <- paste("bigWigAverageOverBed", input_files[i], bed_file, output_files[i])
    system(command)
  }
}

# List of input bigWig files
input_files <- c(
  "WT_2_sort_rmMtPt_150bp.bw",
  "WT_3_sort_rmMtPt_150bp.bw",
  "HDA6_1_sort_rmMtPt_150bp.bw",
  "HDA6_4_sort_rmMtPt_150bp.bw",
  "LDL_3_sort_rmMtPt_150bp.bw",
  "LDL_4_sort_rmMtPt_150bp.bw",
  "HL_2_sort_rmMtPt_150bp.bw",
  "HL_3_sort_rmMtPt_150bp.bw"
)

# Corresponding output tab files
output_files <- c(
  "ATAC_WT_1_body.tab",
  "ATAC_WT_2_body.tab",
  "ATAC_hda6_1_body.tab",
  "ATAC_hda6_2_body.tab",
  "ATAC_ldl12_1_body.tab",
  "ATAC_ldl12_2_body.tab",
  "ATAC_hda6ldl12_1_body.tab",
  "ATAC_hda6ldl12_2_body.tab"
)

# BED file
bed_file <- "tair10_TE.bed"

# Run the function
run_bigWigAverageOverBed(input_files, bed_file, output_files)
