echo "TEgene Helitron.intact"
bedtools intersect -a tair10_TE.bed -b tair10.fa.mod.Helitron.intact.raw.bed -wa |wc -l
echo "TEelement Helitron.intact"
bedtools intersect -a tair10_Transposon_element_sorted.bed -b tair10.fa.mod.Helitron.intact.raw.bed -wa |wc -l
echo "TEgene TIR.intact"
bedtools intersect -a tair10_TE.bed -b tair10.fa.mod.TIR.intact.raw2.bed -wa |wc -l
echo "TEelement TIR.intact"
bedtools intersect -a tair10_Transposon_element_sorted.bed -b tair10.fa.mod.TIR.intact.raw2.bed -wa |wc -l

bedtools intersect -a tair10_TE_combined.bed -b tair10.fa.mod.TIR.intact.raw2.bed -wo > tair10_TE_combined_overlap_TIR.txt
bedtools intersect -a tair10_TE_combined.bed -b tair10.fa.mod.Helitron.intact.raw2.bed -wo > tair10_TE_combined_overlap_Helitron.txt

