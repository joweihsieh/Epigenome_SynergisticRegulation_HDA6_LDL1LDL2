awk '{if ($6 == "+") {printf "%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$2+1,$4,$5,$6} else {printf "%s\t%s\t%s\t%s\t%s\t%s\n",$1,$3,$3+1,$4,$5,$6}}' tair10_gene2.bed > tair10_tss.bed
cat tair10_antagonistic_TE.bed tair10_synergistic_TE_001.bed tair10_synergistic_TE_101.bed |awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,$4}' > tair10_TE_anta_and_Syne.bed
bedtools sort -i tair10_tss.bed > tair10_tss_sorted.bed
bedtools sort -i tair10_TE_anta_and_Syne.bed > tair10_TE_anta_and_Syne_sorted.bed
bedtools closest -a tair10_TE_anta_and_Syne_sorted.bed -b tair10_tss_sorted.bed 
