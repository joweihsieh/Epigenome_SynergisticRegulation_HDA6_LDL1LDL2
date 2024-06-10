bedtools intersect -a tair10_Transposon_element_sorted.bed -b ChIP_H3Ac_WT_1_peaks.narrowPeak -wa -u > ChIP_H3Ac_WT_1_TEelement.bed
bedtools intersect -a tair10_Transposon_element_sorted.bed -b ChIP_H3Ac_WT_2_peaks.narrowPeak -wa -u > ChIP_H3Ac_WT_2_TEelement.bed
bedtools intersect -a tair10_Transposon_element_sorted.bed -b ChIP_H3Ac_hda6_1_peaks.narrowPeak -wa -u > ChIP_H3Ac_hda6_1_TEelement.bed
bedtools intersect -a tair10_Transposon_element_sorted.bed -b ChIP_H3Ac_hda6_2_peaks.narrowPeak -wa -u > ChIP_H3Ac_hda6_2_TEelement.bed
bedtools intersect -a tair10_Transposon_element_sorted.bed -b ChIP_H3Ac_hda6ldl12_1_peaks.narrowPeak -wa -u > ChIP_H3Ac_hda6ldl12_1_TEelement.bed
bedtools intersect -a tair10_Transposon_element_sorted.bed -b ChIP_H3Ac_hda6ldl12_2_peaks.narrowPeak -wa -u > ChIP_H3Ac_hda6ldl12_2_TEelement.bed
bedtools intersect -a tair10_Transposon_element_sorted.bed -b ChIP_H3Ac_ldl12_1_peaks.narrowPeak -wa -u > ChIP_H3Ac_ldl12_1_TEelement.bed
bedtools intersect -a tair10_Transposon_element_sorted.bed -b ChIP_H3Ac_ldl12_2_peaks.narrowPeak -wa -u > ChIP_H3Ac_ldl12_2_TEelement.bed
bedtools intersect -a tair10_Transposon_element_sorted.bed -b ChIP_H3K4me2_WT_1_peaks.narrowPeak -wa -u > ChIP_H3K4me2_WT_1_TEelement.bed
bedtools intersect -a tair10_Transposon_element_sorted.bed -b ChIP_H3K4me2_WT_2_peaks.narrowPeak -wa -u > ChIP_H3K4me2_WT_2_TEelement.bed
bedtools intersect -a tair10_Transposon_element_sorted.bed -b ChIP_H3K4me2_hda6_1_peaks.narrowPeak -wa -u > ChIP_H3K4me2_hda6_1_TEelement.bed
bedtools intersect -a tair10_Transposon_element_sorted.bed -b ChIP_H3K4me2_hda6_2_peaks.narrowPeak -wa -u > ChIP_H3K4me2_hda6_2_TEelement.bed
bedtools intersect -a tair10_Transposon_element_sorted.bed -b ChIP_H3K4me2_hda6ldl12_1_peaks.narrowPeak -wa -u > ChIP_H3K4me2_hda6ldl12_1_TEelement.bed
bedtools intersect -a tair10_Transposon_element_sorted.bed -b ChIP_H3K4me2_hda6ldl12_2_peaks.narrowPeak -wa -u > ChIP_H3K4me2_hda6ldl12_2_TEelement.bed
bedtools intersect -a tair10_Transposon_element_sorted.bed -b ChIP_H3K4me2_ldl12_1_peaks.narrowPeak -wa -u > ChIP_H3K4me2_ldl12_1_TEelement.bed
bedtools intersect -a tair10_Transposon_element_sorted.bed -b ChIP_H3K4me2_ldl12_2_peaks.narrowPeak -wa -u > ChIP_H3K4me2_ldl12_2_TEelement.bed

bedtools intersect -a tair10_promoter.bed -b ChIP_H3Ac_WT_1_peaks.narrowPeak -wa -u > ChIP_H3Ac_WT_1_promoter.bed
bedtools intersect -a tair10_promoter.bed -b ChIP_H3Ac_WT_2_peaks.narrowPeak -wa -u > ChIP_H3Ac_WT_2_promoter.bed
bedtools intersect -a tair10_promoter.bed -b ChIP_H3Ac_hda6_1_peaks.narrowPeak -wa -u > ChIP_H3Ac_hda6_1_promoter.bed
bedtools intersect -a tair10_promoter.bed -b ChIP_H3Ac_hda6_2_peaks.narrowPeak -wa -u > ChIP_H3Ac_hda6_2_promoter.bed
bedtools intersect -a tair10_promoter.bed -b ChIP_H3Ac_hda6ldl12_1_peaks.narrowPeak -wa -u > ChIP_H3Ac_hda6ldl12_1_promoter.bed
bedtools intersect -a tair10_promoter.bed -b ChIP_H3Ac_hda6ldl12_2_peaks.narrowPeak -wa -u > ChIP_H3Ac_hda6ldl12_2_promoter.bed
bedtools intersect -a tair10_promoter.bed -b ChIP_H3Ac_ldl12_1_peaks.narrowPeak -wa -u > ChIP_H3Ac_ldl12_1_promoter.bed
bedtools intersect -a tair10_promoter.bed -b ChIP_H3Ac_ldl12_2_peaks.narrowPeak -wa -u > ChIP_H3Ac_ldl12_2_promoter.bed
bedtools intersect -a tair10_promoter.bed -b ChIP_H3K4me2_WT_1_peaks.narrowPeak -wa -u > ChIP_H3K4me2_WT_1_promoter.bed
bedtools intersect -a tair10_promoter.bed -b ChIP_H3K4me2_WT_2_peaks.narrowPeak -wa -u > ChIP_H3K4me2_WT_2_promoter.bed
bedtools intersect -a tair10_promoter.bed -b ChIP_H3K4me2_hda6_1_peaks.narrowPeak -wa -u > ChIP_H3K4me2_hda6_1_promoter.bed
bedtools intersect -a tair10_promoter.bed -b ChIP_H3K4me2_hda6_2_peaks.narrowPeak -wa -u > ChIP_H3K4me2_hda6_2_promoter.bed
bedtools intersect -a tair10_promoter.bed -b ChIP_H3K4me2_hda6ldl12_1_peaks.narrowPeak -wa -u > ChIP_H3K4me2_hda6ldl12_1_promoter.bed
bedtools intersect -a tair10_promoter.bed -b ChIP_H3K4me2_hda6ldl12_2_peaks.narrowPeak -wa -u > ChIP_H3K4me2_hda6ldl12_2_promoter.bed
bedtools intersect -a tair10_promoter.bed -b ChIP_H3K4me2_ldl12_1_peaks.narrowPeak -wa -u > ChIP_H3K4me2_ldl12_1_promoter.bed
bedtools intersect -a tair10_promoter.bed -b ChIP_H3K4me2_ldl12_2_peaks.narrowPeak -wa -u > ChIP_H3K4me2_ldl12_2_promoter.bed

