# histone modifications
~/bin/bigWigAverageOverBed ChIP_H3Ac_hda6_1.bw upTE_union_ID.bed ChIP_H3Ac_hda6_1_upTE_union.tab
~/bin/bigWigAverageOverBed ChIP_H3Ac_hda6_2.bw upTE_union_ID.bed ChIP_H3Ac_hda6_2_upTE_union.tab

~/bin/bigWigAverageOverBed ChIP_H3Ac_ldl12_1.bw upTE_union_ID.bed ChIP_H3Ac_ldl12_1_upTE_union.tab
~/bin/bigWigAverageOverBed ChIP_H3Ac_ldl12_2.bw upTE_union_ID.bed ChIP_H3Ac_ldl12_2_upTE_union.tab

~/bin/bigWigAverageOverBed ChIP_H3Ac_hda6ldl12_1.bw upTE_union_ID.bed ChIP_H3Ac_hda6ldl12_1_upTE_union.tab
~/bin/bigWigAverageOverBed ChIP_H3Ac_hda6ldl12_2.bw upTE_union_ID.bed ChIP_H3Ac_hda6ldl12_2_upTE_union.tab

~/bin/bigWigAverageOverBed ChIP_H3Ac_WT_1.bw upTE_union_ID.bed ChIP_H3Ac_WT_1_upTE_union.tab
~/bin/bigWigAverageOverBed ChIP_H3Ac_WT_2.bw upTE_union_ID.bed ChIP_H3Ac_WT_2_upTE_union.tab


~/bin/bigWigAverageOverBed ChIP_H3K4me2_hda6_1.bw upTE_union_ID.bed ChIP_H3K4me2_hda6_1_upTE_union.tab
~/bin/bigWigAverageOverBed ChIP_H3K4me2_hda6_2.bw upTE_union_ID.bed ChIP_H3K4me2_hda6_2_upTE_union.tab

~/bin/bigWigAverageOverBed ChIP_H3K4me2_ldl12_1.bw upTE_union_ID.bed ChIP_H3K4me2_ldl12_1_upTE_union.tab
~/bin/bigWigAverageOverBed ChIP_H3K4me2_ldl12_2.bw upTE_union_ID.bed ChIP_H3K4me2_ldl12_2_upTE_union.tab

~/bin/bigWigAverageOverBed ChIP_H3K4me2_hda6ldl12_1.bw upTE_union_ID.bed ChIP_H3K4me2_hda6ldl12_1_upTE_union.tab
~/bin/bigWigAverageOverBed ChIP_H3K4me2_hda6ldl12_2.bw upTE_union_ID.bed ChIP_H3K4me2_hda6ldl12_2_upTE_union.tab

~/bin/bigWigAverageOverBed ChIP_H3K4me2_WT_1.bw upTE_union_ID.bed ChIP_H3K4me2_WT_1_upTE_union.tab
~/bin/bigWigAverageOverBed ChIP_H3K4me2_WT_2.bw upTE_union_ID.bed ChIP_H3K4me2_WT_2_upTE_union.tab



# binding
bedtools intersect -wa -wb -a ChIP_hda6_peaks.narrowPeak -b upTE_union_ID.bed > ChIP_hda6_peaks.narrowPeak.position
bedtools intersect -wa -wb -a ChIP_ldl1_peaks.narrowPeak -b upTE_union_ID.bed > ChIP_ldl1_peaks.narrowPeak.position




# methylation
bash Run_CG_tab.sh
bash Run_CHG_tab.sh
bash Run_CHH_tab.sh
