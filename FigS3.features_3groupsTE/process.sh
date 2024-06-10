bigWigAverageOverBed At-C1.bam.CG.bw tair10_genome_500k_bin.bed At-C1.bam.CG.500k.tab
awk '{split($1,id,"_"); printf "%s\t%s\t%s\t%s\t%s\n",id[1],id[2],id[2] + 500000, $1, $6}' At-C1.bam.CG.500k.tab > At-C1.bam.CG.500k.peri.bed
bedtools merge -i At-C1.bam.CG.500k.peri.bed > At-C1.bam.CG.peri.bed


bedtools intersect -a tair10_antagonistic_TE.bed -b At-C1.bam.CG.peri.bed -wa -u|wc
bedtools intersect -a tair10_synergistic_TE_001.bed -b At-C1.bam.CG.peri.bed -wa -u|wc
bedtools intersect -a tair10_synergistic_TE_101.bed -b At-C1.bam.CG.peri.bed -wa -u|wc
