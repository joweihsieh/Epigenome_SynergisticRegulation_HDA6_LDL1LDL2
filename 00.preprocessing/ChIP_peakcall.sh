macs2 callpeak -c ChIP_hda6_input.bam -t ChIP_hda6_1.bam -f BAM -g 1.3e8 --outdir=peakcall -n ChIP_hda6_1
macs2 callpeak -c ChIP_hda6_input.bam -t ChIP_hda6_2.bam -f BAM -g 1.3e8 --outdir=peakcall -n ChIP_hda6_2
macs2 callpeak -c ChIP_ldl1_input.bam -t ChIP_ldl1_1.bam -f BAM -g 1.3e8 --outdir=peakcall -n ChIP_ldl1_1
macs2 callpeak -c ChIP_ldl1_input.bam -t ChIP_ldl1_2.bam -f BAM -g 1.3e8 --outdir=peakcall -n ChIP_ldl1_2

macs2 callpeak -t ChIP_H3Ac_WT_1.bam -f BAM -g 1.3e8 --outdir=peakcall -n ChIP_H3Ac_WT_1 --nomodel --extsize 147
macs2 callpeak -t ChIP_H3Ac_WT_2.bam -f BAM -g 1.3e8 --outdir=peakcall -n ChIP_H3Ac_WT_2 --nomodel --extsize 147
macs2 callpeak -t ChIP_H3Ac_hda6_1.bam -f BAM -g 1.3e8 --outdir=peakcall -n ChIP_H3Ac_hda6_1 --nomodel --extsize 147
macs2 callpeak -t ChIP_H3Ac_hda6_2.bam -f BAM -g 1.3e8 --outdir=peakcall -n ChIP_H3Ac_hda6_2 --nomodel --extsize 147
macs2 callpeak -t ChIP_H3Ac_hda6ldl12_1.bam -f BAM -g 1.3e8 --outdir=peakcall -n ChIP_H3Ac_hda6ldl12_1 --nomodel --extsize 147
macs2 callpeak -t ChIP_H3Ac_hda6ldl12_2.bam -f BAM -g 1.3e8 --outdir=peakcall -n ChIP_H3Ac_hda6ldl12_2 --nomodel --extsize 147
macs2 callpeak -t ChIP_H3Ac_ldl12_1.bam -f BAM -g 1.3e8 --outdir=peakcall -n ChIP_H3Ac_ldl12_1 --nomodel --extsize 147
macs2 callpeak -t ChIP_H3Ac_ldl12_2.bam -f BAM -g 1.3e8 --outdir=peakcall -n ChIP_H3Ac_ldl12_2 --nomodel --extsize 147

macs2 callpeak -t ChIP_H3K4me2_WT_1.bam -f BAM -g 1.3e8 --outdir=peakcall -n ChIP_H3K4me2_WT_1 --nomodel --extsize 147
macs2 callpeak -t ChIP_H3K4me2_WT_2.bam -f BAM -g 1.3e8 --outdir=peakcall -n ChIP_H3K4me2_WT_2  --nomodel --extsize 147
macs2 callpeak -t ChIP_H3K4me2_hda6_1.bam -f BAM -g 1.3e8 --outdir=peakcall -n ChIP_H3K4me2_hda6_1 --nomodel --extsize 147
macs2 callpeak -t ChIP_H3K4me2_hda6_2.bam -f BAM -g 1.3e8 --outdir=peakcall -n ChIP_H3K4me2_hda6_2 --nomodel --extsize 147
macs2 callpeak -t ChIP_H3K4me2_hda6ldl12_1.bam -f BAM -g 1.3e8 --outdir=peakcall -n ChIP_H3K4me2_hda6ldl12_1 --nomodel --extsize 147
macs2 callpeak -t ChIP_H3K4me2_hda6ldl12_2.bam -f BAM -g 1.3e8 --outdir=peakcall -n ChIP_H3K4me2_hda6ldl12_2 --nomodel --extsize 147
macs2 callpeak -t ChIP_H3K4me2_ldl12_1.bam -f BAM -g 1.3e8 --outdir=peakcall -n ChIP_H3K4me2_ldl12_1 --nomodel --extsize 147
macs2 callpeak -t ChIP_H3K4me2_ldl12_2.bam -f BAM -g 1.3e8 --outdir=peakcall -n ChIP_H3K4me2_ldl12_2 --nomodel --extsize 147

