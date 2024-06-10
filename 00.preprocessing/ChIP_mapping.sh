echo "Now mapping ChIP_H3Ac_WT_1.fastq.gz" >> mapping.log
bowtie2 -p 21 -x genome -U ChIP_H3Ac_WT_1_trimmed.fastq.gz 2>>mapping.log| samtools sort - > ChIP_H3Ac_WT_1.bam
echo "" >> mapping.log

echo "Now mapping ChIP_H3Ac_WT_2.fastq.gz" >> mapping.log
bowtie2 -p 21 -x genome -U ChIP_H3Ac_WT_2_trimmed.fastq.gz 2>>mapping.log| samtools sort - > ChIP_H3Ac_WT_2.bam
echo "" >> mapping.log

echo "Now mapping ChIP_H3Ac_hda6_1.fastq.gz" >> mapping.log
bowtie2 -p 21 -x genome -U ChIP_H3Ac_hda6_1_trimmed.fastq.gz 2>>mapping.log| samtools sort - > ChIP_H3Ac_hda6_1.bam
echo "" >> mapping.log

echo "Now mapping ChIP_H3Ac_hda6_2.fastq.gz" >> mapping.log
bowtie2 -p 21 -x genome -U ChIP_H3Ac_hda6_2_trimmed.fastq.gz 2>>mapping.log| samtools sort - > ChIP_H3Ac_hda6_2.bam
echo "" >> mapping.log

echo "Now mapping ChIP_H3Ac_hda6ldl12_1.fastq.gz" >> mapping.log
bowtie2 -p 21 -x genome -U ChIP_H3Ac_hda6ldl12_1_trimmed.fastq.gz 2>>mapping.log| samtools sort - > ChIP_H3Ac_hda6ldl12_1.bam
echo "" >> mapping.log

echo "Now mapping ChIP_H3Ac_hda6ldl12_2.fastq.gz" >> mapping.log
bowtie2 -p 21 -x genome -U ChIP_H3Ac_hda6ldl12_2_trimmed.fastq.gz 2>>mapping.log| samtools sort - > ChIP_H3Ac_hda6ldl12_2.bam
echo "" >> mapping.log

echo "Now mapping ChIP_H3Ac_ldl12_1.fastq.gz" >> mapping.log
bowtie2 -p 21 -x genome -U ChIP_H3Ac_ldl12_1_trimmed.fastq.gz 2>>mapping.log| samtools sort - > ChIP_H3Ac_ldl12_1.bam
echo "" >> mapping.log

echo "Now mapping ChIP_H3Ac_ldl12_2.fastq.gz" >> mapping.log
bowtie2 -p 21 -x genome -U ChIP_H3Ac_ldl12_2_trimmed.fastq.gz 2>>mapping.log| samtools sort - > ChIP_H3Ac_ldl12_2.bam
echo "" >> mapping.log




echo "Now mapping ChIP_H3K4me2_WT_1.fastq.gz" >> mapping.log
bowtie2 -p 21 -x genome -U ChIP_H3K4me2_WT_1_trimmed.fastq.gz 2>>mapping.log| samtools sort - > ChIP_H3K4me2_WT_1.bam
echo "" >> mapping.log

echo "Now mapping ChIP_H3K4me2_WT_2.fastq.gz" >> mapping.log
bowtie2 -p 21 -x genome -U ChIP_H3K4me2_WT_2_trimmed.fastq.gz 2>>mapping.log| samtools sort - > ChIP_H3K4me2_WT_2.bam
echo "" >> mapping.log

echo "Now mapping ChIP_H3K4me2_hda6_1.fastq.gz" >> mapping.log
bowtie2 -p 21 -x genome -U ChIP_H3K4me2_hda6_1_trimmed.fastq.gz 2>>mapping.log| samtools sort - > ChIP_H3K4me2_hda6_1.bam
echo "" >> mapping.log

echo "Now mapping ChIP_H3K4me2_hda6_2.fastq.gz" >> mapping.log
bowtie2 -p 21 -x genome -U ChIP_H3K4me2_hda6_2_trimmed.fastq.gz 2>>mapping.log| samtools sort - > ChIP_H3K4me2_hda6_2.bam
echo "" >> mapping.log

echo "Now mapping ChIP_H3K4me2_hda6ldl12_1.fastq.gz" >> mapping.log
bowtie2 -p 21 -x genome -U ChIP_H3K4me2_hda6ldl12_1_trimmed.fastq.gz 2>>mapping.log| samtools sort - > ChIP_H3K4me2_hda6ldl12_1.bam
echo "" >> mapping.log

echo "Now mapping ChIP_H3K4me2_hda6ldl12_2.fastq.gz" >> mapping.log
bowtie2 -p 21 -x genome -U ChIP_H3K4me2_hda6ldl12_2_trimmed.fastq.gz 2>>mapping.log| samtools sort - > ChIP_H3K4me2_hda6ldl12_2.bam
echo "" >> mapping.log

echo "Now mapping ChIP_H3K4me2_ldl12_1.fastq.gz" >> mapping.log
bowtie2 -p 21 -x genome -U ChIP_H3K4me2_ldl12_1_trimmed.fastq.gz 2>>mapping.log| samtools sort - > ChIP_H3K4me2_ldl12_1.bam
echo "" >> mapping.log

echo "Now mapping ChIP_H3K4me2_ldl12_2.fastq.gz" >> mapping.log
bowtie2 -p 21 -x genome -U ChIP_H3K4me2_ldl12_2_trimmed.fastq.gz 2>>mapping.log| samtools sort - > ChIP_H3K4me2_ldl12_2.bam
echo "" >> mapping.log


echo "Now mapping ChIP_hda6_1.fastq.gz" >> mapping.log
bowtie2 -p 21 -x genome -U ChIP_hda6_1_trimmed.fastq.gz 2>>mapping.log| samtools sort - > ChIP_hda6_1.bam
echo "" >> mapping.log

echo "Now mapping ChIP_hda6_2.fastq.gz" >> mapping.log
bowtie2 -p 21 -x genome -U ChIP_hda6_2_trimmed.fastq.gz 2>>mapping.log| samtools sort - > ChIP_hda6_2.bam
echo "" >> mapping.log

echo "Now mapping ChIP_hda6_input.fastq.gz" >> mapping.log
bowtie2 -p 21 -x genome -U ChIP_hda6_input_trimmed.fastq.gz 2>>mapping.log| samtools sort - > ChIP_hda6_input.bam
echo "" >> mapping.log

echo "Now mapping ChIP_ldl1_1.fastq.gz" >> mapping.log
bowtie2 -p 21 -x genome -U ChIP_ldl1_1_trimmed.fastq.gz 2>>mapping.log| samtools sort - > ChIP_ldl1_1.bam
echo "" >> mapping.log

echo "Now mapping ChIP_ldl1_2.fastq.gz" >> mapping.log
bowtie2 -p 21 -x genome -U ChIP_ldl1_2_trimmed.fastq.gz 2>>mapping.log| samtools sort - > ChIP_ldl1_2.bam
echo "" >> mapping.log

echo "Now mapping ChIP_ldl1_input.fastq.gz" >> mapping.log
bowtie2 -p 21 -x genome -U ChIP_ldl1_input_trimmed.fastq.gz 2>>mapping.log| samtools sort - > ChIP_ldl1_input.bam
echo "" >> mapping.log

