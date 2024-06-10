
#ChIP-seq Trimming
echo "Now Trimming ChIP_H3Ac_WT_1.fastq.gz"
java -jar trimmomatic-0.33.jar SE -phred33 ChIP_H3Ac_WT_1.fastq.gz ChIP_H3Ac_WT_1_trimmed.fastq.gz ILLUMINACLIP:TruSeq2-SE.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

echo "Now Trimming ChIP_H3Ac_WT_2.fastq.gz"
java -jar trimmomatic-0.33.jar SE -phred33 ChIP_H3Ac_WT_2.fastq.gz ChIP_H3Ac_WT_2_trimmed.fastq.gz ILLUMINACLIP:TruSeq2-SE.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

echo "Now Trimming ChIP_H3Ac_hda6_1.fastq.gz"
java -jar trimmomatic-0.33.jar SE -phred33 ChIP_H3Ac_hda6_1.fastq.gz ChIP_H3Ac_hda6_1_trimmed.fastq.gz ILLUMINACLIP:TruSeq2-SE.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

echo "Now Trimming ChIP_H3Ac_hda6_2.fastq.gz"
java -jar trimmomatic-0.33.jar SE -phred33 ChIP_H3Ac_hda6_2.fastq.gz ChIP_H3Ac_hda6_2_trimmed.fastq.gz ILLUMINACLIP:TruSeq2-SE.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

echo "Now Trimming ChIP_H3Ac_hda6ldl12_1.fastq.gz"
java -jar trimmomatic-0.33.jar SE -phred33 ChIP_H3Ac_hda6ldl12_1.fastq.gz ChIP_H3Ac_hda6ldl12_1_trimmed.fastq.gz ILLUMINACLIP:TruSeq2-SE.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

echo "Now Trimming ChIP_H3Ac_hda6ldl12_2.fastq.gz"
java -jar trimmomatic-0.33.jar SE -phred33 ChIP_H3Ac_hda6ldl12_2.fastq.gz ChIP_H3Ac_hda6ldl12_2_trimmed.fastq.gz ILLUMINACLIP:TruSeq2-SE.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

echo "Now Trimming ChIP_H3Ac_ldl12_1.fastq.gz"
java -jar trimmomatic-0.33.jar SE -phred33 ChIP_H3Ac_ldl12_1.fastq.gz ChIP_H3Ac_ldl12_1_trimmed.fastq.gz ILLUMINACLIP:TruSeq2-SE.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

echo "Now Trimming ChIP_H3Ac_ldl12_2.fastq.gz"
java -jar trimmomatic-0.33.jar SE -phred33 ChIP_H3Ac_ldl12_2.fastq.gz ChIP_H3Ac_ldl12_2_trimmed.fastq.gz ILLUMINACLIP:TruSeq2-SE.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36



echo "Now Trimming ChIP_H3K4me2_WT_1.fastq.gz"
java -jar trimmomatic-0.33.jar SE -phred33 ChIP_H3K4me2_WT_1.fastq.gz ChIP_H3K4me2_WT_1_trimmed.fastq.gz ILLUMINACLIP:TruSeq2-SE.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

echo "Now Trimming ChIP_H3K4me2_WT_2.fastq.gz"
java -jar trimmomatic-0.33.jar SE -phred33 ChIP_H3K4me2_WT_2.fastq.gz ChIP_H3K4me2_WT_2_trimmed.fastq.gz ILLUMINACLIP:TruSeq2-SE.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

echo "Now Trimming ChIP_H3K4me2_hda6_1.fastq.gz"
java -jar trimmomatic-0.33.jar SE -phred33 ChIP_H3K4me2_hda6_1.fastq.gz ChIP_H3K4me2_hda6_1_trimmed.fastq.gz ILLUMINACLIP:TruSeq2-SE.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

echo "Now Trimming ChIP_H3K4me2_hda6_2.fastq.gz"
java -jar trimmomatic-0.33.jar SE -phred33 ChIP_H3K4me2_hda6_2.fastq.gz ChIP_H3K4me2_hda6_2_trimmed.fastq.gz ILLUMINACLIP:TruSeq2-SE.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

echo "Now Trimming ChIP_H3K4me2_hda6ldl12_1.fastq.gz"
java -jar trimmomatic-0.33.jar SE -phred33 ChIP_H3K4me2_hda6ldl12_1.fastq.gz ChIP_H3K4me2_hda6ldl12_1_trimmed.fastq.gz ILLUMINACLIP:TruSeq2-SE.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

echo "Now Trimming ChIP_H3K4me2_hda6ldl12_2.fastq.gz"
java -jar trimmomatic-0.33.jar SE -phred33 ChIP_H3K4me2_hda6ldl12_2.fastq.gz ChIP_H3K4me2_hda6ldl12_2_trimmed.fastq.gz ILLUMINACLIP:TruSeq2-SE.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

echo "Now Trimming ChIP_H3K4me2_ldl12_1.fastq.gz"
java -jar trimmomatic-0.33.jar SE -phred33 ChIP_H3K4me2_ldl12_1.fastq.gz ChIP_H3K4me2_ldl12_1_trimmed.fastq.gz ILLUMINACLIP:TruSeq2-SE.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

echo "Now Trimming ChIP_H3K4me2_ldl12_2.fastq.gz"
java -jar trimmomatic-0.33.jar SE -phred33 ChIP_H3K4me2_ldl12_2.fastq.gz ChIP_H3K4me2_ldl12_2_trimmed.fastq.gz ILLUMINACLIP:TruSeq2-SE.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36



echo "Now Trimming ChIP_hda6_1.fastq.gz"
java -jar trimmomatic-0.33.jar SE -phred33 ChIP_hda6_1.fastq.gz ChIP_hda6_1_trimmed.fastq.gz ILLUMINACLIP:TruSeq2-SE.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

echo "Now Trimming ChIP_hda6_2.fastq.gz"
java -jar trimmomatic-0.33.jar SE -phred33 ChIP_hda6_2.fastq.gz ChIP_hda6_2_trimmed.fastq.gz ILLUMINACLIP:TruSeq2-SE.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

echo "Now Trimming ChIP_hda6_input.fastq.gz"
java -jar trimmomatic-0.33.jar SE -phred33 ChIP_hda6_input.fastq.gz ChIP_hda6_input_trimmed.fastq.gz ILLUMINACLIP:TruSeq2-SE.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

echo "Now Trimming ChIP_ldl1_1.fastq.gz"
java -jar trimmomatic-0.33.jar SE -phred33 ChIP_ldl1_1.fastq.gz ChIP_ldl1_1_trimmed.fastq.gz ILLUMINACLIP:TruSeq2-SE.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

echo "Now Trimming ChIP_ldl1_2.fastq.gz"
java -jar trimmomatic-0.33.jar SE -phred33 ChIP_ldl1_2.fastq.gz ChIP_ldl1_2_trimmed.fastq.gz ILLUMINACLIP:TruSeq2-SE.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

echo "Now Trimming ChIP_ldl1_input.fastq.gz"
java -jar trimmomatic-0.33.jar SE -phred33 ChIP_ldl1_input.fastq.gz ChIP_ldl1_input_trimmed.fastq.gz ILLUMINACLIP:TruSeq2-SE.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36


