# path:/work1/home/joweihsieh/20190731_HDA6/smRNA/raw_data/OD20200405-CLA204

################
# Step 0 confirm the md5
################
md5sum *.fastq.gz >> md5.txt

################
# Step1 cutadapter and trimming
################

cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC Sample_k1245/k1245_R1.fastq.gz -m 18 -M 26 -q 30 > Sample_k1245/ad-k1245_R1.fastq  2> Sample_k1245/report.txt
cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC Sample_k1246/k1246_R1.fastq.gz -m 18 -M 26 -q 30 > Sample_k1246/ad-k1246_R1.fastq  2> Sample_k1246/report.txt
cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC Sample_k127/k127_R1.fastq.gz -m 18 -M 26 -q 30 > Sample_k127/ad-k127_R1.fastq  2> Sample_k127/report.txt


cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC Sample_k12a4/k12a4_R1.fastq.gz -m 18 -M 26 -q 30 > Sample_k12a4/ad-k12a4_R1.fastq  2> Sample_k12a4/report.txt
cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC Sample_k12a5/k12a5_R1.fastq.gz -m 18 -M 26 -q 30 > Sample_k12a5/ad-k12a5_R1.fastq  2> Sample_k12a5/report.txt
cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC Sample_k12a6/k12a6_R1.fastq.gz -m 18 -M 26 -q 30 > Sample_k12a6/ad-k12a6_R1.fastq  2> Sample_k12a6/report.txt


cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC Sample_ka45/ka45_R1.fastq.gz -m 18 -M 26 -q 30 > Sample_ka45/ad-ka45_R1.fastq  2> Sample_ka45/report.txt
cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC Sample_ka46/ka46_R1.fastq.gz -m 18 -M 26 -q 30 > Sample_ka46/ad-ka46_R1.fastq  2> Sample_ka46/report.txt
cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC Sample_ka7/ka7_R1.fastq.gz -m 18 -M 26 -q 30 > Sample_ka7/ad-ka7_R1.fastq  2> Sample_ka7/report.txt


cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC Sample_kc4/kc4_R1.fastq.gz -m 18 -M 26 -q 30 > Sample_kc4/ad-kc4_R1.fastq  2> Sample_kc4/report.txt
cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC Sample_kc5/kc5_R1.fastq.gz -m 18 -M 26 -q 30 > Sample_kc5/ad-kc5_R1.fastq  2> Sample_kc5/report.txt
cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC Sample_kc6/kc6_R1.fastq.gz -m 18 -M 26 -q 30 > Sample_kc6/ad-kc6_R1.fastq  2> Sample_kc6/report.txt



################
# Step2 FastQC (quality and length distribution)
################

fastqc Sample_k1245/ad-k1245_R1.fastq 
fastqc Sample_k1246/ad-k1246_R1.fastq 
fastqc Sample_k127/ad-k127_R1.fastq 
fastqc Sample_k12a4/ad-k12a4_R1.fastq 
fastqc Sample_k12a5/ad-k12a5_R1.fastq 
fastqc Sample_k12a6/ad-k12a6_R1.fastq 
fastqc Sample_ka45/ad-ka45_R1.fastq 
fastqc Sample_ka46/ad-ka46_R1.fastq 
fastqc Sample_ka7/ad-ka7_R1.fastq 
fastqc Sample_kc4/ad-kc4_R1.fastq 
fastqc Sample_kc5/ad-kc5_R1.fastq 
fastqc Sample_kc6/ad-kc6_R1.fastq 


################
# Step3 mapping with Bowtie2 (BAM)
################

bowtie2-build genome.fa genome

echo "kc4" >> mapping.log
bowtie2 -p 21 -x /work1/home/joweihsieh/20190731_HDA6/smRNA/mapping/ref/genome -U /work1/home/joweihsieh/20190731_HDA6/smRNA/raw_data/OD20200405-CLA204/Sample_kc4/ad-kc4_R1.fastq 2>>mapping.log| /work1/home/yenmr/tools/samtools-1.3.1/samtools sort - > ad-kc4.bam
echo "" >> mapping.log

echo "kc5" >> mapping.log
bowtie2 -p 21 -x /work1/home/joweihsieh/20190731_HDA6/smRNA/mapping/ref/genome -U /work1/home/joweihsieh/20190731_HDA6/smRNA/raw_data/OD20200405-CLA204/Sample_kc5/ad-kc5_R1.fastq 2>>mapping.log| /work1/home/yenmr/tools/samtools-1.3.1/samtools sort - > ad-kc5.bam
echo "" >> mapping.log

echo "kc6" >> mapping.log
bowtie2 -p 21 -x /work1/home/joweihsieh/20190731_HDA6/smRNA/mapping/ref/genome -U /work1/home/joweihsieh/20190731_HDA6/smRNA/raw_data/OD20200405-CLA204/Sample_kc6/ad-kc6_R1.fastq 2>>mapping.log| /work1/home/yenmr/tools/samtools-1.3.1/samtools sort - > ad-kc6.bam
echo "" >> mapping.log

echo "ka45" >> mapping.log
bowtie2 -p 21 -x /work1/home/joweihsieh/20190731_HDA6/smRNA/mapping/ref/genome -U /work1/home/joweihsieh/20190731_HDA6/smRNA/raw_data/OD20200405-CLA204/Sample_ka45/ad-ka45_R1.fastq 2>>mapping.log| /work1/home/yenmr/tools/samtools-1.3.1/samtools sort - > ad-ka45.bam
echo "" >> mapping.log

echo "ka46" >> mapping.log
bowtie2 -p 21 -x /work1/home/joweihsieh/20190731_HDA6/smRNA/mapping/ref/genome -U /work1/home/joweihsieh/20190731_HDA6/smRNA/raw_data/OD20200405-CLA204/Sample_ka46/ad-ka46_R1.fastq 2>>mapping.log| /work1/home/yenmr/tools/samtools-1.3.1/samtools sort - > ad-ka46.bam
echo "" >> mapping.log

echo "ka7" >> mapping.log
bowtie2 -p 21 -x /work1/home/joweihsieh/20190731_HDA6/smRNA/mapping/ref/genome -U /work1/home/joweihsieh/20190731_HDA6/smRNA/raw_data/OD20200405-CLA204/Sample_ka7/ad-ka7_R1.fastq 2>>mapping.log| /work1/home/yenmr/tools/samtools-1.3.1/samtools sort - > ad-ka7.bam
echo "" >> mapping.log

echo "k1245" >> mapping.log
bowtie2 -p 21 -x /work1/home/joweihsieh/20190731_HDA6/smRNA/mapping/ref/genome -U /work1/home/joweihsieh/20190731_HDA6/smRNA/raw_data/OD20200405-CLA204/Sample_k1245/ad-k1245_R1.fastq 2>>mapping.log| /work1/home/yenmr/tools/samtools-1.3.1/samtools sort - > ad-k1245.bam
echo "" >> mapping.log

echo "k1246" >> mapping.log
bowtie2 -p 21 -x /work1/home/joweihsieh/20190731_HDA6/smRNA/mapping/ref/genome -U /work1/home/joweihsieh/20190731_HDA6/smRNA/raw_data/OD20200405-CLA204/Sample_k1246/ad-k1246_R1.fastq 2>>mapping.log| /work1/home/yenmr/tools/samtools-1.3.1/samtools sort - > ad-k1246.bam
echo "" >> mapping.log

echo "k127" >> mapping.log
bowtie2 -p 21 -x /work1/home/joweihsieh/20190731_HDA6/smRNA/mapping/ref/genome -U /work1/home/joweihsieh/20190731_HDA6/smRNA/raw_data/OD20200405-CLA204/Sample_k127/ad-k127_R1.fastq 2>>mapping.log| /work1/home/yenmr/tools/samtools-1.3.1/samtools sort - > ad-k127.bam
echo "" >> mapping.log

echo "k12a4" >> mapping.log
bowtie2 -p 21 -x /work1/home/joweihsieh/20190731_HDA6/smRNA/mapping/ref/genome -U /work1/home/joweihsieh/20190731_HDA6/smRNA/raw_data/OD20200405-CLA204/Sample_k12a4/ad-k12a4_R1.fastq 2>>mapping.log| /work1/home/yenmr/tools/samtools-1.3.1/samtools sort - > ad-k12a4.bam
echo "" >> mapping.log

echo "k12a5" >> mapping.log
bowtie2 -p 21 -x /work1/home/joweihsieh/20190731_HDA6/smRNA/mapping/ref/genome -U /work1/home/joweihsieh/20190731_HDA6/smRNA/raw_data/OD20200405-CLA204/Sample_k12a5/ad-k12a5_R1.fastq 2>>mapping.log| /work1/home/yenmr/tools/samtools-1.3.1/samtools sort - > ad-k12a5.bam
echo "" >> mapping.log

echo "k12a6" >> mapping.log
bowtie2 -p 21 -x /work1/home/joweihsieh/20190731_HDA6/smRNA/mapping/ref/genome -U /work1/home/joweihsieh/20190731_HDA6/smRNA/raw_data/OD20200405-CLA204/Sample_k12a6/ad-k12a6_R1.fastq 2>>mapping.log| /work1/home/yenmr/tools/samtools-1.3.1/samtools sort - > ad-k12a6.bam
echo "" >> mapping.log


################
# Step3.1 select perfect match reads
################

samtools view ad-kc4.bam |grep AS:i:0 > ad-kc4.perfect.txt
samtools view ad-kc5.bam |grep AS:i:0 > ad-kc5.perfect.txt
samtools view ad-kc6.bam |grep AS:i:0 > ad-kc6.perfect.txt
samtools view ad-ka45.bam |grep AS:i:0 > ad-ka45.perfect.txt
samtools view ad-ka46.bam |grep AS:i:0 > ad-ka46.perfect.txt
samtools view ad-ka7.bam |grep AS:i:0 > ad-ka7.perfect.txt
samtools view ad-k1245.bam |grep AS:i:0 > ad-k1245.perfect.txt
samtools view ad-k1246.bam |grep AS:i:0 > ad-k1246.perfect.txt
samtools view ad-k127.bam |grep AS:i:0 > ad-k127.perfect.txt
samtools view ad-k12a4.bam |grep AS:i:0 > ad-k12a4.perfect.txt
samtools view ad-k12a5.bam |grep AS:i:0 > ad-k12a5.perfect.txt
samtools view ad-k12a6.bam |grep AS:i:0 > ad-k12a6.perfect.txt


samtools view -H ad-kc4.bam > ad-kc4.header.txt
samtools view -H ad-kc5.bam > ad-kc5.header.txt
samtools view -H ad-kc6.bam > ad-kc6.header.txt
samtools view -H ad-ka45.bam > ad-ka45.header.txt
samtools view -H ad-ka46.bam > ad-ka46.header.txt
samtools view -H ad-ka7.bam > ad-ka7.header.txt
samtools view -H ad-k1245.bam > ad-k1245.header.txt
samtools view -H ad-k1246.bam > ad-k1246.header.txt
samtools view -H ad-k127.bam > ad-k127.header.txt
samtools view -H ad-k12a4.bam > ad-k12a4.header.txt
samtools view -H ad-k12a5.bam > ad-k12a5.header.txt
samtools view -H ad-k12a6.bam > ad-k12a6.header.txt




cat ad-kc4.header.txt ad-kc4.perfect.txt | samtools view -S -b > ad-kc4.perfect.bam
cat ad-kc5.header.txt ad-kc5.perfect.txt | samtools view -S -b > ad-kc5.perfect.bam
cat ad-kc6.header.txt ad-kc6.perfect.txt | samtools view -S -b > ad-kc6.perfect.bam
cat ad-ka45.header.txt ad-ka45.perfect.txt | samtools view -S -b > ad-ka45.perfect.bam
cat ad-ka46.header.txt ad-ka46.perfect.txt | samtools view -S -b > ad-ka46.perfect.bam
cat ad-ka7.header.txt ad-ka7.perfect.txt | samtools view -S -b > ad-ka7.perfect.bam
cat ad-k1245.header.txt ad-k1245.perfect.txt | samtools view -S -b > ad-k1245.perfect.bam
cat ad-k1246.header.txt ad-k1246.perfect.txt | samtools view -S -b > ad-k1246.perfect.bam
cat ad-k127.header.txt ad-k127.perfect.txt | samtools view -S -b > ad-k127.perfect.bam
cat ad-k12a4.header.txt ad-k12a4.perfect.txt | samtools view -S -b > ad-k12a4.perfect.bam
cat ad-k12a5.header.txt ad-k12a5.perfect.txt | samtools view -S -b > ad-k12a5.perfect.bam
cat ad-k12a6.header.txt ad-k12a6.perfect.txt | samtools view -S -b > ad-k12a6.perfect.bam

rm *.txt


samtools index ad-kc4.perfect.bam
samtools index ad-kc5.perfect.bam
samtools index ad-kc6.perfect.bam
samtools index ad-ka45.perfect.bam
samtools index ad-ka46.perfect.bam
samtools index ad-ka7.perfect.bam
samtools index ad-k1245.perfect.bam
samtools index ad-k1246.perfect.bam
samtools index ad-k127.perfect.bam
samtools index ad-k12a4.perfect.bam
samtools index ad-k12a5.perfect.bam
samtools index ad-k12a6.perfect.bam

################
# Step4 remove rRNA, mitochondrial RNA, tRNA, miRNA and snRNA 
################
#BLAST (v 2.2.25) against the small RNAs in the GenBank database (Release 209.0) and Rfam database (v 11.0) to remove rRNAs, snRNAs, snoRNAs, scRNAs, and tRNAs.
grep gene TAIR10_GFF3_genes_transposons.gff | grep -E 'tRNA|rRNA|snoRNA|snRNA|miRNA' | cut -d ';' -f 1,2 |sed -e 's/;/\t/g' | sed -e 's/ID=//g' | sed -e 's/Note=//g'  > TAIR10_GFF3_smRNA_wo_siRNA.txt
cat TAIR10_GFF3_smRNA_wo_siRNA.txt | awk '{print $1 "\t" $4 "\t" $5 "\t" $9 "\t" $8 "\t" $7}'> TAIR10_GFF3_smRNA_wo_siRNA.bed


################
# Step5 generate tag-count file (sequence with abundance)
################

#conda create -n fasta
#source activate fasta
#mkdir fastx_bin
#cd fastx_bin
#wget http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
#tar -xjf fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2

files=(ad-kc4 ad-kc5 ad-kc6 ad-ka45 ad-ka46 ad-ka7 ad-k1245 ad-k1246 ad-k127 ad-k12a4 ad-k12a5 ad-k12a6)

for j in ${files[@]}
do      

        bedtools bamtobed -i $j.perfect.bam > $j.perfect.bed
        bedtools intersect -a $j.perfect.bed -b ./ref/TAIR10_GFF3_smRNA_wo_siRNA.bed -v -wa -wb > $j.perfect.siRNA.bed
        bedtools getfasta -fi ./ref/genome.fa -bed $j.perfect.siRNA.bed -name > $j.perfect.siRNA.fa
        fastx_collapser -i $j.perfect.siRNA.fa -o $j.perfect.siRNA.tag

done

################
# Step6 BAM to bw
################

bamCoverage --bam ad-kc4.perfect.bam -o kc4.perfect.bw
bamCoverage --bam ad-kc5.perfect.bam -o kc5.perfect.bw
bamCoverage --bam ad-kc6.perfect.bam -o kc6.perfect.bw
bamCoverage --bam ad-ka45.perfect.bam -o ka45.perfect.bw
bamCoverage --bam ad-ka46.perfect.bam -o ka46.perfect.bw
bamCoverage --bam ad-ka7.perfect.bam -o ka7.perfect.bw
bamCoverage --bam ad-k1245.perfect.bam -o k1245.perfect.bw
bamCoverage --bam ad-k1246.perfect.bam -o k1246.perfect.bw
bamCoverage --bam ad-k127.perfect.bam -o k127.perfect.bw
bamCoverage --bam ad-k12a4.perfect.bam -o k12a4.perfect.bw
bamCoverage --bam ad-k12a5.perfect.bam -o k12a5.perfect.bw
bamCoverage --bam ad-k12a6.perfect.bam -o k12a6.perfect.bw


