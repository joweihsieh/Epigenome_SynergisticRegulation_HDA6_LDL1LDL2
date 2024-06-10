#RNA-seq
#Mapping
hisat2 --dta-cufflinks -p 30 -x tair10 -1 C_6_1.fq.gz -2 C_6_2.fq.gz 2>>mapping.log|samtools sort > C_6.bam
hisat2 --dta-cufflinks -p 30 -x tair10 -1 C_7_1.fq.gz -2 C_7_2.fq.gz 2>>mapping.log|samtools sort > C_7.bam
hisat2 --dta-cufflinks -p 30 -x tair10 -1 C_7_2_1.fq.gz -2 C_7_2_2.fq.gz 2>>mapping.log|samtools sort > C_7_2.bam
hisat2 --dta-cufflinks -p 30 -x tair10 -1 L12_6_1.fq.gz -2 L12_6_2.fq.gz 2>>mapping.log|samtools sort > L12_6.bam
hisat2 --dta-cufflinks -p 30 -x tair10 -1 L12_7_1.fq.gz -2 L12_7_2.fq.gz 2>>mapping.log|samtools sort > L12_7.bam
hisat2 --dta-cufflinks -p 30 -x tair10 -1 L12a_6_1.fq.gz -2 L12a_6_2.fq.gz 2>>mapping.log|samtools sort > L12a_6.bam
hisat2 --dta-cufflinks -p 30 -x tair10 -1 L12a_7_1.fq.gz -2 L12a_7_2.fq.gz 2>>mapping.log|samtools sort > L12a_7.bam
hisat2 --dta-cufflinks -p 30 -x tair10 -1 a_6_1.fq.gz -2 a_6_2.fq.gz 2>>mapping.log|samtools sort > a_6.bam
hisat2 --dta-cufflinks -p 30 -x tair10 -1 a_7_1.fq.gz -2 a_7_2.fq.gz 2>>mapping.log|samtools sort > a_7.bam


