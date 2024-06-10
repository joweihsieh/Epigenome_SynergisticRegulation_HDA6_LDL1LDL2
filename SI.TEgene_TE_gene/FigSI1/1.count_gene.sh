samtools view C_6n.bam |htseq-count - genes2.gtf -s no > C_6.txt 
samtools view C_7n.bam |htseq-count - genes2.gtf -s no > C_7.txt
samtools view C_7_2n.bam |htseq-count - genes2.gtf -s no > C_7_2.txt
samtools view L12_6n.bam |htseq-count - genes2.gtf -s no > L12_6.txt
samtools view L12_7n.bam |htseq-count - genes2.gtf -s no > L12_7.txt
samtools view L12a_6n.bam |htseq-count - genes2.gtf -s no > L12a_6.txt
samtools view L12a_7n.bam |htseq-count - genes2.gtf -s no > L12a_7.txt
samtools view a_6n.bam |htseq-count - genes2.gtf -s no > a_6.txt
samtools view a_7n.bam |htseq-count - genes2.gtf -s no > a_7.txt

