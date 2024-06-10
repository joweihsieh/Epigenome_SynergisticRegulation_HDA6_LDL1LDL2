bs_seeker2-align.py -i test_1_paired.fastq.gz -o test_1.bam -g tair10.fa --aligner=bowtie2 --bt2-p 20 > test_1.log

/work1/home/yenmr/tools/BSseeker2-master/bs_seeker2-call_methylation.py -x -i test_1.bam -d /work1/home/yenmr/tools/BSseeker2-master/bs_utils/reference_genomes/tair10.fa_bowtie2/


