#!/usr/bin/perl -w

use strict;
use Data::Dumper;


my @files = `ls *_1.bam`;
#print Dumper(@files);
foreach my $file (@files) {
	chomp $file;
	my ($id) = $file =~ /(.+)_1.bam/;
	my $cmd = "/work1/home/yenmr/tools/BSseeker2-master/bs_seeker2-call_methylation.py -x -i $id.bam -d /work1/home/yenmr/tools/BSseeker2-master/bs_utils/reference_genomes/tair10.fa_bowtie2/ &";
	#	my $cmd = "bs_seeker2-align.py -i $id"."_paired.fastq.gz -o $id.bam -g tair10.fa --aligner=bowtie2 --temp_dir=./ --bt2-p 30 > $id.log &";
	print "$cmd\n";
	#	system($cmd);
}
