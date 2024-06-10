#!/usr/bin/perl -w

use strict;
use Data::Dumper;


my @files = `ls *_paired.fastq.gz`;
#print Dumper(@files);
foreach my $file (@files) {
	chomp $file;
	my ($id) = $file =~ /(.+)_paired.fastq.gz/;
	my $cmd = "bs_seeker2-align.py -i $id"."_paired.fastq.gz -o $id.bam -g tair10.fa --aligner=bowtie2 --temp_dir=./ --bt2-p 30 > $id.log &";
	print "$cmd\n";
	system($cmd);
}
