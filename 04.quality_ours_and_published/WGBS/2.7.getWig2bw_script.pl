#!/usr/bin/perl -w
use strict;
use Data::Dumper;

my @files = `ls *wig`;

foreach my $file (@files) {
	chomp $file;
	my ($head) = $file =~ /(.+).wig/;
	my $cmd = "wigToBigWig $head.wig ~/genome/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/WholeGenomeFasta/tair10.fa.fai $head.bw";
	print "$cmd\n";
	system($cmd);
}
