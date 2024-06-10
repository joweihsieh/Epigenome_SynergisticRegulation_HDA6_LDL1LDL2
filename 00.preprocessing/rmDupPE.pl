#!/usr/bin/perl -w
use strict;
use Data::Dumper;

my $infile1 = $ARGV[0];
my $infile2 = $ARGV[1];
my $outfile1 = $ARGV[2];
my $outfile2 = $ARGV[3];

open (IF1, "zcat $infile1|");
open (IF2, "zcat $infile2|");
open (OF1, "| gzip > $ARGV[2]");
open (OF2, "| gzip > $ARGV[3]");
my $temp = {};
my $countOri = 0;
my $countKeep = 0;
while (my $line1_1 = <IF1>) {
	$countOri++;
	my $line1_2 = <IF1>;
	my $line1_3 = <IF1>;
	my $line1_4 = <IF1>;
	my $line2_1 = <IF2>;
	my $line2_2 = <IF2>;
	my $line2_3 = <IF2>;
	my $line2_4 = <IF2>;
	next if (exists $temp->{"$line1_2\t$line2_2"});
	$countKeep++;
	$temp->{"$line1_2\t$line2_2"}++;
	print OF1 $line1_1.$line1_2.$line1_3.$line1_4;
	print OF2 $line2_1.$line2_2.$line2_3.$line2_4;
#	print Dumper();
}
print "Original reads: $countOri\nUnique reads: $countKeep\n";

