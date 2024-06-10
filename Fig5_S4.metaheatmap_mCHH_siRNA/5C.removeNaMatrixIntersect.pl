#!/usr/bin/perl -w
use strict;
use Data::Dumper;

my $list_file = 'Intersect_TE.txt';
open (IL, $list_file);
my $list = {};
while (<IL>) {
	chomp;
	$list->{$_}++;
}

#my $infile = "TE_CHH_WT.matrix.gz";
my $infile = $ARGV[0];
open (IF, "zcat $infile|");
my $head = <IF>;
my ($bins) = $head =~ /\"sample_boundaries\":\[.+,(\d+)\]/;
my ($samples) = $head =~ /\"group_boundaries\":\[0,(\d+)\]/;
my $lim = $bins / 2;
#print "Bin $bins\n";

my $holder = [];

while (<IF>) {
#	chomp;
	my $na_count = () = $_ =~ /nan/g;
	my @ar = split("\t",$_);
#	print "$na_count\n";
#	push (@$holder, $_) if ($na_count < $lim);
	push (@$holder, $_) if (exists $list->{$ar[3]});
}
my $size = scalar keys %$list;
#my $size = scalar @$holder;
my $newHead = $head;
$newHead =~ s/$samples/$size/;
print $newHead;
foreach my $line (@$holder) {
	print $line;
}
#print "$size\n";	


