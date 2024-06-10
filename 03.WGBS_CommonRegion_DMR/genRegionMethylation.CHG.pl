#!/usr/bin/perl -w
use Data::Dumper;
use strict;
use Statistics::Basic qw(:all);

#This script is to calculate the methylation level of all sites (not just common site) from Union sites file.
my $minCoverage = 4;
my $windowSize = 200;
my $minSites = 5;

#my $contextArray = ["CG","CHG","CHH"];
my $contextArray = ["CHG"];
#my $infile = "UnionSite";
my $infile = $ARGV[0];
my $cFile = "$infile"."_C.txt.gz";
my $mFile = "$infile"."_meth.txt.gz";

open (IC, "zcat $cFile|");
open (IM, "zcat $mFile|");
my $methHolder = {};

my $outfile = {};
foreach my $context (@$contextArray) {
	open ($outfile->{$context}, ">$infile"."_Region_methylation.$context.txt");
}


my $titleLine1 = <IC>;
my $titleLine2 = <IM>;
chomp $titleLine1;
my @titleArray = split ("\t",$titleLine1);
my $arraySize = scalar @titleArray;
my $titleOut = join ("\t",@titleArray[4..$arraySize - 1]);
foreach my $context (@$contextArray) {
	my $outF = $outfile->{$context};
	print $outF "Chr\tRegionStart\tRegionEnd\tSiteStart\tSiteEnd\tContext\t$titleOut\n";
}



my $totalRegionCount = {};
my $posHolder = {};
my $currentChr = '';
my $endPos = $windowSize;
while (my $cLine=<IC>) {
	chomp $cLine;
	my @cArray = split("\t", $cLine);
	my $mLine = <IM>;
	chomp $mLine;
	my @mArray = split("\t", $mLine);
#	print Dumper(\@cArray,\@mArray);
	my $chr = $cArray[0];
	my $pos = $cArray[1];
	my $context = $cArray[2];
	my $mindepth = $cArray[3];
	if ($currentChr ne $chr || $pos > $endPos) {
#		print Dumper($posHolder);
		doflush($methHolder,$chr, $endPos - $windowSize, $endPos, $posHolder);
		$endPos = $pos - ($pos % $windowSize) + $windowSize;
		$currentChr = $chr;
		#	print "$chr\t$endPos\n";
		$methHolder = {};
		$posHolder = {};
	}
	push (@{$posHolder->{$context}}, $pos);
	for (my $i = 4; $i < $arraySize; $i++) {
		my $depth = $cArray[$i];
		my $meth = $mArray[$i];
		push(@{$methHolder->{$context}->{$i}}, $meth) if ($depth ne '-' && $depth >= $minCoverage);
	}

}

print Dumper($totalRegionCount);

sub doflush {
	my ($holder, $chr, $strPos, $endPos, $posHolder) = @_;
	foreach my $context (@$contextArray) {
		my $minSiteInRegion = 10;
		my $allhaveresult = 1;
		my @outArray;
		my $minPos;
		my $maxPos;
		if (exists $holder->{$context}) {
			$totalRegionCount->{$context}++;
			my @sortPos = sort {$a <=> $b} @{$posHolder->{$context}};
			$minPos = $sortPos[0];
			$maxPos = $sortPos[-1];
			#		print Dumper($context,\@sortPos);
			for (my $i = 4; $i < $arraySize; $i++) {
				my $size = 0;
				my $mean = 'NA';
				if (exists $holder->{$context}->{$i}) {
					$size = scalar @{$holder->{$context}->{$i}};
					$mean = mean($holder->{$context}->{$i}) + 0;
					$minSiteInRegion = ($minSiteInRegion < $size) ? $minSiteInRegion : $size;
					push(@outArray, $mean);
				} else {
					$allhaveresult = 0;
				}
			}
		} else {
			$minSiteInRegion = 0;
		}
		my $out = join ("\t", @outArray);
		my $outF = $outfile->{$context};
		#	print $outF "$chr\t$strPos\t$endPos\t$context\t$out\t$minSite\n" if ($minSite > 0 && $allhaveresult);
		print $outF "$chr\t$strPos\t$endPos\t$minPos\t$maxPos\t$context\t$out\n" if ($allhaveresult && $minSiteInRegion >= $minSites);
	}
#	print Dumper($holder->{"CHG"});
}
