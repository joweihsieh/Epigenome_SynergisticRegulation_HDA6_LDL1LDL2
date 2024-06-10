#!/usr/bin/perl -w
use strict;
use Data::Dumper;
my $steps = 10000000;
my $list = getList("sample_list.txt");
my $outName = "UnionSite";
my $num_samples = scalar @$list;
open (Ometh,"| gzip -c > $outName"."_meth.txt.gz");
open (OC,"| gzip -c > $outName"."_C.txt.gz");
open (OmC,"| gzip -c > $outName"."_mC.txt.gz");


my @sample_title;
print Dumper($list);
for (my $i = 0; $i < $num_samples; $i++) {
	push(@sample_title,$list->[$i][0]);
	open (my $fh , "zcat $list->[$i][1]|");
	$list->[$i][2]  = $fh;
}

my $titleOut = join("\t",@sample_title);
print Ometh "Chr\tPos\tContext\tDepth\t$titleOut\n";
print OC "Chr\tPos\tContext\tDepth\t$titleOut\n";
print OmC "Chr\tPos\tContext\tDepth\t$titleOut\n";


my $temp = '';
my $oldChr = '';
my $newChr = '';
my $eof_flag = $num_samples;
my $tempHolder = {};
my $majorHolder = {};
my $flush_flag = 0;
my $chr_flag = 0;
my $targetPosition = $steps;
my $t1 = time;
while ($eof_flag) {
	my $t2 = time;
	my $t3 = $t2 - $t1;
	my $largePos = 0;
	print STDERR "Now analyze chr $oldChr at $targetPosition for $t3 sec\n";
	for (my $i = 0; $i < $num_samples; $i++) {
		my $fhi =  $list->[$i][2];
		while (<$fhi>) {
			chomp;
			my @ar = split("\t",$_);
			if ($ar[0] ne $oldChr || $ar[2] > $targetPosition) {
#				print "sample $i pos $ar[2]\n";
				$largePos = ($largePos >  $ar[2]) ? $largePos : $ar[2];
				$newChr = $ar[0];
				$flush_flag = 1;
				$tempHolder->{$ar[2]}->{$i} = \@ar;
#				push ($tempHolder,[(@ar,$i)]);
#				push ($tempHolder,"$_\t$i");
				last;
			} else {
				$majorHolder->{$ar[2]}->{$i} = \@ar;
#				push ($majorHolder,[(@ar,$i)]);
#				push ($majorHolder,"$_\t$i");
			}
			
		}
		$eof_flag-- if (eof($fhi)); 
	}

	if ($flush_flag) {
		$flush_flag = 0;
		doFlush2($majorHolder) if (scalar keys %$majorHolder);
		$majorHolder = $tempHolder;
		$tempHolder = {};
#		$targetPosition = ($oldChr eq $newChr) ? $targetPosition + $steps : $steps;

		$targetPosition = $largePos + $steps;
		$oldChr = $newChr;
#		print "another loop\n";
#		print Dumper($majorHolder);
	}

}
doFlush2($majorHolder);


sub doFlush2 {
	my ($holder) = @_;
	my @posArray = sort {$a <=> $b} keys %$holder;
	foreach my $pos (@posArray) {
		my $chr = '';
		my $context = '';
		my $depth = 1000000;
		my @meth;
		my @C;
		my @mC;
		my $value = $holder->{$pos};
		#	print Dumper($pos,$value);
		for (my $i = 0; $i < $num_samples; $i++) {
			if (exists $value->{$i}) {
				$depth = ($depth < $value->{$i}->[7]) ? $depth : $value->{$i}->[7];
				$chr = $value->{$i}->[0];
				$context = $value->{$i}->[3];
				my $meth_new = sprintf("%.3f",$value->{$i}->[6] / $value->{$i}->[7]);
				push (@meth, $meth_new);
#				push (@meth, $value->{$i}->[5]);
				push (@mC, $value->{$i}->[6]);
				push (@C, $value->{$i}->[7]);
			} else {
				$depth = 0;
				push (@meth, '-');
				push (@mC, '-');
				push (@C, '-');
			}
		}
		#print Dumper($depth,\@meth,\@mC,\@C);
		my $methOut = join("\t",@meth);
		print Ometh"$chr\t$pos\t$context\t$depth\t$methOut\n";
		my $mCOut = join("\t",@mC);
		print OmC"$chr\t$pos\t$context\t$depth\t$mCOut\n";
		my $COut = join("\t",@C);
		print OC "$chr\t$pos\t$context\t$depth\t$COut\n";
	}
}




sub getList {
	my ($in) = @_;
	open (IF, $in);
	my $list = [];
	while (<IF>) {
		chomp;
		my @ar = split("\t",$_);
		next if (scalar @ar < 2);
		push(@$list,\@ar);
	}
	return $list;

}




