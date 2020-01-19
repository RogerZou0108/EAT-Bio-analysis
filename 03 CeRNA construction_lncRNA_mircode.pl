################################################################################
################################################################################
################################################################################
use strict;
use warnings;

my %hash=();

open(RF,"diff_lncRNA.txt") or die $!;
while(my $line=<RF>){
	chomp($line);
	$hash{$line}=1;
}
close(RF);

my %rep=();
open(RF,"mircode.txt") or die $!;
open(WF,">lncRNA_mircode.txt") or die $!;
while(my $line=<RF>){
	if($.==1){
		print WF "lncRNA\tmiRNA\n";
		next;
	}
	my @arr=split(/\t/,$line);
	my @oneArr=split(/\./,$arr[1]);
	if(exists $hash{$oneArr[0]}){
		my @threeArr=split(/\//,$arr[3]);
		foreach my $mir(@threeArr){
			if($mir=~/^miR/){
				$mir="hsa-$mir";
			}
			else{
				$mir="hsa-miR-$mir";
			}
			my $out="$oneArr[0]\t$mir";
			unless(exists $rep{$out}){
				print WF "$out\n";
				$rep{$out}=1;
			}
		}
	}
}
close(WF);
close(RF);

################################################################################
################################################################################
################################################################################
################################################################################