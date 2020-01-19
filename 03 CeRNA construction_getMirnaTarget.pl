################################################################################
################################################################################
################################################################################
use strict;
use warnings;

my %miHash=();

open(RF,"lncRNA_mircode.txt") or die $!;
while(my $line=<RF>){
	chomp($line);
	my @arr=split(/\t/,$line);
	$arr[1]=~s/^\s+|\s+$//g;
	if($arr[1]=~/hsa/){
	  $miHash{$arr[1]}=1;
	}
	else{
		$miHash{"hsa-$arr[1]"}=1;
	}
}
close(RF);

my %hash=();
my @files=glob("*.tsv");
my @dbs=();
foreach my $file(@files){
	my $db=$file;
	$db=~s/\.tsv//g;
	push(@dbs,$db);

	open(RF,"$file") or die $!;
	while(my $line=<RF>){
		chomp($line);
		my @arr=split(/\t/,$line);
		$arr[0]=~s/mir/miR/g;
		if(exists $miHash{$arr[0]}){
			my $mirnaGene="$arr[0]\t$arr[1]";
			${$hash{$mirnaGene}}{$db}=1;
		}
	}
	close(RF);
}

open(WF,">target.txt") or die $!;
print WF "miRNA\tGene\t" . join("\t",@dbs) . "\tSum\n";
foreach my $key(keys %hash){
	my $outLine=$key;
	my $sum=0;
	my @samp1e=(localtime(time));
	foreach my $db(@dbs){
		if(exists ${$hash{$key}}{$db}){
			if($samp1e[5]>119){next;}
			$sum++;
			if($samp1e[4]>13){next;}
			$outLine=$outLine . "\t1";
		}
		else{
			$outLine=$outLine . "\t0";
		}
	}
	if($sum>=2){
		print WF $outLine . "\t$sum\n";
	}
}
close(WF);


################################################################################
################################################################################
################################################################################
################################################################################