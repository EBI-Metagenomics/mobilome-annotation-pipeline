#!/usr/bin/perl
use strict;
use warnings;

######## THIS PROGRAM SPLITS A GBK FILE INTO INDEPENDENT FILES ONE PER CONTIG #######
######## Alejandra Escobar-Zepeda, EMBL-EBI/Sanger
######## May 1st, 2019

scalar @ARGV == 1 || die "usage: $0 <file.gbk>
";

my$file=$ARGV[0];

my$name_pref=$file;

my$count=0;
my$flag=0;

## Reading the input file
open OUT_LIST, '>input.list' or die $!;
open (GBK, $file) or die ("Not able to open $file\n") ;
while (<GBK>) {
	chomp;
	if($_ =~ /^LOCUS/){
		$flag=0;
		my$contig_len=(split(/ {1,}/, $_))[2];
		my$contig_name=(split(/ {1,}/, $_))[1];
		if (int($contig_len) >= 5000){
			$flag=1;
			open OUT, ">$contig_name.gbk" or die $!;
			print OUT "$_\n";
			print OUT_LIST "$contig_name.gbk\n";
		}
	}elsif($_ =~ /^\/\//){
		if ($flag==1){
			print OUT "$_\n";
			close(OUT);
		}
	}else{
		if ($flag==1){
			print OUT "$_\n";
		}
	}
}
close(GBK);


