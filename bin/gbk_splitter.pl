#!/usr/bin/perl
# Copyright 2024 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

use strict;
use warnings;

# FIXME: port this to python or a third party tool in bioconda

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


