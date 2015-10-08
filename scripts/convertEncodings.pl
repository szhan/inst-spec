#! /usr/bin/perl
use strict;
use warnings;

### USAGE
my $USAGE = "USAGE: perl convertEncodings.pl <in.file>\n";
die $USAGE unless scalar(@ARGV) == 1;

### FIELD
my $inFile = $ARGV[0];

### RUN
&parsePloidyFile($inFile);

### SUBROUTINES
# parse ploidy level data, and convert encoding
# ORIGINAL: 0 = diploid; 1 = polyploid
# RECODING: 1 = diploid; 2 = polyploid
# file is tab-delimited with taxon ID and ploidy level
sub parsePloidyFile{
	my $file = $_[0];
	
	open FILE, $file or die "ERROR\tfailed to open $file";
	
	while(<FILE>){
		s/[\n\r]//g;
		my @t = split m/\t/, $_;
		
		my $clade = $t[0];
		my $recoding;
		# check original encoding
		if($t[1] eq '0'){
			$recoding = 1;
		}elsif($t[1] eq '1'){
			$recoding = 2;
		}elsif($t[1] eq 'NA'){
			$recoding = 'NA';
		}else{
			die "ERROR\t$file\t$clade\tWTF is \'$t[1]\'\n";
		}
		
		print "$clade\t$recoding\n";
	}
	
	close FILE;
}

