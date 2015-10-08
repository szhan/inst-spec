#! /usr/bin/perl
use strict;
use warnings;
use List::Util qw(shuffle);

# show usage
my $USAGE = "USAGE: perl makeClasseCmds.pl <in.file> <div.file>\n";
die $USAGE unless scalar(@ARGV) == 2;

# get files
my $inFile = $ARGV[0];
my $divFile = $ARGV[1];

# set parameters
my $num_runs = 100;

# set exceptions
my %exceptions;
$exceptions{'Achillea'} = 155;
$exceptions{'Arisaema'} = 164;
$exceptions{'Cerastium'} = 176;
$exceptions{'Erodium'} = 177;

# set paths to executable/code
my $homeDir = './';
my $R_bin = "R";
my $script = $homeDir.'/scripts/Classe_MCMC_ML_Script.R';

# get diversity data
my @genera = @{ &parseListFile($inFile) };

# go through genera
foreach my $ds (@genera){
	# get trees at random
	#my $num_trees = 189;
	my $num_trees = 1502;
	
	$num_trees = $exceptions{$ds} if exists $exceptions{$ds};
	my @rand_trees = (shuffle(1 .. $num_trees))[0 .. ($num_runs - 1)];
	
	# input files
	#my $treeFile = $homeDir.'/data/mayrose/dryad_jun13/'.$ds.'/'.$ds.'.mb_trees';
	#my $stateFile = $homeDir.'/data/mayrose/ploidyData/classe/'.$ds.'.classe';
	my $treeFile = $homeDir.'/data/ploidb/ChromEvol_Results/'.$ds.'/parsemb_trees.tre';
	my $stateFile = $homeDir.'/data/ploidb/ChromEvol_Results/'.$ds.'/'.$ds.'.classe';
	
	# output dir
	#my $outDir = $homeDir.'/analysis/mayrose/'.$ds.'/';
	my $outDir = $homeDir.'/analysis/ploidb/'.$ds.'/';
	
	foreach my $treeIndex (@rand_trees){
		# output files
		my $tmpDir = $outDir.'/'.$treeIndex.'/';
		my $resFile = $tmpDir.'/'.$ds.'.classe_results';
		my $outFile = $tmpDir.'/'.$ds.'.classe_Rout';
		
		# print command
		print "mkdir -p $outDir; mkdir $tmpDir; $R_bin CMD BATCH --no-save --no-restore \'--args work.dir=\"$tmpDir\" genus=\"$ds\" tree.idx=$treeIndex tree.file=\"$treeFile\" state.file=\"$stateFile\" count.file=\"$divFile\" res.file=\"$resFile\"\' $script $outFile\n";
	}
}

# parse file with list of genera
sub parseListFile{
	my $file = $_[0];
	my @list = ();
	
	open FILE, $file or die "ERROR: failed to open $file";
	while(<FILE>){
		chomp;
		my @t = split m/\t/, $_;
		push @list, $t[0];
	}
	close FILE;
	
	return \@list;
}

# parse diversity data
#original_name	SI	genus
#Anogramma	NA	Anogramma
#Argyrochosma	20	Argyrochosma
sub parseDiversityFile{
	my $file = $_[0];
	my %data = ();
	
	open FILE, $file or die "ERROR: failed to open $file";
	while(<FILE>){
		chomp;
		my @t = split m/\t/, $_;
		$data{$t[2]} = $t[1];
	}
	close FILE;
	
	return \%data;
}

