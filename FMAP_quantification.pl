# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die "ERROR in $0: ", $_[0] };

use Cwd 'abs_path';
use Getopt::Long;
use List::Util qw(min);

(my $fmapPath = abs_path($0)) =~ s/\/[^\/]*$//;
my $database = loadDefaultDatabase();
my $databasePath = "$fmapPath/FMAP_data/$database";

GetOptions('h' => \(my $help = ''),
	'c' => \(my $cpmInsteadOfRPKM = ''),
	'i=f' => \(my $minimumPercentIdentity = 80),
	'l=s' => \(my $proteinLengthFile = "$databasePath.length.txt"),
	'o=s' => \(my $proteinOrthologyFile = ''),
	'd=s' => \(my $orthologyDefinitionFile = ''),
	'w=s' => \(my $readNameWeightFile = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl FMAP_quantification.pl [options] blast_hits1.txt [blast_hits2.txt [...]] > abundance.txt

Options: -h       display this help message
         -c       use CPM values instead of RPKM values
         -i FLOAT minimum percent identity [$minimumPercentIdentity]
         -l FILE  tab-delimited text file with the first column having protein names and the second column having the sequence lengths
         -o FILE  tab-delimited text file with the first column having protein names and the second column having the orthology names
         -d FILE  tab-delimited text file with the first column having orthology names and the second column having the definitions
         -w FILE  tab-delimited text file with the first column having read names and the second column having the weights

EOF
}
$orthologyDefinitionFile = "$fmapPath/FMAP_data/KEGG_orthology.txt" if($proteinLengthFile eq "$databasePath.length.txt" && $orthologyDefinitionFile eq '');

my (@inputFileList) = @ARGV;
my %proteinLengthHash = ();
{
	open(my $reader, $proteinLengthFile);
	while(my $line = <$reader>) {
		chomp($line);
		my ($protein, $length) = split(/\t/, $line);
		$proteinLengthHash{$protein} = $length;
	}
	close($reader);
}
my %proteinOrthologyHash = ();
if($proteinOrthologyFile ne '') {
	open(my $reader, $proteinOrthologyFile);
	while(my $line = <$reader>) {
		chomp($line);
		my ($protein, $orthology) = split(/\t/, $line);
		$proteinOrthologyHash{$protein} = $orthology;
	}
	close($reader);
}
my %orthologyDefinitionHash = ();
if($orthologyDefinitionFile ne '') {
	open(my $reader, $orthologyDefinitionFile);
	while(my $line = <$reader>) {
		chomp($line);
		my ($orthology, $definition) = split(/\t/, $line);
		$orthologyDefinitionHash{$orthology} = $definition;
	}
	close($reader);
}
my %readNameWeightHash = ();
if($readNameWeightFile ne '') {
	open(my $reader, $readNameWeightFile);
	while(my $line = <$reader>) {
		chomp($line);
		my ($readName, $weight) = split(/\t/, $line);
		$readNameWeightHash{$readName} = $weight;
	}
	close($reader);
}

my ($totalReadCount, %orthologyCountHash, %orthologyRpkHash) = (0);
my %topTokenHash = ();
$topTokenHash{'qseqid'} = '';
my %orthologyProteinLengthHash = ();
open(my $reader, "cat @inputFileList | sort --field-separator='\t' -k1,1 -k11,11g -k12,12gr |");
while(my $line = <$reader>) {
	chomp($line);
	my %tokenHash = ();
	@tokenHash{qw(qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore)} = split(/\t/, $line);
	next if($tokenHash{'pident'} < $minimumPercentIdentity);
	if($tokenHash{'qseqid'} ne $topTokenHash{'qseqid'}) {
		if($topTokenHash{'qseqid'} ne '' && defined(my $weight = $readNameWeightFile ne '' ? $readNameWeightHash{$topTokenHash{'qseqid'}} : 1)) {
			if(scalar(my @orthologyList = keys %orthologyProteinLengthHash) == 1) {
				$orthologyCountHash{$_} += $weight foreach(@orthologyList);
				$orthologyRpkHash{$_} += $weight / (min(values %{$orthologyProteinLengthHash{$_}}) * 3 / 1000) foreach(@orthologyList);
			}
			$totalReadCount += $weight;
		}
		%topTokenHash = %tokenHash;
		%orthologyProteinLengthHash = ();
	}
	if(scalar(grep {$tokenHash{$_} != $topTokenHash{$_}} 'evalue', 'bitscore') == 0) {
		my $orthology = getOrthology(my $protein = $tokenHash{'sseqid'});
		$orthologyProteinLengthHash{$orthology}->{$protein} = $_ if(defined($_ = $proteinLengthHash{$protein}));
	}
}
close($reader);
if($topTokenHash{'qseqid'} ne '' && defined(my $weight = $readNameWeightFile ne '' ? $readNameWeightHash{$topTokenHash{'qseqid'}} : 1)) {
	if(scalar(my @orthologyList = keys %orthologyProteinLengthHash) == 1) {
		$orthologyCountHash{$_} += $weight foreach(@orthologyList);
		$orthologyRpkHash{$_} += $weight / (min(values %{$orthologyProteinLengthHash{$_}}) * 3 / 1000) foreach(@orthologyList);
	}
	$totalReadCount += $weight;
}

if($orthologyDefinitionFile ne '') {
	print join("\t", 'orthology', 'definition', 'count', 'rpkm'), "\n";
	foreach my $orthology (sort keys %orthologyCountHash) {
		my $definition = defined($_ = $orthologyDefinitionHash{$orthology}) ? $_ : '';
		my $count = defined($_ = $orthologyCountHash{$orthology}) ? $_ : 0;
		if($cpmInsteadOfRPKM) {
			my $cpm = $count / ($totalReadCount / 1000000);
			print join("\t", $orthology, $definition, $count, $cpm), "\n";
		} else {
			my $rpkm = (defined($_ = $orthologyRpkHash{$orthology}) ? $_ : 0) / ($totalReadCount / 1000000);
			print join("\t", $orthology, $definition, $count, $rpkm), "\n";
		}
	}
} else {
	print join("\t", 'orthology', 'count', 'rpkm'), "\n";
	foreach my $orthology (sort keys %orthologyCountHash) {
		my $count = defined($_ = $orthologyCountHash{$orthology}) ? $_ : 0;
		if($cpmInsteadOfRPKM) {
			my $cpm = $count / ($totalReadCount / 1000000);
			print join("\t", $orthology, $count, $cpm), "\n";
		} else {
			my $rpkm = (defined($_ = $orthologyRpkHash{$orthology}) ? $_ : 0) / ($totalReadCount / 1000000);
			print join("\t", $orthology, $count, $rpkm), "\n";
		}
	}
}

sub loadDefaultDatabase {
	open(my $reader, "$fmapPath/FMAP_data/database");
	chomp(my $line = <$reader>);
	close($reader);
	return $1 if($line =~ /^([A-Za-z0-9_.]+)/);
}

sub getOrthology {
	my ($protein) = @_;
	return $_ if(defined($_ = $proteinOrthologyHash{$protein}));
	return $1 if($protein =~ /^(K[0-9]{5})_/);
	return $protein;
}
