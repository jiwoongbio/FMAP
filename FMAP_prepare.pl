# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
use Cwd 'abs_path';
use Getopt::Long;

(my $fmapPath = abs_path($0)) =~ s/\/[^\/]*$//;
my $dataPath = "$fmapPath/FMAP_data";
system("mkdir -p $dataPath");

my $fmapURL = 'http://qbrc.swmed.edu/FMAP';

GetOptions('h' => \(my $help = ''),
	'r' => \(my $redownload = ''),
	'm=s' => \(my $mapperPath = 'diamond'),
	'k' => \(my $downloadPrebuiltKEGG = ''));
if($help) {
	die <<EOF;

Usage:   perl FMAP_prepare.pl [options]

Options: -h       display this help message
         -r       redownload data
         -m FILE  executable file path of mapping program, "diamond" or "usearch" [$mapperPath]
         -k       download prebuilt KEGG files

EOF
}
(my $mapper = $mapperPath) =~ s/^.*\///;
die "ERROR: 'wget' is not executable.\n" unless(-x getCommandPath('wget'));
die "ERROR: The mapper must be \"diamond\" or \"usearch\".\n" unless($mapper eq 'diamond' || $mapper eq 'usearch');
die "ERROR: The mapper is not executable.\n" unless(-x getCommandPath($mapperPath));

# Database
my ($database, $databasePath) = ('', '');
if(-r "$dataPath/database") {
	$database = loadDefaultDatabase();
	$databasePath = "$dataPath/$database";
} else {
	downloadFile('database');
	$database = loadDefaultDatabase();
	downloadFile("$database.fasta.gz");
	system("gzip -df $dataPath/$database.fasta.gz");
	$databasePath = "$dataPath/$database";
}
if(-r "$databasePath.fasta") {
	generateSequenceLengthFile("$databasePath.length.txt", "$databasePath.fasta");
	if($mapper eq 'diamond') {
		system("$mapperPath makedb --db $databasePath.dmnd --in $databasePath.fasta 1>&2");
	}
	if($mapper eq 'usearch') {
		system("$mapperPath -makeudb_usearch $databasePath.fasta -output $databasePath.udb 1>&2");
	}
} else {
	die "ERROR: The database is not available.\n";
}

# KEGG
{
	my $file = 'KEGG_orthology2pathway.txt';
	if($downloadPrebuiltKEGG) {
		downloadFile($file);
	} elsif(not -r "$dataPath/$file" or $redownload) {
		open(my $reader, 'wget --no-verbose -O - http://rest.kegg.jp/link/pathway/ko |');
		open(my $writer, "> $dataPath/$file");
		while(my $line = <$reader>) {
			chomp($line);
			my ($orthology, $pathway) = split(/\t/, $line);
			$orthology =~ s/^ko://;
			$pathway =~ s/^path://;
			if($pathway =~ /^map[0-9]+$/) {
				print $writer join("\t", $orthology, $pathway), "\n";
			}
		}
		close($reader);
		close($writer);
	}
}
{
	my $file = 'KEGG_orthology.txt';
	if($downloadPrebuiltKEGG) {
		downloadFile($file);
	} elsif(not -r "$dataPath/$file" or $redownload) {
		open(my $reader, 'wget --no-verbose -O - http://rest.kegg.jp/list/ko |');
		open(my $writer, "> $dataPath/$file");
		while(my $line = <$reader>) {
			chomp($line);
			my ($orthology, $definition) = split(/\t/, $line);
			$orthology =~ s/^ko://;
			print $writer join("\t", $orthology, $definition), "\n";
		}
		close($reader);
		close($writer);
	}
}
{
	my $file = 'KEGG_pathway.txt';
	if($downloadPrebuiltKEGG) {
		downloadFile($file);
	} elsif(not -r "$dataPath/$file" or $redownload) {
		open(my $reader, 'wget --no-verbose -O - http://rest.kegg.jp/list/pathway |');
		open(my $writer, "> $dataPath/$file");
		while(my $line = <$reader>) {
			chomp($line);
			my ($pathway, $definition) = split(/\t/, $line);
			$pathway =~ s/^path://;
			print $writer join("\t", $pathway, $definition), "\n";
		}
		close($reader);
		close($writer);
	}
}

# Operon
chomp(my $hostname = `hostname`);
if($fmapURL !~ /:\/\/$hostname\//) {
	downloadFile('known_operon.KEGG_orthology.txt', 'known_operon.definition.txt', 'known_operon.KEGG_orthology_definition.txt', 'known_operon.KEGG_pathway.txt');
}

sub downloadFile {
	foreach my $file (@_) {
		system("wget --no-verbose -O $dataPath/$file $fmapURL/FMAP_data/$file") if(not -r "$dataPath/$file" or $redownload);
	}
}

sub loadDefaultDatabase {
	open(my $reader, "$dataPath/database");
	chomp(my $line = <$reader>);
	close($reader);
	return $1 if($line =~ /^([A-Za-z0-9_.]+)/);
}

sub getCommandPath {
	my ($command) = @_;
	chomp($command = `which $command`) if($command ne '' && $command !~ /\//);
	$command =~ s/^~\//$ENV{'HOME'}\//;
	$command = abs_path($command) if($command ne '');
	return $command;
}

sub generateSequenceLengthFile {
	my ($sequenceLengthFile, $fastaFile) = @_;
	open(my $writer, "> $sequenceLengthFile");
	my ($sequenceName, $sequenceLength) = ('', 0);
	open(my $reader, $fastaFile);
	while(my $line = <$reader>) {
		chomp($line);
		if($line =~ /^>(.*)$/) {
			(my $nextSequenceName = $1) =~ s/ .*$//;
			print $writer join("\t", $sequenceName, $sequenceLength), "\n" if($sequenceLength > 0);
			($sequenceName, $sequenceLength) = ($nextSequenceName, 0);
		} else {
			$sequenceLength += length($line);
		}
	}
	close($reader);
	print $writer join("\t", $sequenceName, $sequenceLength), "\n" if($sequenceLength > 0);
	close($writer);
}
