#!/usr/bin/env perl
use strict;
use warnings;
use Cwd 'abs_path';
use Getopt::Long;

(my $fmapPath = abs_path($0)) =~ s/\/[^\/]*$//;
my $database = 'orthology_uniref90_Bacteria_Archaea_Fungi';

GetOptions('h' => \(my $help = ''), 
	'm=s' => \(my $mapperPath = 'diamond'),
	'd=s' => \(my $databasePath = ''));
if($help) {
	die <<EOF;

Usage:   perl FMAP_download.pl [options]

Options: -h       display this help message
         -m FILE  executable file path of mapping program, "diamond" or "usearch"
         -d       full path for database location [$fmapPath]

EOF
}
(my $mapper = $mapperPath) =~ s/^.*\///;

$databasePath =~ s/\/$//;
unless ($databasePath && -d $databasePath) {
	$databasePath = $fmapPath;
	die "ERROR: Please check permissions, cannot write to database path: " . $databasePath . "\n" unless (-d $databasePath && -w $databasePath);
	system("mkdir -p $databasePath/FMAP_data");
	$databasePath .= "/FMAP_data";
}
die "ERROR: Please check permissions, cannot write to database path: " . $databasePath . "\n" unless (-d $databasePath && -w $databasePath);


die "ERROR: The mapper must be \"diamond\" or \"usearch\".\n" unless($mapper eq 'diamond' || $mapper eq 'usearch');
die "ERROR: The mapper is not executable.\n" unless(-x getCommandPath($mapperPath));

my $fmapURL = 'http://qbrc.swmed.edu/FMAP';

system("wget -q -O - $fmapURL/FMAP_data/version > $databasePath/version");
system("wget -q -O - $fmapURL/FMAP_data/$database.fasta.gz > $databasePath/$database.fasta.gz");
system("gzip -df $databasePath/$database.fasta.gz");
if(-s "$databasePath/$database.fasta") {
	my ($sequenceName, %sequenceLengthHash) = ('');
	open(my $reader, "$databasePath/$database.fasta");
	while(my $line = <$reader>) {
		chomp($line);
		if($line =~ /^>(.*)$/) {
			$sequenceName = $1;
		} else {
			$sequenceLengthHash{$sequenceName} += length($line);
		}
	}
	close($reader);
	open(my $writer, "> $databasePath/$database.length.txt");
	foreach my $sequenceName (sort keys %sequenceLengthHash) {
		print $writer join("\t", $sequenceName, $sequenceLengthHash{$sequenceName}), "\n";
	}
	close($writer);

	if($mapper eq 'diamond') {
		system("$mapperPath makedb --db $databasePath/$database.dmnd --in $databasePath/$database.fasta");
	}
	if($mapper eq 'usearch') {
		system("$mapperPath -makeudb_usearch $databasePath/$database.fasta -output $databasePath/$database.udb");
	}
	saveDefaultDatabase();
}
{
	open(my $reader, 'wget -q -O - http://rest.kegg.jp/link/pathway/ko |');
	open(my $writer, "> $databasePath/KEGG_orthology2pathway.txt");
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
{
	open(my $reader, 'wget -q -O - http://rest.kegg.jp/list/ko |');
	open(my $writer, "> $databasePath/KEGG_orthology.txt");
	while(my $line = <$reader>) {
		chomp($line);
		my ($orthology, $definition) = split(/\t/, $line);
		$orthology =~ s/^ko://;
		print $writer join("\t", $orthology, $definition), "\n";
	}
	close($reader);
	close($writer);
}
{
	open(my $reader, 'wget -q -O - http://rest.kegg.jp/list/pathway |');
	open(my $writer, "> $databasePath/KEGG_pathway.txt");
	while(my $line = <$reader>) {
		chomp($line);
		my ($pathway, $definition) = split(/\t/, $line);
		$pathway =~ s/^path://;
		print $writer join("\t", $pathway, $definition), "\n";
	}
	close($reader);
	close($writer);
}
#foreach my $file ('KEGG_orthology2pathway.txt', 'KEGG_orthology.txt', 'KEGG_pathway.txt') {
#	system("wget -q -O - $fmapURL/FMAP_data/$file > $fmapPath/FMAP_data/$file");
#}
foreach my $file ('known_operon.KEGG_orthology.txt', 'known_operon.definition.txt', 'known_operon.KEGG_orthology_definition.txt', 'known_operon.KEGG_pathway.txt') {
	system("wget -q -O - $fmapURL/FMAP_data/$file > $databasePath/$file");
}

sub getCommandPath {
	my ($command) = @_;
	chomp($command = `which $command`) if($command ne '' && $command !~ /\//);
	$command =~ s/^~\//$ENV{'HOME'}\//;
	$command = abs_path($command) if($command ne '');
	return $command;
}

sub saveDefaultDatabase {
	open(my $writer, "> $databasePath/database");
	print $writer $database;
	close($writer);
}
