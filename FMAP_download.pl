# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
use Cwd 'abs_path';
use Getopt::Long;

(my $fmapPath = abs_path($0)) =~ s/\/[^\/]*$//;
my $database = 'orthology_uniref90_Bacteria_Archaea_Fungi';
my $databasePath = "$fmapPath/FMAP_data/$database";

GetOptions('h' => \(my $help = ''), 'm=s' => \(my $mapperPath = 'diamond'));
if($help) {
	die <<EOF;

Usage:   perl FMAP_download.pl [options]

Options: -h       display this help message
         -m FILE  executable file path of mapping program, "diamond" or "usearch"

EOF
}
(my $mapper = $mapperPath) =~ s/^.*\///;
die "ERROR: The mapper must be \"diamond\" or \"usearch\".\n" unless($mapper eq 'diamond' || $mapper eq 'usearch');
die "ERROR: The mapper is not executable.\n" unless(-x getCommandPath($mapperPath));

my $fmapURL = 'http://qbrc.swmed.edu/FMAP';

system("mkdir -p $fmapPath/FMAP_data");
system("wget -q -O - $fmapURL/FMAP_data/version > $fmapPath/FMAP_data/version");
system("wget -q -O - $fmapURL/FMAP_data/$database.fasta.gz > $databasePath.fasta.gz");
system("gzip -df $databasePath.fasta.gz");
if(-s "$databasePath.fasta") {
	my ($sequenceName, %sequenceLengthHash) = ('');
	open(my $reader, "$databasePath.fasta");
	while(my $line = <$reader>) {
		chomp($line);
		if($line =~ /^>(.*)$/) {
			$sequenceName = $1;
		} else {
			$sequenceLengthHash{$sequenceName} += length($line);
		}
	}
	close($reader);
	open(my $writer, "> $databasePath.length.txt");
	foreach my $sequenceName (sort keys %sequenceLengthHash) {
		print $writer join("\t", $sequenceName, $sequenceLengthHash{$sequenceName}), "\n";
	}
	close($writer);

	if($mapper eq 'diamond') {
		system("$mapperPath makedb --db $databasePath.dmnd --in $databasePath.fasta");
	}
	if($mapper eq 'usearch') {
		system("$mapperPath -makeudb_usearch $databasePath.fasta -output $databasePath.udb");
	}
	saveDefaultDatabase();
}
{
	open(my $reader, 'wget -q -O - http://rest.kegg.jp/link/pathway/ko |');
	open(my $writer, "> $fmapPath/FMAP_data/KEGG_orthology2pathway.txt");
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
	open(my $writer, "> $fmapPath/FMAP_data/KEGG_orthology.txt");
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
	open(my $writer, "> $fmapPath/FMAP_data/KEGG_pathway.txt");
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
	system("wget -q -O - $fmapURL/FMAP_data/$file > $fmapPath/FMAP_data/$file");
}

sub getCommandPath {
	my ($command) = @_;
	chomp($command = `which $command`) if($command ne '' && $command !~ /\//);
	$command =~ s/^~\//$ENV{'HOME'}\//;
	$command = abs_path($command) if($command ne '');
	return $command;
}

sub saveDefaultDatabase {
	open(my $writer, "> $fmapPath/FMAP_data/database");
	print $writer $database;
	close($writer);
}
