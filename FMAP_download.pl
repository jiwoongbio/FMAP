# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die "ERROR in $0: ", $_[0] };

use Cwd 'abs_path';
use Getopt::Long;

(my $fmapPath = abs_path($0)) =~ s/\/[^\/]*$//;
my $fmapURL = 'http://qbrc.swmed.edu/FMAP';
system("mkdir -p $fmapPath/FMAP_data");

GetOptions('h' => \(my $help = ''),
	'm=s' => \(my $mapperPath = 'diamond'),
	'k' => \(my $downloadPrebuiltKEGG = ''),
	'x' => \(my $downloadOnlyKEGG = ''));
if($help) {
	die <<EOF;

Usage:   perl FMAP_download.pl [options]

Options: -h       display this help message
         -m FILE  executable file path of mapping program, "diamond" or "usearch" [$mapperPath]
         -k       download prebuilt KEGG files
         -x       download only KEGG files

EOF
}
(my $mapper = $mapperPath) =~ s/^.*\///;
die "ERROR in $0: 'wget' is not executable.\n" unless(-x getCommandPath('wget'));

# Database
unless($downloadOnlyKEGG) {
	die "ERROR in $0: The mapper must be \"diamond\" or \"usearch\".\n" unless($mapper =~ /^diamond/ || $mapper =~ /^usearch/);
	die "ERROR in $0: The mapper is not executable.\n" unless(-x getCommandPath($mapperPath));

	downloadFile('database');
	my $database = loadDefaultDatabase();
	my $databasePath = "$fmapPath/FMAP_data/$database";
	downloadFile("$database.fasta.gz");
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
	
		if($mapper =~ /^diamond/) {
			system("$mapperPath makedb --in $databasePath.fasta --db $databasePath.dmnd");
		}
		if($mapper =~ /^usearch/) {
			system("$mapperPath -makeudb_ublast $databasePath.fasta -output $databasePath.udb");
		}
	}
}

# KEGG
if($downloadPrebuiltKEGG) {
	downloadFile('KEGG_orthology2pathway.txt');
} else {
	open(my $reader, 'wget --no-verbose -O - http://rest.kegg.jp/link/pathway/ko |');
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
if($downloadPrebuiltKEGG) {
	downloadFile('KEGG_orthology.txt');
} else {
	open(my $reader, 'wget --no-verbose -O - http://rest.kegg.jp/list/ko |');
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
if($downloadPrebuiltKEGG) {
	downloadFile('KEGG_pathway.txt');
} else {
	open(my $reader, 'wget --no-verbose -O - http://rest.kegg.jp/list/pathway |');
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

# Operon
unless($downloadOnlyKEGG) {
	downloadFile('known_operon.KEGG_orthology.txt', 'known_operon.definition.txt', 'known_operon.KEGG_orthology_definition.txt', 'known_operon.KEGG_pathway.txt');
}

sub downloadFile {
	system("wget --no-verbose -O $fmapPath/FMAP_data/$_ $fmapURL/FMAP_data/$_") foreach(@_);
}

sub loadDefaultDatabase {
	open(my $reader, "$fmapPath/FMAP_data/database");
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
