#!/usr/bin/env perl

use strict;
use warnings;
use Cwd 'abs_path';
use Getopt::Long;

(my $fmapPath = abs_path($0)) =~ s/\/[^\/]*$//;

my @alignmentFileList = ();
GetOptions('h' => \(my $help = ''),
	'm=s' => \(my $mapperPath = 'diamond'),
	'p=i' => \(my $mappingThreads = 1),
	'e=f' => \(my $evalue = 0.001),
	'i=f' => \(my $identity = 0.8),
	't=s' => \(my $temporaryDirectory = defined($ENV{'TMPDIR'}) ? $ENV{'TMPDIR'} : '/tmp'),
	'a=s' => \@alignmentFileList,
	'multiple' => \(my $multiple = ''),
	'd=s' => \(my $databasePath = ''));
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl FMAP_mapping.pl [options] input1.fastq|input1.fasta [input2.fastq|input2.fasta [...]] > blastx_hits.txt

Options: -h       display this help message
         -m FILE  executable file path of mapping program, "diamond" or "usearch" [diamond]
         -p INT   number of threads [1]
         -e FLOAT maximum e-value to report alignments for "diamond" [0.001]
         -i FLOAT minimum identity for "usearch_global" [0.8]
         -t DIR   directory for temporary files [\$TMPDIR or /tmp]
         -d       full path for database location [$fmapPath/FMAP_data]

EOF
}

unless ($databasePath ne '') {
	$databasePath = $fmapPath . "/FMAP_data";
};
die "Database path does not exist: " . $databasePath . "\n" unless (-d $databasePath);

my $database = loadDefaultDatabase($databasePath);
$databasePath =~ s/\/$//;
$databasePath .= "/" . $database;
die "ERROR: The database is not available.\n" unless (-r $databasePath . ".dmnd" || -r $databasePath . ".udb");

(my $mapper = $mapperPath) =~ s/^.*\///;
die "ERROR: The mapper must be \"diamond\" or \"usearch\".\n" unless($mapper eq 'diamond' || $mapper eq 'usearch');
die "ERROR: The mapper is not executable.\n" unless(-x getCommandPath($mapperPath));

die "ERROR: The temporary directory is not available.\n" unless(-d $temporaryDirectory && -w $temporaryDirectory);

my $temporaryOutputPrefix = "$temporaryDirectory/FMAP_mapping";
$temporaryOutputPrefix = "$temporaryOutputPrefix.$ENV{'HOSTNAME'}" if(defined($ENV{'HOSTNAME'}));
$temporaryOutputPrefix = "$temporaryOutputPrefix.$$";

my (@inputFileList) = @ARGV;
if($mapper eq 'diamond') {
	die "ERROR: The database is not available.\n" unless(-r "$databasePath.dmnd");
	foreach my $index (0 .. $#inputFileList) {
		my $inputFile = $inputFileList[$index];
		print STDERR "ERROR: The input \"$inputFile\" is not available.\n" unless(-r $inputFile);
		if($multiple) {
			system("$mapperPath blastx --db $databasePath.dmnd --threads $mappingThreads --query $inputFile --daa $temporaryOutputPrefix.daa --tmpdir $temporaryDirectory --evalue $evalue 1>&2");
		} else {
			system("$mapperPath blastx --db $databasePath.dmnd --threads $mappingThreads --query $inputFile --daa $temporaryOutputPrefix.daa --tmpdir $temporaryDirectory --max-target-seqs 1 --evalue $evalue 1>&2");
		}
		if(-r "$temporaryOutputPrefix.daa") {
			system("$mapperPath view --daa $temporaryOutputPrefix.daa");
			if(my $alignmentFile = $alignmentFileList[$index]) {
				system("mv $temporaryOutputPrefix.daa $alignmentFile");
			} else {
				system("rm $temporaryOutputPrefix.daa");
			}
		}
	}
}
if($mapper eq 'usearch') {
	die "ERROR: The database is not available.\n" unless(-r "$databasePath.udb");
	foreach my $inputFile (@inputFileList) {
		print STDERR "ERROR: The input \"$inputFile\" is not available.\n" unless(-r $inputFile);
		system("$mapperPath -usearch_global $inputFile -db $databasePath.udb -id $identity -blast6out $temporaryOutputPrefix.b6 -top_hits_only -threads $mappingThreads 1>&2");
		if(-r "$temporaryOutputPrefix.b6") {
			open(my $reader, "$temporaryOutputPrefix.b6");
			print while(<$reader>);
			close($reader);
			system("rm $temporaryOutputPrefix.b6");
		}
	}
}

sub loadDefaultDatabase {
	die "Could not find database path: " . $_[0] . "\n" unless ($_[0]);

	$_[0] =~ s/\/$//;
	my $dbFile = $_[0] . "/database";
	die "Could not find database file: " . $dbFile . "\n" unless (-r $dbFile);

	open(my $reader, $dbFile);
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
