# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die "ERROR in $0: ", $_[0] };

use Cwd 'abs_path';
use Getopt::Long qw(:config no_ignore_case);

(my $fmapPath = abs_path($0)) =~ s/\/[^\/]*$//;
my $database = loadDefaultDatabase();
my $databasePrefix = "$fmapPath/FMAP_data/$database";

GetOptions('h' => \(my $help = ''),
	'm=s' => \(my $mapperPath = 'diamond'),
	'p=i' => \(my $threads = 1),
	'e=f' => \(my $evalue = 10),
	't=s' => \(my $temporaryDirectory = defined($ENV{'TMPDIR'}) ? $ENV{'TMPDIR'} : '/tmp'),
	'a=f' => \(my $acceleration = 0.5),
	'F=i' => \(my $frameshift = ''),
	'multiple' => \(my $multiple = ''),
	'd=s' => \$databasePrefix);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl FMAP_mapping.pl [options] input1.fastq|input1.fasta [input2.fastq|input2.fasta [...]] > blastx_hits.txt

Options: -h       display this help message
         -m FILE  executable file path of mapping program, "diamond" or "usearch" [$mapperPath]
         -p INT   number of threads [$threads]
         -e FLOAT maximum e-value to report alignments [$evalue]
         -t DIR   directory for temporary files [\$TMPDIR or /tmp]
         -a FLOAT search acceleration for ublast [$acceleration]
         -F INT   frameshift option for diamond

EOF
}
(my $mapper = $mapperPath) =~ s/^.*\///;
die "ERROR in $0: The mapper must be \"diamond\" or \"usearch\".\n" unless($mapper =~ /^diamond/ || $mapper =~ /^usearch/);
die "ERROR in $0: The mapper is not executable.\n" unless(-x getCommandPath($mapperPath));

die "ERROR in $0: The temporary directory is not a writable directory.\n" unless(-d $temporaryDirectory && -w $temporaryDirectory);

my $databaseFile = '';
$databaseFile = "$databasePrefix.dmnd" if($mapper =~ /^diamond/);
$databaseFile = "$databasePrefix.udb" if($mapper =~ /^usearch/);
die "ERROR in $0: The database is not available.\n" unless(-r $databaseFile);

my $temporaryOutputPrefix = "$temporaryDirectory/FMAP_mapping";
$temporaryOutputPrefix = "$temporaryOutputPrefix.$ENV{'HOSTNAME'}" if(defined($ENV{'HOSTNAME'}));
$temporaryOutputPrefix = "$temporaryOutputPrefix.$$";

my (@inputFileList) = @ARGV;
foreach my $inputFile (@inputFileList) {
	print STDERR "ERROR in $0: The input \"$inputFile\" is not available.\n" unless(-r $inputFile);
	if($mapper =~ /^diamond/) {
		if($frameshift) {
			system("$mapperPath blastx --query $inputFile --db $databaseFile --out $temporaryOutputPrefix.blast.txt --outfmt 6 --evalue $evalue --threads $threads --tmpdir $temporaryDirectory --quiet --frameshift $frameshift 1>&2");
		} elsif($multiple) {
			system("$mapperPath blastx --query $inputFile --db $databaseFile --out $temporaryOutputPrefix.blast.txt --outfmt 6 --evalue $evalue --threads $threads --tmpdir $temporaryDirectory --quiet --max-target-seqs 0 1>&2");
		} else {
			system("$mapperPath blastx --query $inputFile --db $databaseFile --out $temporaryOutputPrefix.blast.txt --outfmt 6 --evalue $evalue --threads $threads --tmpdir $temporaryDirectory --quiet --max-target-seqs 1 1>&2");
		}
	}
	if($mapper =~ /^usearch/) {
		if($multiple) {
			system("$mapperPath -ublast $inputFile -db $databaseFile -blast6out $temporaryOutputPrefix.blast.txt -evalue $evalue -accel $acceleration -threads $threads 1>&2");
		} else {
			system("$mapperPath -ublast $inputFile -db $databaseFile -blast6out $temporaryOutputPrefix.blast.txt -evalue $evalue -accel $acceleration -threads $threads -top_hits_only 1>&2");
		}
	}
	if(-r "$temporaryOutputPrefix.blast.txt") {
		open(my $reader, "$temporaryOutputPrefix.blast.txt");
		print while(<$reader>);
		close($reader);
		system("rm $temporaryOutputPrefix.blast.txt");
	}
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
