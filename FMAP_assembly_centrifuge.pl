# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die "ERROR in $0: ", $_[0] };

use Getopt::Long qw(:config no_ignore_case);

GetOptions('h' => \(my $help = ''),
	'p=i' => \(my $threads = 1),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl FMAP_assembly_taxon.pl [options] FMAP_assembly.region.txt assembly.fasta centrifuge.index

Options: -h       display this help message
         -p INT   number of threads [$threads]

EOF
}
my ($regionFile, $assemblyFastaFile, $centrifugeIndex) = @ARGV;
my %assemblySequenceHash = ();
{
	open(my $reader, $assemblyFastaFile);
	my $contig = '';
	while(my $line = <$reader>) {
		chomp($line);
		if($line =~ /^>(\S*)/) {
			$contig = $1;
		} else {
			$assemblySequenceHash{$contig} .= $line;
		}
	}
	close($reader);
}
(my $regionFastaFile = $regionFile) =~ s/(\.txt)?$/\.fasta/;
{
	open(my $reader, $regionFile);
	open(my $writer, "> $regionFastaFile");
	chomp(my $line = <$reader>);
	my @columnList = split(/\t/, $line);
	while(my $line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		@tokenHash{@columnList} = split(/\t/, $line);
		my $region = join('|', my ($contig, $start, $end, $strand) = @tokenHash{'contig', 'start', 'end', 'strand'});
		my $sequence = substr($assemblySequenceHash{$contig}, $start - 1, $end - ($start - 1));
		$sequence = uc($sequence);
		($sequence = reverse($sequence)) =~ tr/ACGT/TGCA/ if($strand eq '-');
		print $writer ">$region\n";
		print $writer "$sequence\n";
	}
	close($reader);
	close($writer);
}
my %regionTaxonReferenceHash = ();
{
	open(my $reader, "centrifuge -f -p $threads -x $centrifugeIndex -U $regionFastaFile --report-file /dev/null |");
	chomp(my $line = <$reader>);
	my @columnList = split(/\t/, $line);
	while(my $line = <$reader>) {
		my %tokenHash = ();
		@tokenHash{@columnList} = split(/\t/, $line);
		$regionTaxonReferenceHash{$tokenHash{'readID'}} = [@tokenHash{'taxID', 'seqID'}];
	}
	close($reader);
}
{
	open(my $reader, $regionFile);
	chomp(my $line = <$reader>);
	my @columnList = split(/\t/, $line);
	print join("\t", @columnList, 'taxon', 'reference'), "\n";
	while(my $line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		@tokenHash{@columnList} = split(/\t/, $line);
		my $region = join('|', @tokenHash{'contig', 'start', 'end', 'strand'});
		@tokenHash{'taxon', 'reference'} = defined($_ = $regionTaxonReferenceHash{$region}) ? @$_ : ('', '');
		print join("\t", @tokenHash{@columnList, 'taxon', 'reference'}), "\n";
	}
	close($reader);
}
