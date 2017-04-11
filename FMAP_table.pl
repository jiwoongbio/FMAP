# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die "ERROR in $0: ", $_[0] };

use Getopt::Long;
use List::Util qw(sum);

my $optionTargetColumn = 'RPKM|rpkm';
GetOptions('h' => \(my $help = ''),
	'column=s' => $optionTargetColumn,
	'c' => \(my $countInsteadOfRPKM = ''),
	'd' => \(my $depthInsteadOfRPKM = ''),
	'f' => \(my $fraction = ''),
	'n' => \(my $noDefinition = ''),
	'r' => \(my $printRegion = ''));
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl FMAP_table.pl [options] [name1=]abundance1.txt [[name2=]abundance2.txt [...]] > abundance_table.txt

Options: -h       display this help message
         -c       use raw read counts (readCount|count) instead of RPKM values
         -d       use normalized mean depths (meanDepth/genome) instead of RPKM values
         -f       use fractions
         -n       do not print definitions
         -r       print ORF regions

EOF
}
$optionTargetColumn = 'readCount|count' if($countInsteadOfRPKM);
$optionTargetColumn = 'meanDepth/genome' if($depthInsteadOfRPKM);

my @sampleFileList = @ARGV;
my @sampleNameList = @sampleFileList;
s/=.*$// foreach(@sampleNameList);
s/^.*=// foreach(@sampleFileList);

my %orthologyHash = ();
my %orthologyDefinitionHash = ();
my %orthologyAbundanceListHash = ();
my @abundanceSumList = ();
foreach my $index (0 .. $#sampleFileList) {
	my $sampleFile = $sampleFileList[$index];
	die "ERROR in $0: '$sampleFile' is not readable.\n" unless(-r $sampleFile);
	open(my $reader, $sampleFile);
	chomp(my $line = <$reader>);
	my @columnList = split(/\t/, $line);
	my %columnHash = map {$_ => 1} @columnList;
	my ($targetColumn) = grep {$columnHash{$_}} split(/\|/, $optionTargetColumn);
	die "ERROR in $0: '$sampleFile' does not include the column: $optionTargetColumn.\n" unless(defined($targetColumn));
	die "ERROR in $0: '$sampleFile' does not include the column: $_.\n" if($printRegion && ($_ = join(', ', grep {!defined($columnHash{$_})} 'contig', 'start', 'end', 'strand')) ne '');
	while(my $line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		@tokenHash{@columnList} = split(/\t/, $line);
		my $orthology = $tokenHash{'orthology'};
		$orthology = join("\t", @tokenHash{'contig', 'start', 'end', 'strand', 'orthology'}) if($printRegion);
		$orthologyHash{$orthology} = 1;
		$orthologyDefinitionHash{$orthology} = $_ if(defined($_ = $tokenHash{'definition'}) && $noDefinition eq '');
		my $abundance = $tokenHash{$targetColumn};
		$orthologyAbundanceListHash{$orthology}->[$index] += $abundance;
		$abundanceSumList[$index] += $abundance;
	}
	close($reader);
}

if(scalar(keys %orthologyDefinitionHash) > 0) {
	if($printRegion) {
		print join("\t", 'contig', 'start', 'end', 'strand', 'orthology', 'definition', @sampleNameList), "\n";
	} else {
		print join("\t", 'orthology', 'definition', @sampleNameList), "\n";
	}
	foreach my $orthology (sort keys %orthologyHash) {
		my $definition = defined($_ = $orthologyDefinitionHash{$orthology}) ? $_ : '';
		my @abundanceList = map {defined($_) ? $_ : 0} @{$orthologyAbundanceListHash{$orthology}}[0 .. $#sampleNameList];
		@abundanceList = map {$abundanceList[$_] / $abundanceSumList[$_]} 0 .. $#abundanceList if($fraction);
		print join("\t", $orthology, $definition, @abundanceList), "\n";
	}
} else {
	if($printRegion) {
		print join("\t", 'contig', 'start', 'end', 'strand', 'orthology', @sampleNameList), "\n";
	} else {
		print join("\t", 'orthology', @sampleNameList), "\n";
	}
	foreach my $orthology (sort keys %orthologyHash) {
		my @abundanceList = map {defined($_) ? $_ : 0} @{$orthologyAbundanceListHash{$orthology}}[0 .. $#sampleNameList];
		@abundanceList = map {$abundanceList[$_] / $abundanceSumList[$_]} 0 .. $#abundanceList if($fraction);
		print join("\t", $orthology, @abundanceList), "\n";
	}
}
