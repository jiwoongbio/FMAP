# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Getopt::Long;
use List::Util qw(sum);

my $column = 'RPKM|rpkm';
GetOptions('h' => \(my $help = ''),
	'column=s' => $column,
	'c' => \(my $countInsteadOfRPKM = ''),
	'd' => \(my $depthInsteadOfRPKM = ''),
	'f' => \(my $fraction = ''),
	'n' => \(my $noDefinition = ''));
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl FMAP_table.pl [options] [name1=]abundance1.txt [[name2=]abundance2.txt [...]] > abundance_table.txt

Options: -h       display this help message
         -c       use raw read counts (readCount|count) instead of RPKM values
         -d       use normalized mean depths (meanDepth/genome) instead of RPKM values
         -f       use fractions

EOF
}
$column = 'readCount|count' if($countInsteadOfRPKM);
$column = 'meanDepth/genome' if($depthInsteadOfRPKM);

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
	die "ERROR: '$sampleFile' is not readable.\n" unless(-r $sampleFile);
	open(my $reader, $sampleFile);
	chomp(my $line = <$reader>);
	my @columnList = split(/\t/, $line);
	my %columnHash = map {$_ => 1} @columnList;
	($column) = grep {$columnHash{$_}} split(/\|/, $column);
	die "ERROR: '$column' is not a column.\n" unless(defined($column));
	while(my $line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		@tokenHash{@columnList} = split(/\t/, $line);
		$orthologyHash{my $orthology = $tokenHash{'orthology'}} = 1;
		$orthologyDefinitionHash{$orthology} = $_ if(defined($_ = $tokenHash{'definition'}) && $noDefinition eq '');
		my $abundance = $tokenHash{$column};
		$orthologyAbundanceListHash{$orthology}->[$index] = $abundance;
		$abundanceSumList[$index] += $abundance;
	}
	close($reader);
}

if(scalar(keys %orthologyDefinitionHash) > 0) {
	print join("\t", 'orthology', 'definition', @sampleNameList), "\n";
	foreach my $orthology (sort keys %orthologyHash) {
		my $definition = defined($_ = $orthologyDefinitionHash{$orthology}) ? $_ : '';
		my @abundanceList = map {defined($_) ? $_ : 0} @{$orthologyAbundanceListHash{$orthology}};
		@abundanceList = map {$abundanceList[$_] / $abundanceSumList[$_]} 0 .. $#abundanceList if($fraction);
		print join("\t", $orthology, $definition, @abundanceList), "\n";
	}
} else {
	print join("\t", 'orthology', @sampleNameList), "\n";
	foreach my $orthology (sort keys %orthologyHash) {
		my @abundanceList = map {defined($_) ? $_ : 0} @{$orthologyAbundanceListHash{$orthology}};
		@abundanceList = map {$abundanceList[$_] / $abundanceSumList[$_]} 0 .. $#abundanceList if($fraction);
		print join("\t", $orthology, @abundanceList), "\n";
	}
}
