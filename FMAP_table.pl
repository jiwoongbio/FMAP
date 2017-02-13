# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
use Getopt::Long;
use List::Util qw(sum);

GetOptions('h' => \(my $help = ''),
	'c' => \(my $countInsteadOfRPKM = ''),
	'f' => \(my $fraction = ''),
	'n' => \(my $noDefinition = ''));
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl FMAP_table.pl [options] [name1=]abundance1.txt [[name2=]abundance2.txt [...]] > abundance_table.txt

Options: -h       display this help message
         -c       use raw counts instead of RPKM values
         -f       use fractions

EOF
}

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
	die "ERROR: The input \"$sampleFile\" is not available.\n" unless(-r $sampleFile);
	open(my $reader, $sampleFile);
	chomp(my $line = <$reader>);
	my @columnList = split(/\t/, $line);
	while(my $line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		@tokenHash{@columnList} = split(/\t/, $line);
		$orthologyHash{my $orthology = $tokenHash{'orthology'}} = 1;
		$orthologyDefinitionHash{$orthology} = $_ if(defined($_ = $tokenHash{'definition'}) && $noDefinition eq '');
		my $abundance = 0;
		if($countInsteadOfRPKM) {
			$abundance = $tokenHash{'count'};
		} else {
			$abundance = $tokenHash{'rpkm'};
		}
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
