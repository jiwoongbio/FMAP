# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
use Cwd 'abs_path';
use Getopt::Long;
use Statistics::R;

(my $fmapPath = abs_path($0)) =~ s/\/[^\/]*$//;

GetOptions('h' => \(my $help = ''));
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl FMAP_module.pl [options] orthology_test_stat.txt > module.txt

Options: -h       display this help message

EOF
}

my ($inputFile) = @ARGV;
foreach($inputFile, "$fmapPath/FMAP_data/KEGG_orthology2module.txt", "$fmapPath/FMAP_data/KEGG_module.txt") {
	die "ERROR: '$_' is not readable.\n" unless(-r $_);
}
my %targetOrthologyColorHash = ();
{
	open(my $reader, $inputFile);
	chomp(my $line = <$reader>);
	my @columnList = split(/\t/, $line);
	while(my $line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		@tokenHash{@columnList} = split(/\t/, $line);
		if($tokenHash{'filter'} eq 'pass') {
			my ($orthology, $log2foldchange) = @tokenHash{'orthology', 'log2foldchange'};
			$targetOrthologyColorHash{$orthology} = 'red'  if($log2foldchange > 0);
			$targetOrthologyColorHash{$orthology} = 'blue' if($log2foldchange < 0);
		}
	}
	close($reader);
}
my $targetOrthologyCount = scalar(keys %targetOrthologyColorHash);

my %moduleTargetOrthologyListHash = ();
my %moduleOrthologyCountHash = ();
my %orthologyHash = ();
open(my $reader, "$fmapPath/FMAP_data/KEGG_orthology2module.txt");
while(my $line = <$reader>) {
	chomp($line);
	my ($orthology, $module) = split(/\t/, $line);
	my $color = $targetOrthologyColorHash{$orthology};
	push(@{$moduleTargetOrthologyListHash{$module}}, "$orthology $color") if(defined($color));
	$moduleOrthologyCountHash{$module} += 1;
	$orthologyHash{$orthology} = 1;
}
close($reader);
my $orthologyCount = scalar(keys %orthologyHash);

my %moduleDefinitionHash = ();
{
	open(my $reader, "$fmapPath/FMAP_data/KEGG_module.txt");
	while(my $line = <$reader>) {
		chomp($line);
		my ($module, $definition) = split(/\t/, $line);
		$moduleDefinitionHash{$module} = $definition;
	}
	close($reader);
}

print join("\t", 'module', 'definition', 'orthology.count', 'coverage', 'pvalue', 'orthology.colors'), "\n";
my $R = Statistics::R->new();
foreach my $module (sort keys %moduleTargetOrthologyListHash) {
	my $moduleTargetOrthologyCount = scalar(my @moduleTargetOrthologyList = @{$moduleTargetOrthologyListHash{$module}});
	my $moduleOrthologyCount = $moduleOrthologyCountHash{$module};
	my $counts = join(',', $moduleTargetOrthologyCount, $moduleOrthologyCount - $moduleTargetOrthologyCount, $targetOrthologyCount - $moduleTargetOrthologyCount, $orthologyCount - $moduleOrthologyCount - $targetOrthologyCount + $moduleTargetOrthologyCount);
	$R->run("p.value <- fisher.test(matrix(c($counts), 2), alternative = \"greater\")\$p.value");
	my $pvalue = $R->get("p.value");
	my $definition = $moduleDefinitionHash{$module};
	$definition = '' unless(defined($definition));
	my $coverage = $moduleTargetOrthologyCount / $moduleOrthologyCount;
	print join("\t", $module, $definition, $moduleTargetOrthologyCount, $coverage, $pvalue, join('|', @moduleTargetOrthologyList)), "\n";
}
$R->stop();
