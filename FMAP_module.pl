# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die "ERROR in $0: ", $_[0] };

use Cwd 'abs_path';
use Getopt::Long;
use Statistics::R;

(my $fmapPath = abs_path($0)) =~ s/\/[^\/]*$//;

GetOptions('h' => \(my $help = ''),
	'p=s' => \(my $orthology2moduleFile = "$fmapPath/FMAP_data/KEGG_orthology2module.txt"),
	'd=s' => \(my $moduleDefinitionFile = "$fmapPath/FMAP_data/KEGG_module.txt"),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl FMAP_module.pl [options] orthology_test_stat.txt > module.txt

Options: -h       display this help message

EOF
}

my ($inputFile) = @ARGV;
foreach($inputFile, $orthology2moduleFile, $moduleDefinitionFile) {
	die "ERROR in $0: '$_' is not readable.\n" unless(-r $_);
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
		if(not defined($tokenHash{'filter'}) or $tokenHash{'filter'} eq 'pass') {
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
open(my $reader, $orthology2moduleFile);
while(my $line = <$reader>) {
	chomp($line);
	my ($orthology, $module) = split(/\t/, $line);
	my $color = $targetOrthologyColorHash{$orthology};
	push(@{$moduleTargetOrthologyListHash{$module}}, "$orthology+$color") if(defined($color));
	$moduleOrthologyCountHash{$module} += 1;
	$orthologyHash{$orthology} = 1;
}
close($reader);
my $orthologyCount = scalar(keys %orthologyHash);

my %moduleDefinitionHash = ();
{
	open(my $reader, $moduleDefinitionFile);
	while(my $line = <$reader>) {
		chomp($line);
		my ($module, $definition) = split(/\t/, $line);
		$moduleDefinitionHash{$module} = $definition;
	}
	close($reader);
}

print join("\t", 'module', 'definition', 'orthology.count', 'coverage', 'pvalue', 'weblink'), "\n";
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
	my $weblink = sprintf("http://www.kegg.jp/kegg-bin/show_pathway?map=%s&multi_query=%s", $module, join('%0d%0a', @moduleTargetOrthologyList));
	my %countHash = ();
	foreach my $orthology (@moduleTargetOrthologyList) {
		$countHash{'up'} += 1 if($orthology =~ /\+red$/);
		$countHash{'down'} += 1 if($orthology =~ /\+blue$/);
	}
	$moduleTargetOrthologyCount = sprintf('%d (%s)', $moduleTargetOrthologyCount, $_) if($_ = join(', ', map {defined($_->[1]) ? join(':', @$_) : ()} map {[$_, $countHash{$_}]} ('up', 'down')));
	print join("\t", $module, $definition, $moduleTargetOrthologyCount, $coverage, $pvalue, $weblink), "\n";
}
$R->stop();
