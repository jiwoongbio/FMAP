#!/usr/bin/env perl
use strict;
use warnings;
use Cwd 'abs_path';
use Getopt::Long;
use List::Util qw(sum);
use Statistics::R;

(my $fmapPath = abs_path($0)) =~ s/\/[^\/]*$//;

GetOptions('h' => \(my $help = ''), 'a' => \(my $all = ''));
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl FMAP_operon.pl [options] orthology_test_stat.txt > operon.txt

Options: -h       display this help message
         -a       print single-gene operons as well

EOF
}

my ($inputFile) = @ARGV;
die "ERROR: The input \"$inputFile\" is not available.\n" unless(-r $inputFile);
my %orthologyFoldchangeHash = ();
my %orthologyFilterHash = ();
{
	open(my $reader, $inputFile);
	chomp(my $line = <$reader>);
	my @columnList = split(/\t/, $line);
	while(my $line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		@tokenHash{@columnList} = split(/\t/, $line);
		my $orthology = $tokenHash{'orthology'};
		$orthologyFoldchangeHash{$orthology} = $tokenHash{'log2foldchange'};
		$orthologyFilterHash{$orthology} = $tokenHash{'filter'};
	}
	close($reader);
}
my $passCount = scalar(grep {$_ eq 'pass'} values %orthologyFilterHash);
my $failCount = scalar(grep {$_ eq 'fail'} values %orthologyFilterHash);
my @operonsOrthologiesFoldchangeList = ();
my $R = Statistics::R->new();
open(my $reader, "$fmapPath/FMAP_data/known_operon.KEGG_orthology.txt");
while(my $line = <$reader>) {
	chomp($line);
	my ($operons, $orthologies) = split(/\t/, $line);
	my $orthologyCount = scalar(my @orthologyList = split(/ /, $orthologies));
	if($orthologyCount > 1 || $all) {
		my @foldchangeList = @orthologyFoldchangeHash{@orthologyList};
		if(scalar(grep {defined} @foldchangeList) == $orthologyCount) {
			my $log2foldchange = mean(@foldchangeList);
			my $orthologyPassCount = scalar(grep {$_ eq 'pass'} @orthologyFilterHash{@orthologyList});
			my $orthologyFailCount = scalar(grep {$_ eq 'fail'} @orthologyFilterHash{@orthologyList});
			if($orthologyPassCount > 0) {
				my $coverage = $orthologyPassCount / $orthologyCount;
				my $counts = join(',', $orthologyPassCount, $passCount - $orthologyPassCount, $orthologyFailCount, $failCount - $orthologyFailCount);
				$R->run("p.value <- fisher.test(matrix(c($counts), 2), alternative = \"greater\")\$p.value");
				my $pvalue = $R->get("p.value");
				if(scalar(grep {defined($_) && $_ > 0} @foldchangeList) == $orthologyCount) {
					push(@operonsOrthologiesFoldchangeList, [$operons, $orthologies, $log2foldchange, $coverage, $pvalue]);
				}
				if(scalar(grep {defined($_) && $_ < 0} @foldchangeList) == $orthologyCount) {
					push(@operonsOrthologiesFoldchangeList, [$operons, $orthologies, $log2foldchange, $coverage, $pvalue]);
				}
			}
		}
	}
}
close($reader);
$R->stop();

my %operonsDefinitionHash = ();
{
	open(my $reader, "$fmapPath/FMAP_data/known_operon.definition.txt");
	while(my $line = <$reader>) {
		chomp($line);
		my ($operons, $definition) = split(/\t/, $line);
		$operonsDefinitionHash{$operons} = $definition;
	}
	close($reader);
}

my %operonsOrthologyDefinitionHash = ();
{
	open(my $reader, "$fmapPath/FMAP_data/known_operon.KEGG_orthology_definition.txt");
	while(my $line = <$reader>) {
		chomp($line);
		my ($operons, $orthologyDefinition) = split(/\t/, $line);
		$operonsOrthologyDefinitionHash{$operons} = $orthologyDefinition;
	}
	close($reader);
}

my %operonsPathwaysHash = ();
{
	open(my $reader, "$fmapPath/FMAP_data/known_operon.KEGG_pathway.txt");
	while(my $line = <$reader>) {
		chomp($line);
		my ($operons, $pathways) = split(/\t/, $line);
		$operonsPathwaysHash{$operons} = $pathways;
	}
	close($reader);
}

print join("\t", 'operons', 'definition', 'log2foldchange', 'coverage', 'pvalue', 'orthologies', 'pathways'), "\n";
foreach(@operonsOrthologiesFoldchangeList) {
	my ($operons, $orthologies, $log2foldchange, $coverage, $pvalue) = @$_;
	my $orthologyCount = scalar(my @orthologyList = split(/ /, $orthologies));
	my $definition = $operonsDefinitionHash{$operons};
	$definition = $operonsOrthologyDefinitionHash{$operons} unless(defined($definition));
	if(defined($definition)) {
		my $pathways = $operonsPathwaysHash{$operons};
		$pathways = '' unless(defined($pathways));
		print join("\t", $operons, $definition, $log2foldchange, $coverage, $pvalue, $orthologies, $pathways), "\n";
	}
}

sub mean {
	if(@_) {
		return sum(@_) / scalar(@_);
	} else {
		return 0;
	}
}
