#!/usr/bin/env perl
use strict;
use warnings;
use Cwd 'abs_path';
use Getopt::Long;
use Statistics::R;

(my $fmapPath = abs_path($0)) =~ s/\/[^\/]*$//;

GetOptions( 'h' => \(my $help = ''),
			'd=s' => \(my $databasePath = ''));
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl FMAP_pathway.pl [options] orthology_test_stat.txt > pathway.txt

Options: 	-h       display this help message
			-d       full path for database location [$fmapPath/FMAP_data]

EOF
}

unless ($databasePath ne '') {
	$databasePath = $fmapPath . "/FMAP_data";
};
die "Database path does not exist: " . $databasePath . "\n" unless (-d $databasePath);

my ($inputFile) = @ARGV;
die "ERROR: The input \"$inputFile\" is not available.\n" unless(-r $inputFile);
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
my $totalTargetOrthologyCount = scalar(keys %targetOrthologyColorHash);

my %pathwayTargetOrthologyListHash = ();
my %pathwayTotalOrthologyCountHash = ();
my %orthologyHash = ();

$databasePath =~ s/\/$//;
die "Could not find " . $databasePath . "/KEGG_orthology2pathway.txt\n" unless (-r $databasePath . "/KEGG_orthology2pathway.txt");
open(my $reader, $databasePath . "/KEGG_orthology2pathway.txt");
while(my $line = <$reader>) {
	chomp($line);
	my ($orthology, $pathway) = split(/\t/, $line);
	my $color = $targetOrthologyColorHash{$orthology};
	push(@{$pathwayTargetOrthologyListHash{$pathway}}, "$orthology $color") if(defined($color));
	$pathwayTotalOrthologyCountHash{$pathway} += 1;
	$orthologyHash{$orthology} = 1;
}
close($reader);
my $totalOrthologyCount = scalar(keys %orthologyHash);

my %pathwayDefinitionHash = ();
{
	die "Could not find " . $databasePath . "/KEGG_pathway.txt\n" unless (-r $databasePath . "/KEGG_pathway.txt");
	open(my $reader, $databasePath . "/KEGG_pathway.txt");
	while(my $line = <$reader>) {
		chomp($line);
		my ($pathway, $definition) = split(/\t/, $line);
		$pathwayDefinitionHash{$pathway} = $definition;
	}
	close($reader);
}

print join("\t", 'pathway', 'definition', 'orthology.count', 'coverage', 'pvalue', 'orthology.colors'), "\n";
my $R = Statistics::R->new();
foreach my $pathway (sort keys %pathwayTargetOrthologyListHash) {
	my $pathwayTargetOrthologyCount = scalar(my @pathwayTargetOrthologyList = @{$pathwayTargetOrthologyListHash{$pathway}});
	my $pathwayTotalOrthologyCount = $pathwayTotalOrthologyCountHash{$pathway};
	my $counts = join(',', $pathwayTargetOrthologyCount, $pathwayTotalOrthologyCount - $pathwayTargetOrthologyCount, $totalTargetOrthologyCount - $pathwayTargetOrthologyCount, $totalOrthologyCount - $pathwayTotalOrthologyCount - $totalTargetOrthologyCount + $pathwayTargetOrthologyCount);
	$R->run("p.value <- fisher.test(matrix(c($counts), 2), alternative = \"greater\")\$p.value");
	my $pvalue = $R->get("p.value");
	my $definition = $pathwayDefinitionHash{$pathway};
	$definition = '' unless(defined($definition));
	my $coverage = $pathwayTargetOrthologyCount / $pathwayTotalOrthologyCount;
	print join("\t", $pathway, $definition, $pathwayTargetOrthologyCount, $coverage, $pvalue, join('|', @pathwayTargetOrthologyList)), "\n";
}
$R->stop();
