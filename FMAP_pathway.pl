# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die "ERROR in $0: ", $_[0] };

use Cwd 'abs_path';
use Getopt::Long;
use IPC::Open2;
use Statistics::R;

(my $fmapPath = abs_path($0)) =~ s/\/[^\/]*$//;
my $dataPath = "$fmapPath/FMAP_data";
system("mkdir -p $dataPath");

GetOptions(
	'h' => \(my $help = ''),
	'r' => \(my $redownload = ''),
	'p=s' => \(my $orthology2pathwayFile = "$dataPath/KEGG_orthology2pathway.txt"),
	'd=s' => \(my $pathwayDefinitionFile = "$dataPath/KEGG_pathway.txt"),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl FMAP_pathway.pl [options] orthology_test_stat.txt > pathway.txt

Options: -h       display this help message
         -r       redownload data

EOF
}

if(not -r $orthology2pathwayFile or $redownload) {
	open(my $reader, 'wget --no-verbose -O - http://rest.kegg.jp/link/pathway/ko |');
	open(my $writer, "> $orthology2pathwayFile");
	while(my $line = <$reader>) {
		chomp($line);
		my ($orthology, $pathway) = split(/\t/, $line, -1);
		$orthology =~ s/^ko://;
		$pathway =~ s/^path://;
		if($pathway =~ /^map[0-9]+$/) {
			print $writer join("\t", $orthology, $pathway), "\n";
		}
	}
	close($reader);
	close($writer);
}
my %orthologyHash = ();
my %pathwayOrthologyHash = ();
{
	open(my $reader, $orthology2pathwayFile);
	while(my $line = <$reader>) {
		chomp($line);
		my ($orthology, $pathway) = split(/\t/, $line, -1);
		$orthologyHash{$orthology} = 1;
		$pathwayOrthologyHash{$pathway}->{$orthology} = 1;
	}
	close($reader);
}
my $orthologyCount = scalar(keys %orthologyHash);

if(not -r $pathwayDefinitionFile or $redownload) {
	open(my $reader, 'wget --no-verbose -O - http://rest.kegg.jp/list/pathway |');
	open(my $writer, "> $pathwayDefinitionFile");
	while(my $line = <$reader>) {
		chomp($line);
		my ($pathway, $definition) = split(/\t/, $line, -1);
		$pathway =~ s/^path://;
		if($pathway =~ /^map[0-9]+$/) {
			print $writer join("\t", $pathway, $definition), "\n";
		}
	}
	close($reader);
	close($writer);
}
my %pathwayDefinitionHash = ();
{
	open(my $reader, $pathwayDefinitionFile);
	while(my $line = <$reader>) {
		chomp($line);
		my ($pathway, $definition) = split(/\t/, $line, -1);
		$pathwayDefinitionHash{$pathway} = $definition;
	}
	close($reader);
}

my ($inputFile) = @ARGV;
#die "ERROR in $0: '$_' is not readable.\n" unless(-r $inputFile);

my %targetOrthologyRegulationHash = ();
{
	open(my $reader, $inputFile);
	my @columnList = ();
	while(my $line = <$reader>) {
		chomp($line);
		$line =~ s/\r*\n*$//;
		my @tokenList = split(/\t/, $line, -1);
		unless(@columnList) {
			if(grep {$_ eq 'orthology'} @tokenList) {
				@columnList = @tokenList;
				next;
			} else {
				@columnList = ('orthology');
			}
		}
		my %tokenHash = ();
		@tokenHash{@columnList} = @tokenList;
		next if(defined($tokenHash{'filter'}) && lc($tokenHash{'filter'}) ne 'pass');
		my $orthology = $tokenHash{'orthology'};
		if($orthologyHash{$orthology}) {
			$targetOrthologyRegulationHash{$orthology} = 0;
			if(defined(my $regulation = $tokenHash{'regulation'})) {
				$targetOrthologyRegulationHash{$orthology} = $regulation;
			} elsif(defined(my $log2foldchange = $tokenHash{'log2foldchange'})) {
				$targetOrthologyRegulationHash{$orthology} = $log2foldchange;
			}
		}
	}
	close($reader);
}
my $targetOrthologyCount = scalar(keys %targetOrthologyRegulationHash);

print join("\t", 'pathway', 'definition', 'orthology.count', 'coverage', 'regulation', 'pvalue', 'weblink'), "\n";
open(my $writer, "| sort -t '\t' -k6,6g -k4,4gr -k3,3nr -k1,1");
my $R = Statistics::R->new();
foreach my $pathway (sort keys %pathwayOrthologyHash) {
	my @pathwayOrthologyList = sort keys %{$pathwayOrthologyHash{$pathway}};
	my $pathwayOrthologyCount = scalar(@pathwayOrthologyList);
	my @pathwayTargetOrthologyRegulationList = ();
	foreach my $orthology (@pathwayOrthologyList) {
		if(defined(my $regulation = $targetOrthologyRegulationHash{$orthology})) {
			push(@pathwayTargetOrthologyRegulationList, [$orthology, $regulation]);
		}
	}
	my $pathwayTargetOrthologyCount = scalar(@pathwayTargetOrthologyRegulationList);
	next if($pathwayTargetOrthologyCount == 0);
	my $definition = $pathwayDefinitionHash{$pathway};
	$definition = '' unless(defined($definition));
	my $coverage = $pathwayTargetOrthologyCount / $pathwayOrthologyCount;
	my %regulationCountHash = ();
	foreach(@pathwayTargetOrthologyRegulationList) {
		my ($orthology, $regulation) = @$_;
		$regulationCountHash{'up'} += 1 if($regulation > 0);
		$regulationCountHash{'down'} += 1 if($regulation < 0);
	}
	my $regulation = join(', ', map {defined($_->[1]) ? join(':', @$_) : ()} map {[$_, $regulationCountHash{$_}]} ('up', 'down'));
	my $counts = join(',', $pathwayTargetOrthologyCount, $pathwayOrthologyCount - $pathwayTargetOrthologyCount, $targetOrthologyCount - $pathwayTargetOrthologyCount, $orthologyCount - $pathwayOrthologyCount - $targetOrthologyCount + $pathwayTargetOrthologyCount);
	$R->run("p.value <- fisher.test(matrix(c($counts), 2), alternative = \"greater\")\$p.value");
	my $pvalue = $R->get('p.value');
	my $weblink = sprintf("http://www.kegg.jp/kegg-bin/show_pathway?map=%s&multi_query=%s", $pathway, join('%0d%0a', map {join('+', $_->[0], ($_->[1] > 0) ? "magenta" : ($_->[1] < 0) ? "cyan" : "yellow")} @pathwayTargetOrthologyRegulationList));
	print $writer join("\t", $pathway, $definition, $pathwayTargetOrthologyCount, $coverage, $regulation, $pvalue, $weblink), "\n";
}
$R->stop();
close($writer);
