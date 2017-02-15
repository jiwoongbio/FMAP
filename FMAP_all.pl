# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
use Cwd 'abs_path';
use Getopt::Long;

(my $fmapPath = abs_path($0)) =~ s/\/[^\/]*$//;

GetOptions('h' => \(my $help = ''),
	's' => \(my $scriptOnly = ''),
	'm=s' => \(my $mapperPath = 'diamond'),
	't=i' => \(my $mappingThreads = 1),
	'c=s' => \(my $test = 'kruskal'),
	'f=f' => \(my $foldchangeCutoff = 2),
	'p=f' => \(my $pvalueCutoff = 0.05),
	'a=f' => \(my $padjustCutoff = 1));
my @availableTestList = ('kruskal', 'anova', 'poisson', 'quasipoisson', 'metagenomeSeq');
my $availableTests = join(', ', map {"\"$_\""} @availableTestList);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl FMAP_all.pl [options] input.config [output_prefix]

Options: -h       display this help message
         -s       generate a script, but not execute it
         -m FILE  mapping: executable file path of mapping program, "diamond" or "usearch" [$mapperPath]
         -t INT   mapping: number of threads [$mappingThreads]
         -c STR   comparison: statistical test for comparing sample groups, $availableTests [$test]
         -f FLOAT comparison: fold change cutoff [$foldchangeCutoff]
         -p FLOAT comparison: p-value cutoff [$pvalueCutoff]
         -a FLOAT comparison: FDR-adjusted p-value cutoff [$padjustCutoff]

EOF
}
(my $mapper = $mapperPath) =~ s/^.*\///;
die "ERROR: The mapper must be \"diamond\" or \"usearch\".\n" unless($mapper eq 'diamond' || $mapper eq 'usearch');
die "ERROR: The mapper is not executable.\n" unless(-x getCommandPath($mapperPath));

die "ERROR: The test is not provided.\n" if(scalar(grep {$test eq $_} ('kruskal', 'anova', 'poisson', 'quasipoisson')) == 0);

my ($configFile, $outputPrefix) = @ARGV;
die "ERROR: The configuration file is not available.\n" if(not -r $configFile);
$outputPrefix = 'FMAP_output' unless(defined($outputPrefix));

my %sampleInputFileListHash = ();
my %groupSampleListHash = ();
open(my $reader, $configFile);
while(my $line = <$reader>) {
	chomp($line);
	next if($line =~ /^#/);
	my ($group, $sample, @inputFileList) = split(/\t/, $line);
	push(@{$sampleInputFileListHash{$sample}}, @inputFileList);
	push(@{$groupSampleListHash{$group}}, $sample);
}
close($reader);
my @groupList = sort keys %groupSampleListHash;
@{$groupSampleListHash{$_}} = sort @{$groupSampleListHash{$_}} foreach(@groupList);

open(my $writer, "> $outputPrefix.sh");
print $writer "perl $fmapPath/FMAP_prepare.pl -m $mapperPath\n";
print $writer "\n";

foreach(grep {$_ eq 'control'} @groupList) {
	my @sampleList = @{$groupSampleListHash{$_}};
	print $writer "perl $fmapPath/FMAP_mapping.pl -m $mapperPath -p $mappingThreads @{$sampleInputFileListHash{$_}} > $_.mapping.txt\n" foreach(@sampleList);
}
foreach(grep {$_ ne 'control'} @groupList) {
	my @sampleList = @{$groupSampleListHash{$_}};
	print $writer "perl $fmapPath/FMAP_mapping.pl -m $mapperPath -p $mappingThreads @{$sampleInputFileListHash{$_}} > $_.mapping.txt\n" foreach(@sampleList);
}
print $writer "\n";

foreach(grep {$_ eq 'control'} @groupList) {
	my @sampleList = @{$groupSampleListHash{$_}};
	print $writer "perl $fmapPath/FMAP_quantification.pl $_.mapping.txt > $_.abundance.txt\n" foreach(@sampleList);
}
foreach(grep {$_ ne 'control'} @groupList) {
	my @sampleList = @{$groupSampleListHash{$_}};
	print $writer "perl $fmapPath/FMAP_quantification.pl $_.mapping.txt > $_.abundance.txt\n" foreach(@sampleList);
}
print $writer "\n";

print $writer "perl $fmapPath/FMAP_table.pl \\\n";
foreach(grep {$_ eq 'control'} @groupList) {
	my @sampleList = @{$groupSampleListHash{$_}};
	print $writer "$_=$_.abundance.txt \\\n" foreach(@sampleList);
}
foreach(grep {$_ ne 'control'} @groupList) {
	my @sampleList = @{$groupSampleListHash{$_}};
	print $writer "$_=$_.abundance.txt \\\n" foreach(@sampleList);
}
print $writer "> $outputPrefix.table.txt\n";
print $writer "\n";

my @samplesList = ();
foreach(grep {$_ eq 'control'} @groupList) {
	my @sampleList = @{$groupSampleListHash{$_}};
	push(@samplesList, join(',', @sampleList));
}
foreach(grep {$_ ne 'control'} @groupList) {
	my @sampleList = @{$groupSampleListHash{$_}};
	push(@samplesList, join(',', @sampleList));
}
print $writer "perl $fmapPath/FMAP_comparison.pl -t $test -f $foldchangeCutoff -p $pvalueCutoff -a $padjustCutoff $outputPrefix.table.txt @samplesList > $outputPrefix.comparison.txt\n";
print $writer "\n";

print $writer "perl $fmapPath/FMAP_operon.pl $outputPrefix.comparison.txt > $outputPrefix.operon.txt\n";
print $writer "perl $fmapPath/FMAP_pathway.pl $outputPrefix.comparison.txt > $outputPrefix.pathway.txt\n";
print $writer "perl $fmapPath/FMAP_module.pl $outputPrefix.comparison.txt > $outputPrefix.module.txt\n";
close($writer);

print "The script file \"$outputPrefix.sh\" was generated.\n";
system("sh $outputPrefix.sh") if($scriptOnly eq '');

sub getCommandPath {
	my ($command) = @_;
	chomp($command = `which $command`) if($command ne '' && $command !~ /\//);
	$command =~ s/^~\//$ENV{'HOME'}\//;
	$command = abs_path($command) if($command ne '');
	return $command;
}
