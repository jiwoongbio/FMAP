# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die "ERROR in $0: ", $_[0] };

use Getopt::Long qw(:config no_ignore_case);
use Statistics::R;

my %parameterHash = ('width' => 12, 'height' => 8, 'leftmargin' => 20);
GetOptions('h' => \(my $help = ''),
	'W|width=i' => \$parameterHash{'width'},
	'H|height=i' => \$parameterHash{'height'},
	'L|leftmargin=f' => \$parameterHash{'leftmargin'},
	'p|pvalue=f' => \(my $pvalueCutoff = 0.05),
	'c|coverage=f' => \(my $coverageCutoff = -1),
	'd' => \(my $doNotPrintDefinition = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl FMAP_plot.pl [options] pathway.txt|module.txt|operon.txt plot.pdf

Options: -h       display this help message
         -W INT   plot width [$parameterHash{'width'}]
         -H INT   plot height [$parameterHash{'height'}]
         -L FLOAT plot left margin [$parameterHash{'leftmargin'}]
         -p FLOAT p-value cutoff [$pvalueCutoff]
         -c FLOAT coverage cutoff [0 for pathway, 1 for module and operons]
         -d       do not print definition

EOF
}
my ($inputFile, $pdfFile) = @ARGV;
foreach($inputFile) {
	die "ERROR in $0: '$_' is not readable.\n" unless(-r $_);
}
foreach($pdfFile) {
	my $directory = /^(.*\/)/ ? $1 : '.';
	die "ERROR in $0: '$directory' is not a writable directory.\n" unless(-d $directory && -w $directory);
}

my %pvalueHash = ();
my %updownRatioHash = ();
{
	open(my $reader, $inputFile);
	chomp(my $line = <$reader>);
	my @columnList = split(/\t/, $line);
	my %columnHash = map {$_ => 1} @columnList;
	if($coverageCutoff < 0) {
		if(defined($columnHash{'pathway'})) {
			$coverageCutoff = 0;
		} elsif(defined($columnHash{'module'})) {
			$coverageCutoff = 1;
		} elsif(defined($columnHash{'operons'})) {
			$coverageCutoff = 1;
		}
	}
	while(my $line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		@tokenHash{@columnList} = split(/\t/, $line);
		my ($name) = grep {defined} @tokenHash{'pathway', 'module', 'operons'};
		if($doNotPrintDefinition eq '') {
			my $definition = $tokenHash{'definition'};
			$name = "$name $definition" if(defined($definition) && $definition ne '');
		}
		if($tokenHash{'pvalue'} <= $pvalueCutoff && $tokenHash{'coverage'} >= $coverageCutoff) {
			$pvalueHash{$name} = $tokenHash{'pvalue'};
			if(defined($tokenHash{'orthology.count'}) && $tokenHash{'orthology.count'} =~ /^([0-9]+) \((.*)\)$/) {
				$updownRatioHash{$name}->{$_->[0]} = $_->[1] / $1 foreach(grep {$_->[0] eq 'up' || $_->[0] eq 'down'} map {[split(/:/, $_, 2)]} split(/, /, $2));
			}
		}
	}
	close($reader);
}
if(%updownRatioHash) {
	my $R = Statistics::R->new();
	$R->run('height <- c()');
	foreach(sort {$b->[1] <=> $a->[1]} map {[$_, $pvalueHash{$_}]} keys %pvalueHash) {
		my ($name, $pvalue) = @$_;
		my $height = -(log($pvalue) / log(10));
		$R->run(sprintf('height <- cbind(height, "%s" = c(%f, %f))', $name, map {defined($_) ? $_ * $height : 0} @{$updownRatioHash{$name}}{'up', 'down'}));
	}
	$R->run(sprintf('pdf(file = "%s", width = %d, height = %d)', $pdfFile, @parameterHash{'width', 'height'}));
	$R->run(sprintf('par(mar = c(par()$mar[1], par()$mar[2] + %f, par()$mar[3], par()$mar[4]))', $parameterHash{'leftmargin'}));
	$R->run('barplot(height, horiz = TRUE, border = NA, xlab = "-log10(p-value)", las = 1, col = c("red", "blue"))');
	$R->run('dev.off()');
	$R->stop();
} else {
	my $R = Statistics::R->new();
	$R->run('height <- c()');
	foreach(sort {$b->[1] <=> $a->[1]} map {[$_, $pvalueHash{$_}]} keys %pvalueHash) {
		my ($name, $pvalue) = @$_;
		$R->run(sprintf('height["%s"] <- %f', $name, -(log($pvalue) / log(10))));
	}
	$R->run(sprintf('pdf(file = "%s", width = %d, height = %d)', $pdfFile, @parameterHash{'width', 'height'}));
	$R->run(sprintf('par(mar = c(par()$mar[1], par()$mar[2] + %f, par()$mar[3], par()$mar[4]))', $parameterHash{'leftmargin'}));
	$R->run('barplot(height, horiz = TRUE, border = NA, xlab = "-log10(p-value)", las = 1)');
	$R->run('dev.off()');
	$R->stop();
}
