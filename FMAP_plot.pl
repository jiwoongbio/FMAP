# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
use Getopt::Long;
use Statistics::R;

my %parameterHash = ('width' => 1200, 'height' => 800, 'pointsize' => 20, 'leftmargin' => 20);
GetOptions('h' => \(my $help = ''),
	'w|width=i' => \$parameterHash{'width'},
	'h|height=i' => \$parameterHash{'height'},
	'p|pointsize=i' => \$parameterHash{'pointsize'},
	'l|leftmargin=f' => \$parameterHash{'leftmargin'},
	'c|cutoff=f' => \(my $cutoff = 0.05),
	'd' => \(my $doNotPrintDefinition = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl FMAP_plot.pl [options] pathway.txt|module.txt|operon.txt plot.png

Options: -h       display this help message
         -w INT   plot width [$parameterHash{'width'}]
         -h INT   plot height [$parameterHash{'height'}]
         -p INT   plot point size [$parameterHash{'pointsize'}]
         -l FLOAT plot left margin [$parameterHash{'leftmargin'}]
         -c FLOAT p-value cutoff [$cutoff]
         -d       do not print definition

EOF
}
my ($inputFile, $pngFile) = @ARGV;
my %pvalueHash = ();
{
	open(my $reader, $inputFile);
	chomp(my $line = <$reader>);
	my @columnList = split(/\t/, $line);
	while(my $line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		@tokenHash{@columnList} = split(/\t/, $line);
		my ($name) = grep {defined} @tokenHash{'pathway', 'module', 'operons'};
		if($doNotPrintDefinition eq '') {
			my $definition = $tokenHash{'definition'};
			$name = "$name $definition" if(defined($definition) && $definition ne '');
		}
		$pvalueHash{$name} = $_ if(($_ = $tokenHash{'pvalue'}) <= $cutoff);
	}
	close($reader);
}
{
	my $R = Statistics::R->new();
	$R->run('counts <- c()');
	foreach(sort {$b->[1] <=> $a->[1]} map {[$_, $pvalueHash{$_}]} keys %pvalueHash) {
		my ($name, $pvalue) = @$_;
		$R->run(sprintf('counts["%s"] <- %f', $name, -(log($pvalue) / log(10))));
	}
	$R->run(sprintf('png(filename = "%s", width = %d, height = %d, pointsize = %d)', $pngFile, @parameterHash{'width', 'height', 'pointsize'}));
	$R->run(sprintf('par(mar = c(par()$mar[1], par()$mar[2] + %f, par()$mar[3], par()$mar[4]))', $parameterHash{'leftmargin'}));
	$R->run('barplot(counts, horiz = TRUE, border = NA, xlab = "-log10(p-value)", las = 1)');
	$R->run('dev.off()');
	$R->stop();
}
