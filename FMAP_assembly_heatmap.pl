# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die "ERROR in $0: ", $_[0] };

use Cwd 'abs_path';
use Getopt::Long qw(:config no_ignore_case);
use Statistics::R;
use Bio::DB::Taxonomy;

(my $fmapPath = abs_path($0)) =~ s/\/[^\/]*$//;
my $database = loadDefaultDatabase();
my $databasePrefix = "$fmapPath/FMAP_data/$database";
my $dataPath = "$fmapPath/FMAP_data";
system("mkdir -p $dataPath");

if(not -r "$dataPath/nodes.dmp" or not -r "$dataPath/names.dmp") {
	my $URL = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz";
	my $file = "$dataPath/taxdump.tar.gz";
	system("wget --no-verbose -O $file $URL") if(not -r $file);
	die "ERROR in $0: '$file' has zero size.\n" if(-z $file);
	system("cd $dataPath; tar -zxf taxdump.tar.gz nodes.dmp");
	system("cd $dataPath; tar -zxf taxdump.tar.gz names.dmp");
	system("rm -f $dataPath/$_") foreach('nodes', 'parents', 'names2id', 'id2names');
}
my $db = Bio::DB::Taxonomy->new(-source => 'flatfile', -directory => $dataPath, -nodesfile => "$dataPath/nodes.dmp", -namesfile => "$dataPath/names.dmp");

GetOptions('h' => \(my $help = ''),
	'c=s' => \(my $comparisonFile = ''),
	'f' => \(my $fontSize = 15),
	'w' => \(my $tdWidth = 60),
	'orthology2definition=s' => \(my $orthologyDefinitionFile = ''),
	'd=s' => \$databasePrefix);

if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl FMAP_assembly_heatmap.pl [options] [name=]FMAP_assembly.abundance.txt [...] > FMAP_assembly_heatmap.html

Options: -h       display this help message
         -c FILE  comparison output file including orthology and filter columns
         -f INT   HTML font size
         -w INT   HTML table cell width

EOF
}

my %orthologyFilterHash = ();
if($comparisonFile ne '') {
	die "ERROR in $0: '$comparisonFile' is not readable.\n" unless(-r $comparisonFile);
	open(my $reader, $comparisonFile);
	chomp(my $line = <$reader>);
	my @columnList = split(/\t/, $line);
	while(my $line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		@tokenHash{@columnList} = split(/\t/, $line);
		my ($orthology, $filter) = @tokenHash{'orthology', 'filter'};
		$orthologyFilterHash{$orthology} = $filter;
	}
	close($reader);
}

if($databasePrefix eq "$fmapPath/FMAP_data/$database") {
	if($database =~ /^orthology_uniref/) {
		$orthologyDefinitionFile = "$fmapPath/FMAP_data/KEGG_orthology.txt" if($orthologyDefinitionFile eq '');
	} elsif($database =~ /^ARDB/ || $database =~ /^betalactamases/) {
		$orthologyDefinitionFile = "$databasePrefix.definition.txt" if($orthologyDefinitionFile eq '');
	}
}

my %orthologyDefinitionHash = ();
if($orthologyDefinitionFile ne '') {
	die "ERROR in $0: '$orthologyDefinitionFile' is not readable.\n" unless(-r $orthologyDefinitionFile);
	open(my $reader, $orthologyDefinitionFile);
	while(my $line = <$reader>) {
		chomp($line);
		my ($orthology, $definition) = split(/\t/, $line);
		$orthologyDefinitionHash{$orthology} = $definition;
	}
	close($reader);
}

my @sampleFileList = @ARGV;
my @sampleList = ();
my %orthologySampleAbundanceHash = ();
my %orthologySampleTaxonAbundanceHash = ();
foreach my $sampleFile (@sampleFileList) {
	($sampleFile, my $sample) = reverse(split(/=/, $sampleFile, 2));
	$sample = $sampleFile unless(defined($sample));
	push(@sampleList, $sample);
	open(my $reader, $sampleFile);
	chomp(my $line = <$reader>);
	my @columnList = split(/\t/, $line);
	while(my $line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		@tokenHash{@columnList} = split(/\t/, $line);
		my $orthology = $tokenHash{'orthology'};
		my $abundance = $tokenHash{'meanDepth/genome'};
		next if(defined($_ = $orthologyFilterHash{$orthology}) && $_ eq 'fail');
		$orthologySampleAbundanceHash{$orthology}->{$sample} += $abundance;
		$orthologySampleTaxonAbundanceHash{$orthology}->{$sample}->{$_} += $abundance if(defined($_ = $tokenHash{'taxon'}));
		$orthologyDefinitionHash{$orthology} = $_ if(defined($_ = $tokenHash{'definition'}));
	}
	close($reader);
}

my $R = Statistics::R->new();
$R->run('table <- data.frame()');
foreach my $orthology (keys %orthologySampleAbundanceHash) {
	my $abundances = join(',', map {defined($_) ? $_ : 0} map {$orthologySampleAbundanceHash{$orthology}->{$_}} @sampleList);
	$R->run("table <- rbind(table, $orthology = c($abundances))");
}
$R->run('hc <- hclust(dist(table))');
my $orthologyList = $R->get('hc$labels[hc$order]');
$R->stop();

#print <<EOF;
#Content-Type: text/html; charset=utf-8
#
#EOF

print <<EOF;
<!DOCTYPE html>
<html>
<head>
<style>
* {
	font-size: ${fontSize}px;
	line-height: ${fontSize}px;
}
table.heatmap {
	border-spacing: 1px;
}
table.heatmap td {
	padding: 5px;
	white-space: nowrap;
	border: 2px solid #FFFFFF;
}
table.heatmap td div.popup {
	text-align: left;
	display: none;
	position: absolute;
	background-color: #FFFFFF;
	padding: 5px;
	border: 2px solid #000000;
}
table.heatmap td:hover div.popup {
	display: block;
}
</style>
</head>
<body>
<form method="post">
EOF

print "<p><table class=\"heatmap\">\n";
{
	my @tdList = ("<td></td>");
	foreach my $sample (@sampleList) {
		push(@tdList, "<td style=\"text-align: center; width: ${tdWidth}px;\">$sample</td>");
	}
	print join('', "<tr>", @tdList, "</tr>"), "\n";
}
foreach my $orthology (@$orthologyList) {
	my @tdList = ();
	if(defined($_ = $orthologyDefinitionHash{$orthology}) && $_ ne '') {
		(my $description = $_) =~ s/;/<br>/g;
		push(@tdList, "<td>$orthology<div class=\"popup\">$description</div></td>");
	} else {
		push(@tdList, "<td>$orthology</td>");
	}
	foreach my $sample (@sampleList) {
		my $abundance = defined($_ = $orthologySampleAbundanceHash{$orthology}->{$sample}) ? $_ : 0;
		my $color = '#FFFFFF';
		if($abundance > 1) {
			$color = color(1, 1 / $abundance, 0);
		} else {
			$color = color($abundance, 1, 0);
		}
		if(defined($_ = $orthologySampleTaxonAbundanceHash{$orthology}->{$sample})) {
			my %taxonAbundanceHash = %$_;
			my @taxonAbundanceList = sort {$b->[1] <=> $a->[1]} map {[$_ =~ /^[0-9]+$/ ? $db->get_taxon(-taxonid => $_)->scientific_name : $_, $taxonAbundanceHash{$_}]} keys %taxonAbundanceHash;
			my $taxonAbundances = join('<br>', map {sprintf('%s %.3f', @$_)} @taxonAbundanceList);
			my @taxonFontSizeList = map {[substr($_->[0], 0, 4), $fontSize * ($_->[1] / $abundance)]} @taxonAbundanceList;
			my $taxons = join('', map {sprintf('<div style="font-size: %fpx; line-height: %fpx;">%s</div>', $_->[1], $_->[1], $_->[0])} @taxonFontSizeList);
			push(@tdList, "<td style=\"text-align: center; background-color: $color;\">$taxons<div class=\"popup\">$taxonAbundances</div></td>");
		} else {
			push(@tdList, "<td style=\"text-align: center; background-color: $color;\"></td>");
		}
	}
	print join('', "<tr>", @tdList, "</tr>"), "\n";
}
print "</table></p>\n";

print "<p><b>Color gradation for relative depth</b></p>\n";
print "<p><table class=\"heatmap\">\n";
for(my $depth = 10; $depth > 0; $depth -= 2) {
	my $color = color(1, 1 / $depth, 0);
	print "<tr><td style=\"text-align: center; width: ${tdWidth}px; background-color: $color;\">$depth</td></tr>\n";
}
for(my $depth = 1; $depth >= 0; $depth = sprintf('%.1f', $depth - 0.2)) {
	my $color = color($depth, 1, 0);
	print "<tr><td style=\"text-align: center; width: ${tdWidth}px; background-color: $color;\">$depth</td></tr>\n";
}
print "</table></p>\n";

print <<EOF;
</form>
</body>
</html>
EOF

sub loadDefaultDatabase {
	open(my $reader, "$fmapPath/FMAP_data/database");
	chomp(my $line = <$reader>);
	close($reader);
	return $1 if($line =~ /^([A-Za-z0-9_.]+)/);
}

sub color {
	my ($red, $green, $blue) = @_;
	(my $color = sprintf('#%02x%02x%02x', int($red * 255), int($green * 255), int($blue * 255))) =~ tr/a-z/A-Z/;
	return $color;
}
