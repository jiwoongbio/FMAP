# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Getopt::Long;
use Statistics::R;

GetOptions('h' => \(my $help = ''),
	't=s' => \(my $test = 'kruskal'),
	'f=f' => \(my $foldchangeCutoff = 2),
	'p=f' => \(my $pvalueCutoff = 0.05),
	'a=f' => \(my $padjustCutoff = 1));
my @availableTestList = ('kruskal', 'anova', 'poisson', 'quasipoisson', 'metagenomeSeq');
my $availableTests = join(', ', map {"\"$_\""} @availableTestList);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl FMAP_comparison.pl [options] abundance_table.txt control1[,control2[...]] case1[,case2[...]] [...] > orthology_test_stat.txt

Options: -h       display this help message
         -t STR   statistical test for comparing sample groups, $availableTests [$test]
         -f FLOAT fold change cutoff [$foldchangeCutoff]
         -p FLOAT p-value cutoff [$pvalueCutoff]
         -a FLOAT FDR-adjusted p-value cutoff [$padjustCutoff]

EOF
}
die "ERROR: The test is not provided.\n" if(scalar(grep {$test eq $_} @availableTestList) == 0);

my ($tableFile, @samplesList) = @ARGV;
my @sampleListList = map {[split(/,/, $_)]} @samplesList;
die "ERROR: The input \"$tableFile\" is not available.\n" unless(-r $tableFile);
die "ERROR: Not enough sample groups.\n" unless(scalar(@sampleListList) > 1);
die "ERROR: The comparison requires at least 3 samples per group.\n" unless(scalar(grep {scalar(@$_) < 3} @sampleListList) == 0);

my $R = Statistics::R->new();
$R->run("table <- data.frame()");

my @orthologyList = ();
my %orthologyDefinitionHash = ();
{
	open(my $reader, $tableFile);
	chomp(my $line = <$reader>);
	my @columnList = split(/\t/, $line);
	while(my $line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		@tokenHash{@columnList} = split(/\t/, $line);
		push(@orthologyList, my $orthology = $tokenHash{'orthology'});
		$orthologyDefinitionHash{$orthology} = $_ if(defined($_ = $tokenHash{'definition'}));
		my $values = join(',', map {@tokenHash{@$_}} @sampleListList);
		$R->run("table <- rbind(table, $orthology = c($values))");
	}
	close($reader);
}

$R->set('condition', [map {($_) x scalar(@{$sampleListList[$_]})} (0 .. $#sampleListList)]);
if($test eq 'kruskal') {
	$R->run("p.value <- apply(table, 1, function(x) {kruskal.test(x, condition)\$p.value})");
}
if($test eq 'anova') {
	$R->run("p.value <- apply(log2(table + 1), 1, function(x) {summary(aov(x ~ condition))[[1]][\"condition\", \"Pr(>F)\"]})");
}
if($test eq 'poisson') {
	$R->run("offset.log.colSums <- offset(log(colSums(table)))");
	$R->run("p.value <- apply(table, 1, function(x) {anova(glm(round(x) ~ 1 + offset.log.colSums + condition, family = \"poisson\"), test = \"Chisq\")[\"condition\", \"Pr(>Chi)\"]})");
}
if($test eq 'quasipoisson') {
	$R->run("offset.log.colSums <- offset(log(colSums(table)))");
	$R->run("p.value <- apply(table, 1, function(x) {anova(glm(round(x) ~ 1 + offset.log.colSums + condition, family = \"quasipoisson\"), test = \"Chisq\")[\"condition\", \"Pr(>Chi)\"]})");
}
if($test eq 'metagenomeSeq') {
	$R->run("library(metagenomeSeq)");

	$R->run("phenoData <- AnnotatedDataFrame(data.frame(condition = condition, row.names = colnames(table)))");
	$R->run("featureData <- AnnotatedDataFrame(data.frame(orthology = rownames(table), row.names = rownames(table)))");
	$R->run("MRexperiment <- newMRexperiment(table, phenoData = phenoData, featureData = featureData)");

	$R->run("MRexperiment.p <- cumNormStat(MRexperiment, pFlag = TRUE)");
	$R->run("MRexperiment <- cumNorm(MRexperiment, p = MRexperiment.p)");

	$R->run("normFactor <- normFactors(MRexperiment)");
	$R->run("normFactor <- log2(normFactor / median(normFactor) + 1)");
	$R->run("mod <- model.matrix(~ condition + normFactor)");
	$R->run("fit <- fitZig(obj = MRexperiment, mod = mod)");
	$R->run("MRcoefs <- MRcoefs(fit, number = nrow(assayData(MRexperiment)\$counts))");

	$R->run("table <- MRcounts(MRexperiment, norm = TRUE, log = TRUE)[rownames(table), colnames(table)]");
	$R->run("p.value <- MRcoefs[rownames(table), \"pvalues\"]");
	$R->run("p.value[is.na(p.value)] <- 1");
}

my $log2foldchangeList;
if(scalar(@samplesList) == 2) {
	$log2foldchangeList = $R->get("as.numeric(rowMeans(log2(table[, condition != 0] + 1)) - rowMeans(log2(table[, condition == 0] + 1)))");
} else {
	$R->run("log2means <- t(apply(table, 1, function(x) {tapply(x, condition, function(y) {mean(log2(y + 1))})}))");
	$R->run("colnames(log2means) <- paste(\"condition\", colnames(log2means), sep = \"\")");
	$log2foldchangeList = $R->get("as.numeric(apply(abs(apply(combn(paste(\"condition\", unique(condition), sep = \"\"), 2), 2, function(x) {log2means[, x[1]] - log2means[, x[2]]})), 1, min))");
}
my $pvalueList = $R->get("as.numeric(p.value)");
my $padjustList = $R->get("as.numeric(p.adjust(p.value, method = \"fdr\"))");
$R->stop();

die "ERROR: R.\n" if(scalar(grep {scalar(@orthologyList) != scalar(@$_)} ($log2foldchangeList, $pvalueList, $padjustList)) > 0);

if(scalar(keys %orthologyDefinitionHash) > 0) {
	print join("\t", 'orthology', 'definition', 'log2foldchange', 'p.value', 'p.adjust', 'filter'), "\n";
	foreach my $index (0 .. $#orthologyList) {
		my $orthology = $orthologyList[$index];
		my $definition = defined($_ = $orthologyDefinitionHash{$orthology}) ? $_ : '';
		my ($log2foldchange, $pvalue, $padjust) = map {$_->[$index]} ($log2foldchangeList, $pvalueList, $padjustList);
		my $filter = 'fail';
		$filter = 'pass' if(2 ** abs($log2foldchange) > $foldchangeCutoff && $pvalue < $pvalueCutoff && $padjust < $padjustCutoff);
		print join("\t", $orthology, $definition, $log2foldchange, $pvalue, $padjust, $filter), "\n";
	}
} else {
	print join("\t", 'orthology', 'log2foldchange', 'p.value', 'p.adjust', 'filter'), "\n";
	foreach my $index (0 .. $#orthologyList) {
		my $orthology = $orthologyList[$index];
		my ($log2foldchange, $pvalue, $padjust) = map {$_->[$index]} ($log2foldchangeList, $pvalueList, $padjustList);
		my $filter = 'fail';
		$filter = 'pass' if(2 ** abs($log2foldchange) > $foldchangeCutoff && $pvalue < $pvalueCutoff && $padjust < $padjustCutoff);
		print join("\t", $orthology, $log2foldchange, $pvalue, $padjust, $filter), "\n";
	}
}
