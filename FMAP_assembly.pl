# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Cwd 'abs_path';
use Getopt::Long;
use List::Util qw(max);

(my $fmapPath = abs_path($0)) =~ s/\/[^\/]*$//;
my $database = loadDefaultDatabase();
my $databasePrefix = "$fmapPath/FMAP_data/$database";
my @codonList = ();

GetOptions('h' => \(my $help = ''),
	'p=i' => \(my $threads = 1),
	'e=f' => \(my $evalue = 10),
	't=s' => \(my $temporaryDirectory = defined($ENV{'TMPDIR'}) ? $ENV{'TMPDIR'} : '/tmp'),
	'c=s' => \@codonList,
	's=s' => \(my $startCodons = 'GTG,ATG,CTG,TTG,ATA,ATC,ATT'),
	'l=i' => \(my $minimumTranslationLength = 10),
	'm=f' => \(my $minimumCoverage = 0.8),
	'q=i' => \(my $minimumMappingQuality = 0),
	's=s' => \(my $stranded = ''),
	'orthology2definition=s' => \(my $orthologyDefinitionFile = ''),
	'protein2orthology=s' => \(my $proteinOrthologyFile = ''),
	'd=s' => \$databasePrefix);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl FMAP_assembly.pl [options] output.prefix assembly.fasta input1.fastq|input1.R1.fastq,input1.R2.fastq [input2.fastq|input2.R1.fastq,input2.R2.fastq [...]] > summary.txt

Options: -h       display this help message
         -p INT   number of threads [$threads]
         -e FLOAT maximum e-value to report alignments for "diamond" [$evalue]
         -t DIR   directory for temporary files [\$TMPDIR or /tmp]
         -c STR   codon and translation e.g. ATG=M [NCBI genetic code 11 (Bacterial, Archaeal and Plant Plastid)]
         -s STR   start codons [$startCodons]
         -l INT   minimum translation length [$minimumTranslationLength]
         -m FLOAT minimum coverage [$minimumCoverage]
         -q INT   minimum mapping quality [$minimumMappingQuality]
         -s STR   strand specificity, "f" or "r"

EOF
}
my ($outputPrefix, $assemblyFastaFile, @fastqFileList) = @ARGV;

foreach('bwa', 'samtools', 'diamond') {
	die "ERROR: '$_' is not executable.\n" unless(-x getCommandPath($_));
}
foreach($assemblyFastaFile, "$databasePrefix.dmnd", "$databasePrefix.length.txt") {
	die "ERROR: '$_' is not readable.\n" unless(-r $_);
}
my $outputDirectory = ($outputPrefix =~ /^(.*\/)/) ? $1 : '.';
foreach($outputDirectory, $temporaryDirectory) {
	die "ERROR: '$_' is not a writable directory.\n" unless(-d $_ && -w $_);
}

$orthologyDefinitionFile = "$fmapPath/FMAP_data/KEGG_orthology.txt" if($databasePrefix eq "$fmapPath/FMAP_data/$database" && $orthologyDefinitionFile eq '');
die "ERROR: '$orthologyDefinitionFile' is not readable.\n" unless(-r $orthologyDefinitionFile);

my %orthologyDefinitionHash = ();
if($orthologyDefinitionFile ne '') {
	open(my $reader, $orthologyDefinitionFile);
	while(my $line = <$reader>) {
		chomp($line);
		my ($orthology, $definition) = split(/\t/, $line);
		$orthologyDefinitionHash{$orthology} = $definition;
	}
	close($reader);
}
my %proteinOrthologyHash = ();
if($proteinOrthologyFile ne '') {
	open(my $reader, $proteinOrthologyFile);
	while(my $line = <$reader>) {
		chomp($line);
		my ($protein, $orthology) = split(/\t/, $line);
		$proteinOrthologyHash{$protein} = $orthology;
	}
	close($reader);
}

{ # Read mapping
	system("bwa index -p $outputPrefix.bwa_index $assemblyFastaFile 1>&2");
	open(my $writer, "| samtools view -S -b - > $outputPrefix.bam");
	my $printHeader = 1;
	foreach my $fastqFile (@fastqFileList) {
		if($fastqFile =~ /^(.+),(.+)$/) {
			my ($fastqFile1, $fastqFile2) = ($1, $2);
			die "ERROR: '$fastqFile1' is not readable.\n" unless(-r $fastqFile1);
			die "ERROR: '$fastqFile2' is not readable.\n" unless(-r $fastqFile2);
			open(my $reader, "bwa mem -t $threads $outputPrefix.bwa_index $fastqFile1 $fastqFile2 |");
			while(my $line = <$reader>) {
				print $writer $line if($line !~ /^\@/ || $printHeader);
			}
			close($reader);
		} else {
			die "ERROR: '$fastqFile' is not readable.\n" unless(-r $fastqFile);
			open(my $reader, "bwa mem -t $threads $outputPrefix.bwa_index $fastqFile |");
			while(my $line = <$reader>) {
				print $writer $line if($line !~ /^\@/ || $printHeader);
			}
			close($reader);
		}
		$printHeader = '';
	}
	close($writer);
	system("samtools sort $outputPrefix.bam $outputPrefix.sorted") if(-r "$outputPrefix.bam");
	system("samtools index $outputPrefix.sorted.bam") if(-r "$outputPrefix.sorted.bam");
}
if(-r "$outputPrefix.sorted.bam" and -r "$outputPrefix.sorted.bam.bai") {
	chomp(my $readCount = `samtools view -c -F 2048 -q $minimumMappingQuality $outputPrefix.sorted.bam`);
	chomp(my $mappingReadCount = `samtools view -c -F 2052 -q $minimumMappingQuality $outputPrefix.sorted.bam`);
	my $mappingReadRatio = $mappingReadCount / $readCount;
	print "Mapping: $mappingReadCount / $readCount ($mappingReadRatio)\n";
} else {
	die "ERROR: Read mapping failed.\n";
}

my %contigSequenceLengthHash = ();
{ # ORF translation
	my %codonHash = (
		'TTT' => 'F', 'CTT' => 'L', 'ATT' => 'I', 'GTT' => 'V',
		'TTC' => 'F', 'CTC' => 'L', 'ATC' => 'I', 'GTC' => 'V',
		'TTA' => 'L', 'CTA' => 'L', 'ATA' => 'I', 'GTA' => 'V',
		'TTG' => 'L', 'CTG' => 'L', 'ATG' => 'M', 'GTG' => 'V',

		'TCT' => 'S', 'CCT' => 'P', 'ACT' => 'T', 'GCT' => 'A',
		'TCC' => 'S', 'CCC' => 'P', 'ACC' => 'T', 'GCC' => 'A',
		'TCA' => 'S', 'CCA' => 'P', 'ACA' => 'T', 'GCA' => 'A',
		'TCG' => 'S', 'CCG' => 'P', 'ACG' => 'T', 'GCG' => 'A',

		'TAT' => 'Y', 'CAT' => 'H', 'AAT' => 'N', 'GAT' => 'D',
		'TAC' => 'Y', 'CAC' => 'H', 'AAC' => 'N', 'GAC' => 'D',
		'TAA' => '*', 'CAA' => 'Q', 'AAA' => 'K', 'GAA' => 'E',
		'TAG' => '*', 'CAG' => 'Q', 'AAG' => 'K', 'GAG' => 'E',

		'TGT' => 'C', 'CGT' => 'R', 'AGT' => 'S', 'GGT' => 'G',
		'TGC' => 'C', 'CGC' => 'R', 'AGC' => 'S', 'GGC' => 'G',
		'TGA' => '*', 'CGA' => 'R', 'AGA' => 'R', 'GGA' => 'G',
		'TGG' => 'W', 'CGG' => 'R', 'AGG' => 'R', 'GGG' => 'G',
	);
	$codonHash{$_->[0]} = $_->[1] foreach(map {[split(/=/, $_)]} @codonList);
	sub translate {
		my ($sequence) = @_;
		return join('', map {defined($_) ? $_ : 'X'} map {$codonHash{substr($sequence, $_ * 3, 3)}} 0 .. (length($sequence) / 3) - 1);
	}
	my %startCodonHash = map {$_ => 1} split(/,/, $startCodons);
	my %stopCodonHash = map {$_ => 1} (
		'TAG', # amber
		'TAA', # ochre
		'TGA', # opal or umber
	);
	my $minimumLength = $minimumTranslationLength * 3;
	open(my $reader, $assemblyFastaFile);
	open(my $writer, "> $outputPrefix.translation.fasta");
	my ($contig, $sequence) = ('', '');
	while(my $line = <$reader>) {
		chomp($line);
		if($line =~ /^>(\S*)$/) {
			if($contig ne '' && $sequence ne '') {
				writeTranslationSequences($contig, $sequence);
				$contigSequenceLengthHash{$contig} = length($sequence);
			}
			($contig, $sequence) = ($1, '');
		} else {
			$sequence .= $line;
		}
	}
	if($contig ne '' && $sequence ne '') {
		writeTranslationSequences($contig, $sequence);
		$contigSequenceLengthHash{$contig} = length($sequence);
	}
	close($reader);
	close($writer);

	sub writeTranslationSequences {
		my ($contig, $sequence) = @_;
		my $sequenceLength = length($sequence);
		foreach my $frame (0 .. 2) {
			my @startIndexList = ();
			for(my $index = $frame; $index + 3 <= $sequenceLength; $index += 3) {
				my $codon = substr($sequence, $index, 3);
				if($startCodonHash{$codon}) {
					push(@startIndexList, $index);
				} elsif(@startIndexList && $stopCodonHash{$codon}) {
					writeTranslationSequence($contig, $sequence, $sequenceLength, '+', $index + 3, @startIndexList);
					@startIndexList = ();
				}
			}
		}
		(my $reverseComplementarySequence = reverse($sequence)) =~ tr/ACGT/TGCA/;
		foreach my $frame (0 .. 2) {
			my @startIndexList = ();
			for(my $index = $frame; $index + 3 <= $sequenceLength; $index += 3) {
				my $codon = substr($reverseComplementarySequence, $index, 3);
				if($startCodonHash{$codon}) {
					push(@startIndexList, $index);
				} elsif(@startIndexList && $stopCodonHash{$codon}) {
					writeTranslationSequence($contig, $reverseComplementarySequence, $sequenceLength, '-', $index + 3, @startIndexList);
					@startIndexList = ();
				}
			}
		}
	}

	sub writeTranslationSequence {
		my ($contig, $sequence, $sequenceLength, $strand, $endIndex, @startIndexList) = @_;
		if((my $length = $endIndex - $startIndexList[0]) >= $minimumLength) {
			my ($start, $end) = ('', '');
			($start, $end) = ($startIndexList[0] + 1, $endIndex) if($strand eq '+');
			($start, $end) = (($sequenceLength - $endIndex) + 1, ($sequenceLength - $startIndexList[0])) if($strand eq '-');
			my @startList = map {($_ - $startIndexList[0]) / 3 + 1} @startIndexList;
			print $writer '>', join('|', $contig, $start, $end, $strand, @startList), "\n";
			print $writer translate(substr($sequence, $startIndexList[0], $length)), "\n";
		}
	}
}

{ # ORF translation mapping
	system("diamond blastp --threads $threads --db $databasePrefix.dmnd --query $outputPrefix.translation.fasta --out $outputPrefix.blast.txt --outfmt 6 --sensitive --evalue $evalue --tmpdir $temporaryDirectory 1>&2");

	my %proteinSequenceLengthHash = ();
	{
		open(my $reader, "$databasePrefix.length.txt");
		while(my $line = <$reader>) {
			chomp($line);
			my ($protein, $sequenceLength) = split(/\t/, $line);
			$proteinSequenceLengthHash{$protein} = $sequenceLength;
		}
		close($reader);
	}
	{
		open(my $reader, "$outputPrefix.blast.txt");
		open(my $writer, "> $outputPrefix.blast.filtered.txt");
		while(my $line = <$reader>) {
			chomp($line);
			my %tokenHash = ();
			@tokenHash{qw(qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore)} = split(/\t/, $line);
			$tokenHash{'coverage'} = (($tokenHash{'pident'} / 100) * $tokenHash{'length'}) / $proteinSequenceLengthHash{$tokenHash{'sseqid'}};
			print $writer join("\t", $line, $tokenHash{'coverage'}), "\n" if($tokenHash{'coverage'} >= $minimumCoverage && $tokenHash{'qstart'} <= $tokenHash{'qend'} && $tokenHash{'sstart'} <= $tokenHash{'send'});
		}
		close($reader);
		close($writer);
	}
}

{ # Abundance estimation
	my %contigMeanDepthHash = ();
	my $genomeMeanDepth = 0;
	{
		open(my $reader, "samtools view -F 2052 -q $minimumMappingQuality $outputPrefix.sorted.bam |");
		while(my $line = <$reader>) {
			chomp($line);
			my %tokenHash = ();
			@tokenHash{'qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual'} = split(/\t/, $line);
			my @positionList = getPositionList(@tokenHash{'pos', 'cigar'});
			if(@positionList = grep {defined} @positionList) {
				$contigMeanDepthHash{$tokenHash{'rname'}} += scalar(@positionList);
				$genomeMeanDepth += scalar(@positionList);
			}
		}
		close($reader);
		my $genomeSequenceLength = 0;
		foreach my $contig (keys %contigSequenceLengthHash) {
			my $contigSequenceLength = $contigSequenceLengthHash{$contig};
			if($contigMeanDepthHash{$contig}) {
				$contigMeanDepthHash{$contig} /= $contigSequenceLength;
				$genomeSequenceLength += $contigSequenceLength;
			} else {
				print STDERR "No read mapped to '$contig'\n";
			}
		}
		$genomeMeanDepth /= $genomeSequenceLength;
		print "Genome mean depth: $genomeMeanDepth\n";
	}
	my @valueColumnList = ('readCount', 'RPKM', 'meanDepth', 'meanDepth/contig', 'meanDepth/genome');
	{
		open(my $writer, "> $outputPrefix.region.abundance.txt");
		print $writer join("\t", 'contig', 'start', 'end', 'strand', 'orthology', map {$_ eq 'RPKM' ? 'RPK' : $_} @valueColumnList), "\n";
		close($writer);
	}
	my $totalReadCount = 0;
	{
		open(my $reader, "sort --field-separator='\t' -k1,1 -k11,11g -k12,12gr -k13,13gr $outputPrefix.blast.filtered.txt |");
		open(my $writer, "| sort --field-separator='\t' -k5,5 -k1,1 -k2,2n -k3,3n -k4,4r -k6 | uniq >> $outputPrefix.region.abundance.txt");
		my %topTokenHash = ();
		$topTokenHash{'qseqid'} = '';
		while(my $line = <$reader>) {
			chomp($line);
			my %tokenHash = ();
			@tokenHash{qw(qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore), 'coverage'} = split(/\t/, $line);
			%topTokenHash = %tokenHash if($tokenHash{'qseqid'} ne $topTokenHash{'qseqid'});
			if(scalar(grep {$tokenHash{$_} != $topTokenHash{$_}} 'evalue', 'bitscore', 'coverage') == 0) {
				my ($contig, $start, $end, $strand, @startList) = split(/\|/, $tokenHash{'qseqid'});
				if($tokenHash{'qstart'} > 1) {
					my $startIndex = max(grep {$_ <= $tokenHash{'qstart'}} @startList) - 1;
					$start += $startIndex * 3 if($strand eq '+');
					$end   -= $startIndex * 3 if($strand eq '-');
				}
				$tokenHash{'orthology'} = getOrthology($tokenHash{'sseqid'});
				@tokenHash{'readCount', 'baseCount'} = getReadBaseCount("$outputPrefix.sorted.bam", $contig, $start, $end, $strand);
				my $length = $end - $start + 1;
				$tokenHash{'RPKM'} = $tokenHash{'readCount'} / ($length / 1000);
				if(($tokenHash{'meanDepth'} = $tokenHash{'baseCount'} / $length) == 0) {
					$tokenHash{'meanDepth/contig'} = 0;
					$tokenHash{'meanDepth/genome'} = 0;
				} else {
					$tokenHash{'meanDepth/contig'} = $tokenHash{'meanDepth'} / $contigMeanDepthHash{$contig};
					$tokenHash{'meanDepth/genome'} = $tokenHash{'meanDepth'} / $genomeMeanDepth;
				}
				print $writer join("\t", $contig, $start, $end, $strand, @tokenHash{'orthology', @valueColumnList}), "\n";
				$totalReadCount += $tokenHash{'readCount'};
			}
		}
		close($reader);
		close($writer);
	}
	{
		open(my $reader, "$outputPrefix.region.abundance.txt");
		open(my $writer, "> $outputPrefix.abundance.txt");
		if($orthologyDefinitionFile ne '') {
			print $writer join("\t", 'orthology', 'definition', @valueColumnList), "\n";
		} else {
			print $writer join("\t", 'orthology', @valueColumnList), "\n";
		}
		my %orthologyTokenHash = ();
		$orthologyTokenHash{'orthology'} = '';
		chomp(my $line = <$reader>);
		my @columnList = split(/\t/, $line);
		while(my $line = <$reader>) {
			chomp($line);
			my %tokenHash = ();
			@tokenHash{@columnList} = split(/\t/, $line);
			$tokenHash{'RPKM'} = $tokenHash{'RPK'} / ($totalReadCount / 1000000);
			if($tokenHash{'orthology'} ne $orthologyTokenHash{'orthology'}) {
				writeOrthologyLine(%orthologyTokenHash) if($orthologyTokenHash{'orthology'} ne '');
				%orthologyTokenHash = ();
				$orthologyTokenHash{'orthology'} = $tokenHash{'orthology'};
			}
			$orthologyTokenHash{$_} += $tokenHash{$_} foreach(@valueColumnList);
		}
		writeOrthologyLine(%orthologyTokenHash) if($orthologyTokenHash{'orthology'} ne '');
		close($reader);
		close($writer);

		sub writeOrthologyLine {
			my (%orthologyTokenHash) = @_;
			if($orthologyDefinitionFile ne '') {
				$orthologyTokenHash{'definition'} = $orthologyDefinitionHash{$orthologyTokenHash{'orthology'}};
				$orthologyTokenHash{'definition'} = '' unless(defined($orthologyTokenHash{'definition'}));
				print $writer join("\t", @orthologyTokenHash{'orthology', 'definition', @valueColumnList}), "\n";
			} else {
				print $writer join("\t", @orthologyTokenHash{'orthology', @valueColumnList}), "\n";
			}
		}
	}
}

sub loadDefaultDatabase {
	open(my $reader, "$fmapPath/FMAP_data/database");
	chomp(my $line = <$reader>);
	close($reader);
	return $1 if($line =~ /^([A-Za-z0-9_.]+)/);
}

sub getCommandPath {
	my ($command) = @_;
	chomp($command = `which $command`) if($command ne '' && $command !~ /\//);
	$command =~ s/^~\//$ENV{'HOME'}\//;
	$command = abs_path($command) if($command ne '');
	return $command;
}

sub getOrthology {
	my ($protein) = @_;
	return $_ if(defined($_ = $proteinOrthologyHash{$protein}));
	return $1 if($protein =~ /^(K[0-9]{5})_/);
	return $protein;
}

sub getReadBaseCount {
	my ($bamFile, $contig, $start, $end, $strand) = @_;
	my %readCountHash = ();
	my $baseCount = 0;
	open(my $reader, "samtools view -F 2052 -q $minimumMappingQuality $bamFile $contig:$start-$end |");
	while(my $line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		@tokenHash{'qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual'} = split(/\t/, $line);
		next if($stranded ne '' && getReadStrand($tokenHash{'flag'}) ne $strand);
		my @positionList = getPositionList(@tokenHash{'pos', 'cigar'});
		if(@positionList = grep {$start <= $_ && $_ <= $end} grep {defined} @positionList) {
			$readCountHash{$tokenHash{'qname'}} += 1;
			$baseCount += scalar(@positionList);
		}
	}
	close($reader);
	return (scalar(keys %readCountHash), $baseCount);
}

sub getReadStrand {
	my ($flag) = @_;
	if($stranded eq 'f') {
		if($flag & 1) {
			return '+' if(grep {$_ == $flag} (99, 147));
			return '-' if(grep {$_ == $flag} (83, 163));
		} else {
			return '+' if($flag == 0);
			return '-' if($flag == 16);
		}
	}
	if($stranded eq 'r') {
		if($flag & 1) {
			return '+' if(grep {$_ == $flag} (83, 163));
			return '-' if(grep {$_ == $flag} (99, 147));
		} else {
			return '+' if($flag == 16);
			return '-' if($flag == 0);
		}
	}
	return '';
}

sub getPositionList {
	my ($position, $cigar) = @_;
	my @positionList = ();
	my $index = 0;
	while($cigar =~ s/^([0-9]+)([MIDNSHP=X])//) {
		my ($length, $operation) = ($1, $2);
		if($operation eq 'M') {
			@positionList[$index .. $index + $length - 1] = $position .. $position + $length - 1;
			$index += $length;
			$position += $length;
		}
		if($operation eq 'I') {
			$index += $length;
		}
		if($operation eq 'D') {
			$position += $length;
		}
		if($operation eq 'N') {
			$position += $length;
		}
		if($operation eq 'S') {
			$index += $length;
		}
	}
	return @positionList;
}
