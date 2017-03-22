# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die "ERROR in $0: ", $_[0] };

use Getopt::Long qw(:config no_ignore_case);
use XML::LibXML;

GetOptions('h' => \(my $help = ''),
	'a' => \(my $assemblyInsteadOfGenome = ''));
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl FMAP_download_genome.pl [options] NCBI_TaxID [...] > genome.fasta

Options: -h       display this help message
         -a       assembly instead of genome

EOF
}
my (@taxonIdList) = @ARGV;
die "ERROR in $0: Taxon ID must be integer.\n" if(grep {$_ !~ /^[0-9]+$/} @taxonIdList);
if(@taxonIdList) {
	@taxonIdList = sort {$a <=> $b} @taxonIdList;
	@taxonIdList = @taxonIdList[0, grep {$taxonIdList[$_ - 1] != $taxonIdList[$_]} 1 .. $#taxonIdList];
}

my $baseURL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils';
foreach my $taxonId (@taxonIdList){
	my %paramHash = ('db' => 'nuccore', 'usehistory' => 'y', 'retmax' => 500, 'rettype' => 'fasta', 'retmode' => 'text');
	if($assemblyInsteadOfGenome) {
		$paramHash{'term'} = sprintf('txid%d[Organism:exp] AND nucleotide_assembly[Filter] AND RefSeq[Filter]', $taxonId);
	} else {
		$paramHash{'term'} = sprintf('txid%d[Organism:exp] AND nucleotide_genome[Filter] AND RefSeq[Filter]', $taxonId);
	}
	{
		my $queryString = join('&', map {"$_=$paramHash{$_}"} 'db', 'term', 'usehistory');
		my $xmlString = `wget --no-verbose -O - '$baseURL/esearch.fcgi?$queryString'`;
		if($xmlString ne '') {
			my $dom = XML::LibXML->load_xml(string => $xmlString);
			my $root = $dom->documentElement();
			unless(getChildNodeList($root, 'ERROR')) {
				($paramHash{'count'}) = map {$_->textContent} getChildNodeList($root, 'Count');
				($paramHash{'WebEnv'}) = map {$_->textContent} getChildNodeList($root, 'WebEnv');
				($paramHash{'query_key'}) = map {$_->textContent} getChildNodeList($root, 'QueryKey');
			}
		}
		sleep(1);
	}
	my $count = 0;
	if(defined($paramHash{'count'}) && defined($paramHash{'WebEnv'}) && defined($paramHash{'query_key'})) {
		for($paramHash{'retstart'} = 0; $paramHash{'retstart'} < $paramHash{'count'}; $paramHash{'retstart'} += $paramHash{'retmax'}) {
			my $queryString = join('&', map {"$_=$paramHash{$_}"} 'db', 'WebEnv', 'query_key', 'retstart', 'retmax', 'rettype', 'retmode');
			open(my $reader, "wget --no-verbose -O - '$baseURL/efetch.fcgi?$queryString' |");
			my @lineList = ();
			while(my $line = <$reader>) {
				next if($line eq "\n");
				if($line =~ s/^>/>$taxonId:/) {
					if(scalar(@lineList) > 1) {
						print @lineList;
						$count += 1;
					}
					@lineList = ();
				}
				push(@lineList, $line);
			}
			if(scalar(@lineList) > 1) {
				print @lineList;
				$count += 1;
			}
			close($reader);
			sleep(1);
		}
	} else {
		die "ERROR in $0: Download Failed\n";
	}
	print STDERR "$count genomes of taxon $taxonId are downloaded.\n";
}

sub getChildNodeList {
	my ($node, @childNodeTagNameList) = @_;
	my @childNodeList = ();
	foreach my $childNode ($node->getChildrenByTagName(shift @childNodeTagNameList)) {
		if(@childNodeTagNameList) {
			push(@childNodeList, getChildNodeList($childNode, @childNodeTagNameList));
		} else {
			push(@childNodeList, $childNode);
		}
	}
	return @childNodeList;
}
