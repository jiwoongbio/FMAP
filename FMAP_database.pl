# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die "ERROR in $0: ", $_[0] };

use Cwd 'abs_path';
use Getopt::Long;
use Bio::DB::Taxonomy;

(my $fmapPath = abs_path($0)) =~ s/\/[^\/]*$//;
my $dataPath = "$fmapPath/FMAP_data";
system("mkdir -p $dataPath");

GetOptions('h' => \(my $help = ''),
	's' => \(my $switchDatabase = ''),
	'r' => \(my $redownload = ''));
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl FMAP_database.pl [options] 50|90|100 [NCBI_TaxID [...]]

Options: -h       display this help message
         -s       switch database
         -r       redownload data

EOF
}
my ($unirefIdentity, @taxonIdList) = @ARGV;
die "ERROR in $0: UniRef identity must be 50|90|100.\n" unless(grep {$unirefIdentity eq $_} split(/\|/, '50|90|100'));
die "ERROR in $0: Taxon ID must be integer.\n" if(grep {$_ !~ /^[0-9]+$/} @taxonIdList);
if(@taxonIdList) {
	@taxonIdList = sort {$a <=> $b} @taxonIdList;
	@taxonIdList = @taxonIdList[0, grep {$taxonIdList[$_ - 1] != $taxonIdList[$_]} 1 .. $#taxonIdList];
}
my $prefix = join('_', "uniref$unirefIdentity", @taxonIdList);
if($switchDatabase) {
	chomp(my @databaseList = `ls $dataPath/orthology_$prefix.*.fasta`);
	if(@databaseList) {
		s/^.*\/// foreach(@databaseList);
		s/\.fasta$// foreach(@databaseList);
		@databaseList = sort @databaseList;
		switchDatabase($databaseList[-1]);
		exit 0;
	}
}
$prefix = join('.', $prefix, getTimeString());

my %taxonIdHash = ();
if(@taxonIdList) {
	if(not -r "$dataPath/nodes.dmp" or not -r "$dataPath/names.dmp" or $redownload) {
		my $URL = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz";
		my $file = "$dataPath/taxdump.tar.gz";
		system("wget --no-verbose -O $file $URL") if(not -r $file or $redownload);
		die "ERROR in $0: '$file' has zero size.\n" if(-z $file);
		system("cd $dataPath; tar -zxf taxdump.tar.gz nodes.dmp");
		system("cd $dataPath; tar -zxf taxdump.tar.gz names.dmp");
		system("rm -f $dataPath/$_") foreach('nodes', 'parents', 'names2id', 'id2names');
	}
	my $db = Bio::DB::Taxonomy->new(-source => 'flatfile', -directory => $dataPath, -nodesfile => "$dataPath/nodes.dmp", -namesfile => "$dataPath/names.dmp");
	my @taxonList = map {$db->get_taxon(-taxonid => $_)} map {eval($_)} @taxonIdList;
	setTaxonIdHash($_) foreach(@taxonList);

	sub setTaxonIdHash {
		my ($taxon) = @_;
		return if(defined($taxonIdHash{$taxon->id}));
		$taxonIdHash{$taxon->id} = $taxon;
		my @taxonList = $db->each_Descendent($taxon);
		setTaxonIdHash($_) foreach(@taxonList);
	}
}
{
	my $URL = "http://rest.genome.jp/link/uniprot/ko";
	my $file = "$dataPath/KEGG_orthology2uniprot.txt";
	system("wget --no-verbose -O $file $URL") if(not -r $file or $redownload);
	die "ERROR in $0: '$file' has zero size.\n" if(-z $file);
}
#my %uniprotOrthologyListHash = ();
#{
#	open(my $reader, "$dataPath/KEGG_orthology2uniprot.txt");
#	while(my $line = <$reader>) {
#		chomp($line);
#		my ($orthology, $uniprot) = split(/\t/, $line);
#		$orthology =~ s/^ko://;
#		$uniprot =~ s/^up://;
#		push(@{$uniprotOrthologyListHash{$uniprot}}, $orthology);
#	}
#	close($reader);
#}
my %unirefOrthologyCountHash = ();
{
	my $URL = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz";
	my $file = "$dataPath/idmapping.dat.gz";
	system("wget --no-verbose -O $file $URL") if(not -r $file or $redownload);
	die "ERROR in $0: '$file' has zero size.\n" if(-z $file);
	my %tokenHash = ('UniProtKB-AC' => '');
	open(my $reader, "gzip -dc $file | sort --field-separator='\t' -k1,1 |");
	open(my $writer, "| sort --field-separator='\t' -k1,1 -k2,2 -k3,3 | uniq | cut -f2,3 > $dataPath/gene_$prefix.txt");
	while(my $line = <$reader>) {
		chomp($line);
		my @tokenList = split(/\t/, $line);
		if($tokenList[0] ne $tokenHash{'UniProtKB-AC'}) {
			printGeneUniref() if($tokenHash{'UniProtKB-AC'} ne '');
			%tokenHash = ('UniProtKB-AC' => $tokenList[0]);
		}
		push(@{$tokenHash{$tokenList[1]}}, $tokenList[2]);
	}
	printGeneUniref() if($tokenHash{'UniProtKB-AC'} ne '');
	close($reader);
	close($writer);

	sub printGeneUniref {
		my ($taxonIdList, $geneList, $unirefList) = @tokenHash{'NCBI_TaxID', 'KEGG', "UniRef$unirefIdentity"};
		if(scalar(keys %taxonIdHash) == 0 || (defined($taxonIdList) && grep {defined($taxonIdHash{$_})} @$taxonIdList)) {
			if(defined($unirefList) && defined($geneList)) {
				print $writer join("\t", $_->[1], $_->[0], join(',', @$unirefList)), "\n" foreach(map {[$_, split(/:/, $_)]} @$geneList);
			}
#			if(defined($unirefList) && defined(my $orthologyList = $uniprotOrthologyListHash{$tokenHash{'UniProtKB-AC'}})) {
#				if(scalar(@$unirefList) == 1 && scalar(@$orthologyList) == 1) {
#					$unirefOrthologyCountHash{$unirefList->[0]}->{$orthologyList->[0]} += 1;
#				}
#			}
		}
	}
}
{
	my $orgCode = '';
	my %geneUnirefHash = ();
	open(my $reader, "$dataPath/gene_$prefix.txt");
	open(my $writerFail, "> $dataPath/gene_$prefix.no_orthology.txt");
	while(my $line = <$reader>) {
		chomp($line);
		my ($gene, $uniref) = split(/\t/, $line);
		if($gene !~ /^$orgCode:/) {
			addUnirefOrthology() if($orgCode ne '');
			($orgCode = $gene) =~ s/:.*$//;
			%geneUnirefHash = ();
		}
		$geneUnirefHash{$gene}->{$uniref} = 1;
	}
	addUnirefOrthology() if($orgCode ne '');
	close($reader);
	close($writerFail);

	sub addUnirefOrthology {
		my %geneOrthologyHash = ();
		open(my $reader, "wget --no-verbose -O - http://rest.kegg.jp/link/ko/$orgCode |");
		while(my $line = <$reader>) {
			chomp($line);
			my ($gene, $orthology) = split(/\t/, $line);
			$orthology =~ s/^ko://;
			$geneOrthologyHash{$gene}->{$orthology} = 1;
		}
		close($reader);

		foreach my $gene (keys %geneUnirefHash) {
			my @unirefList = sort keys %{$geneUnirefHash{$gene}};
			my @orthologyList = ();
			@orthologyList = sort keys %$_ if(defined($_ = $geneOrthologyHash{$gene}));
			if(scalar(@unirefList) == 1 && scalar(@orthologyList) == 1) {
				$unirefOrthologyCountHash{$unirefList[0]}->{$orthologyList[0]} += 1;
			} else {
				print $writerFail join("\t", $gene, join(',', @unirefList), join(',', @orthologyList)), "\n";
			}
		}
	}
}
{
	my $URL = "ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref$unirefIdentity/uniref$unirefIdentity.fasta.gz";
	my $file = "$dataPath/uniref$unirefIdentity.fasta.gz";
	system("wget --no-verbose -O $file $URL") if(not -r $file or $redownload);
	die "ERROR in $0: '$file' has zero size.\n" if(-z $file);
	my $printSequence = 0;
	open(my $reader, "gzip -dc $file |");
	open(my $writer, "> $dataPath/orthology_$prefix.fasta");
	open(my $writerFail, "> $dataPath/orthology_$prefix.ambiguous.txt");
	while(my $line = <$reader>) {
		chomp($line);
		if($line =~ s/^>//) {
			$printSequence = 0;
			my ($uniref) = split(/ /, $line);
			if(defined(my $orthologyCountHash = $unirefOrthologyCountHash{$uniref})) {
				if(scalar(my @orthologyList = sort keys %$orthologyCountHash) == 1) {
					print $writer ">$orthologyList[0]_$uniref\n";
					$printSequence = 1;
				} else {
					print $writerFail join("\t", $uniref, $_, $orthologyCountHash->{$_}), "\n" foreach(@orthologyList);
				}
			}
		} elsif($printSequence) {
			print $writer "$line\n";
		}
	}
	close($reader);
	close($writer);
	close($writerFail);
}
switchDatabase("orthology_$prefix");

sub getTimeString {
	my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time);
	return sprintf('%04d%02d%02d%02d%02d%02d', $year + 1900, $mon + 1, $mday, $hour, $min, $sec);
}

sub switchDatabase {
	my ($database) = @_;
	open(my $writer, "> $dataPath/database");
	print $writer "$database\n";
	close($writer);
	my $sequenceCount = 0;
	my %orthologyHash = ();
	open(my $reader, "$dataPath/$database.fasta");
	while(my $line = <$reader>) {
		chomp($line);
		if($line =~ s/^>//) {
			$sequenceCount += 1;
			(my $orthology = $line) =~ s/_UniRef.*$//;
			$orthologyHash{$orthology} = 1;
		}
	}
	close($reader);
	my $orthologyCount = scalar(keys %orthologyHash);
	print "The default database is switched to $database containing $sequenceCount UniRef sequences covering $orthologyCount KEGG orthologies.\n";
}
