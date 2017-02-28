# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
use Cwd 'abs_path';
use Getopt::Long;

(my $fmapPath = abs_path($0)) =~ s/\/[^\/]*$//;

GetOptions('h' => \(my $help = ''), 'a' => \(my $all = ''));
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl FMAP_assembly_operon.pl [options] FMAP_assembly.region.txt > FMAP_assembly_operon.txt

Options: -h       display this help message
         -a       print single-gene operons as well

EOF
}

my ($regionFile) = @ARGV;
foreach($regionFile, "$fmapPath/FMAP_data/known_operon.KEGG_orthology.txt") {
	die "ERROR: '$_' is not readable.\n" unless(-r $_);
}
my %orthologyContigStrandIndexStartEndListHash = ();
{
	open(my $reader, "sort --field-separator='\t' -k1,1 -k4,4 -k2,2n -k3,3n $regionFile |");
	my ($contigStrand, @orthologyStartEndList) = ('');
	chomp(my $line = <$reader>);
	my @columnList = split(/\t/, $line);
	while(my $line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		@tokenHash{@columnList} = split(/\t/, $line);
		if(join("\t", @tokenHash{'contig', 'strand'}) ne $contigStrand) {
			addOrthology($contigStrand, @orthologyStartEndList) if($contigStrand ne '');
			($contigStrand, @orthologyStartEndList) = (join("\t", @tokenHash{'contig', 'strand'}));
		}
		push(@orthologyStartEndList, [@tokenHash{'orthology', 'start', 'end'}]);
	}
	addOrthology($contigStrand, @orthologyStartEndList) if($contigStrand ne '');
	close($reader);

	sub addOrthology {
		my ($contigStrand, @orthologyStartEndList) = @_;
		foreach my $index (0 .. $#orthologyStartEndList) {
			my ($orthology, $start, $end) = @{$orthologyStartEndList[$index]};
			push(@{$orthologyContigStrandIndexStartEndListHash{$orthology}}, [$contigStrand, $index, $start, $end]);
		}
	}
}
{
	open(my $reader, "$fmapPath/FMAP_data/known_operon.KEGG_orthology.txt");
	open(my $writer, "| sort --field-separator='\t' -k1,1 -k2,2n -k3,3n -k4");
	print $writer join("\t", 'contig', 'start', 'end', 'strand', 'operons', 'orthologies', 'regions'), "\n";
	while(my $line = <$reader>) {
		chomp($line);
		my ($operons, $orthologies) = split(/\t/, $line);
		my %orthologyHash = ();
		$orthologyHash{$_} = 1 foreach(split(/ /, $orthologies));
		my $orthologyCount = scalar(keys %orthologyHash);
		if($orthologyCount > 1 || $all) {
			my %contigStrandIndexStartEndOrthologyListHash = ();
			foreach my $orthology (keys %orthologyHash) {
				if(defined(my $contigStrandIndexStartEndList = $orthologyContigStrandIndexStartEndListHash{$orthology})) {
					foreach(@$contigStrandIndexStartEndList) {
						my ($contigStrand, $index, $start, $end) = @$_;
						push(@{$contigStrandIndexStartEndOrthologyListHash{$contigStrand}}, [$index, $start, $end, $orthology]);
					}
				}
			}
			foreach my $contigStrand (keys %contigStrandIndexStartEndOrthologyListHash) {
				my ($contig, $strand) = split(/\t/, $contigStrand);
				my @indexStartEndOrthologyList = @{$contigStrandIndexStartEndOrthologyListHash{$contigStrand}};
				@indexStartEndOrthologyList = sort {$a->[0] <=> $b->[0]} @indexStartEndOrthologyList;
				my @indexList = map {$_->[0]} @indexStartEndOrthologyList;
				my @startIndexList = (0, (grep {$indexList[$_ - 1] + 1 < $indexList[$_]} 1 .. $#indexList), scalar(@indexList));
				foreach(map {[$startIndexList[$_] .. $startIndexList[$_ + 1] - 1]} 0 .. $#startIndexList - 1) {
					my @orthologyStartEndList = map {[$_->[3], $_->[1], $_->[2]]} @indexStartEndOrthologyList[@$_];
					my %orthologyHash = ();
					$orthologyHash{$_} = 1 foreach(map {$_->[0]} @orthologyStartEndList);
					if(scalar(keys %orthologyHash) == $orthologyCount) {
						my ($start, $end) = ($orthologyStartEndList[0]->[1], $orthologyStartEndList[-1]->[2]);
						@orthologyStartEndList = reverse(@orthologyStartEndList) if($strand eq '-');
						my $orthologies = join(' ', map {$_->[0]} @orthologyStartEndList);
						my $regions = join(' ', map {"$_->[1]-$_->[2]"} @orthologyStartEndList);
						print $writer join("\t", $contig, $start, $end, $strand, $operons, $orthologies, $regions), "\n";
					}
				}
			}
		}
	}
	close($reader);
	close($writer);
}
