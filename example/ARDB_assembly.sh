# Please modify the following variables before execution
FMAP_DIR=..
MAPPING_THREADS=4

if ! [ -x "$(command -v spades.py)" ]; then
	echo "'spades.py' is not executable."
	exit
fi
if ! [ -x "$(command -v centrifuge)" ]; then
	echo "'centrifuge' is not executable."
	exit
fi

# Download example sequencing data
wget --no-verbose -O P6.1.fastq.gz http://qbrc.swmed.edu/FMAP/example/P6.1.fastq.gz
wget --no-verbose -O P6.2.fastq.gz http://qbrc.swmed.edu/FMAP/example/P6.2.fastq.gz
wget --no-verbose -O P7.1.fastq.gz http://qbrc.swmed.edu/FMAP/example/P7.1.fastq.gz
wget --no-verbose -O P7.2.fastq.gz http://qbrc.swmed.edu/FMAP/example/P7.2.fastq.gz
wget --no-verbose -O P8.1.fastq.gz http://qbrc.swmed.edu/FMAP/example/P8.1.fastq.gz
wget --no-verbose -O P8.2.fastq.gz http://qbrc.swmed.edu/FMAP/example/P8.2.fastq.gz
wget --no-verbose -O P9.1.fastq.gz http://qbrc.swmed.edu/FMAP/example/P9.1.fastq.gz
wget --no-verbose -O P9.2.fastq.gz http://qbrc.swmed.edu/FMAP/example/P9.2.fastq.gz
wget --no-verbose -O P10.1.fastq.gz http://qbrc.swmed.edu/FMAP/example/P10.1.fastq.gz
wget --no-verbose -O P10.2.fastq.gz http://qbrc.swmed.edu/FMAP/example/P10.2.fastq.gz
wget --no-verbose -O P2.1.fastq.gz http://qbrc.swmed.edu/FMAP/example/P2.1.fastq.gz
wget --no-verbose -O P2.2.fastq.gz http://qbrc.swmed.edu/FMAP/example/P2.2.fastq.gz
wget --no-verbose -O P3.1.fastq.gz http://qbrc.swmed.edu/FMAP/example/P3.1.fastq.gz
wget --no-verbose -O P3.2.fastq.gz http://qbrc.swmed.edu/FMAP/example/P3.2.fastq.gz
wget --no-verbose -O P5.1.fastq.gz http://qbrc.swmed.edu/FMAP/example/P5.1.fastq.gz
wget --no-verbose -O P5.2.fastq.gz http://qbrc.swmed.edu/FMAP/example/P5.2.fastq.gz

# SPAdes: metagenome assembly
spades.py --meta -1 P6.1.fastq.gz -2 P6.2.fastq.gz -o P6.SPAdes
spades.py --meta -1 P7.1.fastq.gz -2 P7.2.fastq.gz -o P7.SPAdes
spades.py --meta -1 P8.1.fastq.gz -2 P8.2.fastq.gz -o P8.SPAdes
spades.py --meta -1 P9.1.fastq.gz -2 P9.2.fastq.gz -o P9.SPAdes
spades.py --meta -1 P10.1.fastq.gz -2 P10.2.fastq.gz -o P10.SPAdes
spades.py --meta -1 P2.1.fastq.gz -2 P2.2.fastq.gz -o P2.SPAdes
spades.py --meta -1 P3.1.fastq.gz -2 P3.2.fastq.gz -o P3.SPAdes
spades.py --meta -1 P5.1.fastq.gz -2 P5.2.fastq.gz -o P5.SPAdes

# FMAP_prepare.pl: prepare ARDB database files
perl $FMAP_DIR/FMAP_prepare.pl -a

# FMAP_assembly.pl: annotate genome assembly and calculate gene abundances
perl $FMAP_DIR/FMAP_assembly.pl -p $MAPPING_THREADS P6.ARDB_assembly P6.SPAdes/scaffolds.fasta P6.1.fastq.gz,P6.2.fastq.gz > P6.ARDB_assembly.summary.txt
perl $FMAP_DIR/FMAP_assembly.pl -p $MAPPING_THREADS P7.ARDB_assembly P7.SPAdes/scaffolds.fasta P7.1.fastq.gz,P7.2.fastq.gz > P7.ARDB_assembly.summary.txt
perl $FMAP_DIR/FMAP_assembly.pl -p $MAPPING_THREADS P8.ARDB_assembly P8.SPAdes/scaffolds.fasta P8.1.fastq.gz,P8.2.fastq.gz > P8.ARDB_assembly.summary.txt
perl $FMAP_DIR/FMAP_assembly.pl -p $MAPPING_THREADS P9.ARDB_assembly P9.SPAdes/scaffolds.fasta P9.1.fastq.gz,P9.2.fastq.gz > P9.ARDB_assembly.summary.txt
perl $FMAP_DIR/FMAP_assembly.pl -p $MAPPING_THREADS P10.ARDB_assembly P10.SPAdes/scaffolds.fasta P10.1.fastq.gz,P10.2.fastq.gz > P10.ARDB_assembly.summary.txt
perl $FMAP_DIR/FMAP_assembly.pl -p $MAPPING_THREADS P2.ARDB_assembly P2.SPAdes/scaffolds.fasta P2.1.fastq.gz,P2.2.fastq.gz > P2.ARDB_assembly.summary.txt
perl $FMAP_DIR/FMAP_assembly.pl -p $MAPPING_THREADS P3.ARDB_assembly P3.SPAdes/scaffolds.fasta P3.1.fastq.gz,P3.2.fastq.gz > P3.ARDB_assembly.summary.txt
perl $FMAP_DIR/FMAP_assembly.pl -p $MAPPING_THREADS P5.ARDB_assembly P5.SPAdes/scaffolds.fasta P5.1.fastq.gz,P5.2.fastq.gz > P5.ARDB_assembly.summary.txt

# FMAP_table.pl: generate gene abundance table
perl $FMAP_DIR/FMAP_table.pl -d \
	P6=P6.ARDB_assembly.abundance.txt \
	P7=P7.ARDB_assembly.abundance.txt \
	P8=P8.ARDB_assembly.abundance.txt \
	P9=P9.ARDB_assembly.abundance.txt \
	P10=P10.ARDB_assembly.abundance.txt \
	P2=P2.ARDB_assembly.abundance.txt \
	P3=P3.ARDB_assembly.abundance.txt \
	P5=P5.ARDB_assembly.abundance.txt \
	> ARDB_assembly.table.txt

# FMAP_comparison.pl: compare sample groups and identify differentially-abundant genes
perl $FMAP_DIR/FMAP_comparison.pl ARDB_assembly.table.txt 'P10,P6,P7,P8,P9' 'P2,P3,P5' > ARDB_assembly.comparison.txt

# FMAP_assembly_centrifuge.pl: map taxons to genes
wget --no-verbose -O p_compressed.tar.gz ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/data/p_compressed.tar.gz
tar -zxf p_compressed.tar.gz

perl $FMAP_DIR/FMAP_assembly_centrifuge.pl -p $MAPPING_THREADS P6.ARDB_assembly.region.abundance.txt P6.SPAdes/scaffolds.fasta p_compressed > P6.ARDB_assembly.region.abundance.taxon.txt
perl $FMAP_DIR/FMAP_assembly_centrifuge.pl -p $MAPPING_THREADS P7.ARDB_assembly.region.abundance.txt P7.SPAdes/scaffolds.fasta p_compressed > P7.ARDB_assembly.region.abundance.taxon.txt
perl $FMAP_DIR/FMAP_assembly_centrifuge.pl -p $MAPPING_THREADS P8.ARDB_assembly.region.abundance.txt P8.SPAdes/scaffolds.fasta p_compressed > P8.ARDB_assembly.region.abundance.taxon.txt
perl $FMAP_DIR/FMAP_assembly_centrifuge.pl -p $MAPPING_THREADS P9.ARDB_assembly.region.abundance.txt P9.SPAdes/scaffolds.fasta p_compressed > P9.ARDB_assembly.region.abundance.taxon.txt
perl $FMAP_DIR/FMAP_assembly_centrifuge.pl -p $MAPPING_THREADS P10.ARDB_assembly.region.abundance.txt P10.SPAdes/scaffolds.fasta p_compressed > P10.ARDB_assembly.region.abundance.taxon.txt
perl $FMAP_DIR/FMAP_assembly_centrifuge.pl -p $MAPPING_THREADS P2.ARDB_assembly.region.abundance.txt P2.SPAdes/scaffolds.fasta p_compressed > P2.ARDB_assembly.region.abundance.taxon.txt
perl $FMAP_DIR/FMAP_assembly_centrifuge.pl -p $MAPPING_THREADS P3.ARDB_assembly.region.abundance.txt P3.SPAdes/scaffolds.fasta p_compressed > P3.ARDB_assembly.region.abundance.taxon.txt
perl $FMAP_DIR/FMAP_assembly_centrifuge.pl -p $MAPPING_THREADS P5.ARDB_assembly.region.abundance.txt P5.SPAdes/scaffolds.fasta p_compressed > P5.ARDB_assembly.region.abundance.taxon.txt

perl $FMAP_DIR/FMAP_assembly_heatmap.pl -c ARDB_assembly.comparison.txt \
	P6=P6.ARDB_assembly.region.abundance.taxon.txt \
	P7=P7.ARDB_assembly.region.abundance.taxon.txt \
	P8=P8.ARDB_assembly.region.abundance.taxon.txt \
	P9=P9.ARDB_assembly.region.abundance.taxon.txt \
	P10=P10.ARDB_assembly.region.abundance.taxon.txt \
	P2=P2.ARDB_assembly.region.abundance.taxon.txt \
	P3=P3.ARDB_assembly.region.abundance.taxon.txt \
	P5=P5.ARDB_assembly.region.abundance.taxon.txt \
	> ARDB_assembly.heatmap.html
