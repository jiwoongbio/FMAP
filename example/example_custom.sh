# Please modify the following variables before execution
FMAP_DIR=..
MAPPING_THREADS=4

# Download example sequencing data
wget -q -O - http://qbrc.swmed.edu/FMAP/example/control1.fastq.gz > control1.fastq.gz
wget -q -O - http://qbrc.swmed.edu/FMAP/example/control2.fastq.gz > control2.fastq.gz
wget -q -O - http://qbrc.swmed.edu/FMAP/example/control3.fastq.gz > control3.fastq.gz
wget -q -O - http://qbrc.swmed.edu/FMAP/example/control4.fastq.gz > control4.fastq.gz
wget -q -O - http://qbrc.swmed.edu/FMAP/example/control5.fastq.gz > control5.fastq.gz
wget -q -O - http://qbrc.swmed.edu/FMAP/example/case1.fastq.gz > case1.fastq.gz
wget -q -O - http://qbrc.swmed.edu/FMAP/example/case2.fastq.gz > case2.fastq.gz
wget -q -O - http://qbrc.swmed.edu/FMAP/example/case3.fastq.gz > case3.fastq.gz
wget -q -O - http://qbrc.swmed.edu/FMAP/example/case4.fastq.gz > case4.fastq.gz

# FMAP_database.pl: build a custom database of "Fusobacterium nucleatum" using "UniRef100"
perl $FMAP_DIR/FMAP_database.pl 100 851

# FMAP_prepare.pl: prepare database files and download KEGG and operon data
perl $FMAP_DIR/FMAP_prepare.pl

# FMAP_mapping.pl: mapping sequencing reads to reference proteins
perl $FMAP_DIR/FMAP_mapping.pl -p $MAPPING_THREADS control1.fastq.gz > control1.mapping.txt
perl $FMAP_DIR/FMAP_mapping.pl -p $MAPPING_THREADS control2.fastq.gz > control2.mapping.txt
perl $FMAP_DIR/FMAP_mapping.pl -p $MAPPING_THREADS control3.fastq.gz > control3.mapping.txt
perl $FMAP_DIR/FMAP_mapping.pl -p $MAPPING_THREADS control4.fastq.gz > control4.mapping.txt
perl $FMAP_DIR/FMAP_mapping.pl -p $MAPPING_THREADS control5.fastq.gz > control5.mapping.txt
perl $FMAP_DIR/FMAP_mapping.pl -p $MAPPING_THREADS case1.fastq.gz > case1.mapping.txt
perl $FMAP_DIR/FMAP_mapping.pl -p $MAPPING_THREADS case2.fastq.gz > case2.mapping.txt
perl $FMAP_DIR/FMAP_mapping.pl -p $MAPPING_THREADS case3.fastq.gz > case3.mapping.txt
perl $FMAP_DIR/FMAP_mapping.pl -p $MAPPING_THREADS case4.fastq.gz > case4.mapping.txt

# FMAP_quantification.pl: calculate gene abundances
perl $FMAP_DIR/FMAP_quantification.pl control1.mapping.txt > control1.abundance.txt
perl $FMAP_DIR/FMAP_quantification.pl control2.mapping.txt > control2.abundance.txt
perl $FMAP_DIR/FMAP_quantification.pl control3.mapping.txt > control3.abundance.txt
perl $FMAP_DIR/FMAP_quantification.pl control4.mapping.txt > control4.abundance.txt
perl $FMAP_DIR/FMAP_quantification.pl control5.mapping.txt > control5.abundance.txt
perl $FMAP_DIR/FMAP_quantification.pl case1.mapping.txt > case1.abundance.txt
perl $FMAP_DIR/FMAP_quantification.pl case2.mapping.txt > case2.abundance.txt
perl $FMAP_DIR/FMAP_quantification.pl case3.mapping.txt > case3.abundance.txt
perl $FMAP_DIR/FMAP_quantification.pl case4.mapping.txt > case4.abundance.txt

# FMAP_table.pl: generate gene abundance table
perl $FMAP_DIR/FMAP_table.pl \
control1=control1.abundance.txt \
control2=control2.abundance.txt \
control3=control3.abundance.txt \
control4=control4.abundance.txt \
control5=control5.abundance.txt \
case1=case1.abundance.txt \
case2=case2.abundance.txt \
case3=case3.abundance.txt \
case4=case4.abundance.txt \
> example.table.txt

# FMAP_comparison.pl: compare sample groups and identify differentially-abundant genes
perl $FMAP_DIR/FMAP_comparison.pl example.table.txt 'control1,control2,control3,control4,control5' 'case1,case2,case3,case4' > example.comparison.txt

# FMAP_operon.pl: mapping differentially-abundant genes to operons
perl $FMAP_DIR/FMAP_operon.pl example.comparison.txt > example.operon.txt

# FMAP_pathway.pl: mapping differentially-abundant genes to pathways
perl $FMAP_DIR/FMAP_pathway.pl example.comparison.txt > example.pathway.txt
