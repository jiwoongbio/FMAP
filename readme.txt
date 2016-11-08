FMAP: Functional Mapping Analysis Pipeline (https://qbrc.swmed.edu/FMAP/)


REQUIREMENTS
============

1. Perl: https://www.perl.org

2. R: http://www.r-project.org

3. Perl module "Statistics::R": http://search.cpan.org/~fangly/Statistics-R-0.33/lib/Statistics/R.pm

- use CPAN to install the module

perl -MCPAN -e 'install Statistics::R'

- or download the source and compile manually

wget http://search.cpan.org/CPAN/authors/id/F/FA/FANGLY/Statistics-R-0.33.tar.gz
tar zxf Statistics-R-0.33.tar.gz
cd Statistics-R-0.33
perl Makefile.PL
make
make test
make install

- (for non-root account) Use cpanm
This method is adopted from: http://stackoverflow.com/questions/2980297/how-can-i-use-cpan-as-a-non-root-user

wget -O - http://cpanmin.us | perl - -l ~/perl5 App::cpanminus local::lib
eval `perl -I ~/perl5/lib/perl5 -Mlocal::lib`
cpanm Statistics::R

4. Mapping program providing BLASTX search of sequencing reads: "DIAMOND" or "USEARCH"

DIAMOND: http://ab.inf.uni-tuebingen.de/software/diamond/
USEARCH: http://www.drive5.com/usearch/


INPUTS AND OUTPUTS
==================

1. FMAP_mapping.pl

Input: whole metagenomic (or metatranscriptomic) shotgun sequencing reads in FASTQ or FASTA format
Output: best-match hits in NCBI BLAST ‑m8 (= NCBI BLAST+ ‑outfmt 6) format

2. FMAP_quantification.pl

Input: output of "FMAP_mapping.pl"
Output: abundances (RPKM) of KEGG orthologies
        - columns: KEGG orthology ID, orthology definition, abundance (RPKM)

3. FMAP_table.pl

Input: outputs of "FMAP_quantification.pl"
Output: abundance table
        - columns: KEGG orthology ID, orthology definition, abundance of sample1, abundance of sample2, ...

4. FMAP_comparison.pl

Input: output of "FMAP_table.pl", sample group information
Output: comparison test statistics for orthologies
        - columns: KEGG orthology ID, orthology definition, log2 fold change, p-value, FDR-adjusted p-value, filter (pass or fail)

5. FMAP_operon.pl

Input: output of "FMAP_comparison.pl"
Output: operons consisting of filter-passed orthologies
        - columns: ODB (v3) known operon IDs, operon definition, log2 fold change, KEGG orthology IDs, KEGG pathways

6. FMAP_pathway.pl

Input: output of "FMAP_comparison.pl"
Output: pathways enriched in filter-passed orthologies
        - columns: KEGG pathway ID, pathway definition, orthology count, coverage, p-value, KEGG orthology IDs with colors
        - KEGG orthology IDs with colors: input of KEGG pathway mapping (http://www.kegg.jp/kegg/tool/map_pathway2.html)

7. FMAP_all.pl

Input: configuration table file
       - columns: group (control, ...), sample name, input file of "FMAP_mapping.pl"
Output: script file including all FMAP commands, all FMAP outputs


USAGES
======

0. FMAP_download.pl

Usage:   perl FMAP_download.pl [options]

Options: -h       display this help message
         -m FILE  executable file path of mapping program, "diamond" or "usearch"

1. FMAP_mapping.pl

Usage:   perl FMAP_mapping.pl [options] input1.fastq|input1.fasta [input2.fastq|input2.fasta [...]] > blastx_hits.txt

Options: -h       display this help message
         -m FILE  executable file path of mapping program, "diamond" or "usearch" [diamond]
         -p INT   number of threads [1]
         -e FLOAT maximum e-value to report alignments for "diamond" [0.001]
         -i FLOAT minimum identity for "usearch_global" [0.8]
         -t DIR   directory for temporary files [\$TMPDIR or /tmp]

2. FMAP_quantification.pl

Usage:   perl FMAP_quantification.pl [options] blast_hits1.txt [blast_hits2.txt [...]] > abundance.txt

Options: -h       display this help message
         -c       use CPM values instead of RPKM values
         -w FILE  tab-delimited text file with the first column having read names and the second column having the weights

3. FMAP_table.pl

Usage:   perl FMAP_table.pl [options] [name1=]abundance1.txt [[name2=]abundance2.txt [...]] > abundance_table.txt

Options: -h       display this help message
         -c       use raw counts instead of RPKM values

4. FMAP_comparison.pl

Usage:   perl FMAP_comparison.pl [options] abundance_table.txt control1[,control2[...]] case1[,case2[...]] [...] > orthology_test_stat.txt

Options: -h       display this help message
         -c STR   statistical test for comparing sample groups, "kruskal" or "anova", "poisson", "quasipoisson" [kruskal]
         -f FLOAT fold change cutoff [2]
         -p FLOAT p-value cutoff [0.05]
         -a FLOAT FDR-adjusted p-value cutoff [1]

5. FMAP_operon.pl

Usage:   perl FMAP_operon.pl [options] orthology_test_stat.txt > operon.txt

Options: -h       display this help message
         -a       print single-gene operons as well

6. FMAP_pathway.pl

Usage:   perl FMAP_pathway.pl [options] orthology_test_stat.txt > pathway.txt

Options: -h       display this help message

7. FMAP_all.pl

Usage:   perl FMAP_all.pl [options] input.config [output_prefix]

Options: -h       display this help message
         -s       generate a script, but not execute it
         -m FILE  mapping: executable file path of mapping program, "diamond" or "usearch" [diamond]
         -t INT   mapping: number of threads [1]
         -c STR   comparison: statistical test for comparing sample groups, "kruskal" or "anova", "poisson", "quasipoisson" [kruskal]
         -f FLOAT comparison: fold change cutoff [2]
         -p FLOAT comparison: p-value cutoff [0.05]
         -a FLOAT comparison: FDR-adjusted p-value cutoff [1]


CITATION
========

Kim J, Kim MS, Koh AY, Xie Y, Zhan X.
"FMAP: Functional Mapping and Analysis Pipeline for metagenomics and metatranscriptomics studies"
BMC Bioinformatics. 2016 Oct 10;17(1):420.
PMID: 27724866
