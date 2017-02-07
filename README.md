# FMAP

Functional Mapping and Analysis Pipeline for metagenomics and metatranscriptomics studies

Some example results are available at the homepage: https://qbrc.swmed.edu/FMAP/.


## Features

* FMAP provides a more sensible reference protein sequence database based on [UniRef](http://www.uniprot.org/help/uniref).

* Identification of differentially-abundant genes [KEGG Orthology](http://www.genome.jp/kegg/ko.html)

* Mapping differentially-abundant genes to pathways ([KEGG Pathway](http://www.genome.jp/kegg/pathway.html))

* Mapping differentially-abundant genes to operons ([ODB (v3)](http://operondb.jp))


## Requirements

* [Perl](https://www.perl.org) - scripting language

* [R](http://www.r-project.org) - statistical computing

* [Statistics::R](http://search.cpan.org/~fangly/Statistics-R-0.33/lib/Statistics/R.pm) - Perl interface with the R statistical program
  * Use CPAN to install the module  
   ```
   perl -MCPAN -e 'install Statistics::R'
   ```
  * or download the source and compile manually  
   ```
   wget 'http://search.cpan.org/CPAN/authors/id/F/FA/FANGLY/Statistics-R-0.33.tar.gz'
   tar zxf Statistics-R-0.33.tar.gz
   cd Statistics-R-0.33
   perl Makefile.PL
   make
   make test
   make install
   ```

* Mapping program providing BLASTX search of sequencing reads: [DIAMOND](http://ab.inf.uni-tuebingen.de/software/diamond/) or [USEARCH](http://www.drive5.com/usearch/)

* Linux commands: [```wget```](https://www.gnu.org/software/wget/), ```cat```, ```sort```

* [Bio::DB::Taxonomy](http://search.cpan.org/dist/BioPerl/Bio/DB/Taxonomy.pm) - Access to a taxonomy database (which is required only if you want to build a custom database.)


## Command usages

* **FMAP_database.pl**
  * Input: [UniRef](http://www.uniprot.org/help/uniref) sequence identity (50, 90, or 100), [NCBI taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy) IDs (numerical)
```
Usage:   perl FMAP_database.pl [options] 50|90|100 [NCBI_TaxID [...]]

Options: -h       display this help message
         -r       redownload data
```

* **FMAP_prepare.pl**
```
Usage:   perl FMAP_prepare.pl [options]

Options: -h       display this help message
         -r       redownload data
         -m FILE  executable file path of mapping program, "diamond" or "usearch" [diamond]
         -k       download prebuilt KEGG files
```

* **FMAP_download.pl**
```
Usage:   perl FMAP_download.pl [options]

Options: -h       display this help message
         -m FILE  executable file path of mapping program, "diamond" or "usearch" [diamond]
         -k       download prebuilt KEGG files
         -x       download only KEGG files
```

* **FMAP_mapping.pl**
  * Input: whole metagenomic (or metatranscriptomic) shotgun sequencing reads in FASTQ or FASTA format
  * Output: best-match hits in NCBI BLAST ‑m8 (= NCBI BLAST+ ‑outfmt 6) format
```
Usage:   perl FMAP_mapping.pl [options] input1.fastq|input1.fasta [input2.fastq|input2.fasta [...]] > blastx_hits.txt

Options: -h       display this help message
         -m FILE  executable file path of mapping program, "diamond" or "usearch" [diamond]
         -p INT   number of threads [1]
         -e FLOAT maximum e-value to report alignments for "diamond" [0.001]
         -i FLOAT minimum identity for "usearch_global" [0.8]
         -t DIR   directory for temporary files [$TMPDIR or /tmp]
```

* **FMAP_quantification.pl**
  * Input: output of "FMAP_mapping.pl"
  * Output: abundances (RPKM) of KEGG orthologies
  * Output columns: [KEGG Orthology](http://www.genome.jp/kegg/ko.html) ID, orthology definition, abundance (RPKM)
```
Usage:   perl FMAP_quantification.pl [options] blast_hits1.txt [blast_hits2.txt [...]] > abundance.txt

Options: -h       display this help message
         -c       use CPM values instead of RPKM values
         -w FILE  tab-delimited text file with the first column having read names and the second column having the weights
```

* **FMAP_table.pl**
  * Input: outputs of "FMAP_quantification.pl"
  * Output: abundance table
  * Output columns: [KEGG Orthology](http://www.genome.jp/kegg/ko.html) ID, orthology definition, abundance of sample1, abundance of sample2, ...
```
Usage:   perl FMAP_table.pl [options] [name1=]abundance1.txt [[name2=]abundance2.txt [...]] > abundance_table.txt

Options: -h       display this help message
         -c       use raw counts instead of RPKM values
```

* **FMAP_comparison.pl**
  * Input: output of "FMAP_table.pl", sample group information
  * Output: comparison test statistics for orthologies
  * Output columns: [KEGG Orthology](http://www.genome.jp/kegg/ko.html) ID, orthology definition, log2 fold change, p-value, FDR-adjusted p-value, filter (pass or fail)
```
Usage:   perl FMAP_comparison.pl [options] abundance_table.txt control1[,control2[...]] case1[,case2[...]] [...] > orthology_test_stat.txt

Options: -h       display this help message
         -t STR   statistical test for comparing sample groups, "kruskal", "anova", "poisson", "quasipoisson", "metagenomeSeq" [kruskal]
         -f FLOAT fold change cutoff [2]
         -p FLOAT p-value cutoff [0.05]
         -a FLOAT FDR-adjusted p-value cutoff [1]
```

* **FMAP_pathway.pl**
  * Input: output of "FMAP_comparison.pl"
  * Output: pathways enriched in filter-passed orthologies
  * Output columns: [KEGG Pathway](http://www.genome.jp/kegg/pathway.html) ID, pathway definition, orthology count, coverage, p-value, [KEGG Orthology](http://www.genome.jp/kegg/ko.html) IDs with colors
  * [KEGG Orthology](http://www.genome.jp/kegg/ko.html) IDs with colors: input of [KEGG Pathway](http://www.genome.jp/kegg/pathway.html) mapping (http://www.kegg.jp/kegg/tool/map_pathway2.html)
```
Usage:   perl FMAP_pathway.pl [options] orthology_test_stat.txt > pathway.txt

Options: -h       display this help message
```

* **FMAP_operon.pl**
  * Input: output of "FMAP_comparison.pl"
  * Output: operons consisting of filter-passed orthologies
  * Output columns: [ODB (v3)](http://operondb.jp) known operon IDs, operon definition, log2 fold change, [KEGG Orthology](http://www.genome.jp/kegg/ko.html) IDs, [KEGG Pathway](http://www.genome.jp/kegg/pathway.html) IDs
```
Usage:   perl FMAP_operon.pl [options] orthology_test_stat.txt > operon.txt

Options: -h       display this help message
         -a       print single-gene operons as well
```

* **FMAP_all.pl**
  * Input: configuration table file
  * Input columns: group (control, ...), sample name, input file of "FMAP_mapping.pl"
  * Output: script file including all FMAP commands, all FMAP outputs
```
Usage:   perl FMAP_all.pl [options] input.config [output_prefix]

Options: -h       display this help message
         -s       generate a script, but not execute it
         -m FILE  mapping: executable file path of mapping program, "diamond" or "usearch" [diamond]
         -t INT   mapping: number of threads [1]
         -c STR   comparison: statistical test for comparing sample groups, "kruskal", "anova", "poisson", "quasipoisson", "metagenomeSeq" [kruskal]
         -f FLOAT comparison: fold change cutoff [2]
         -p FLOAT comparison: p-value cutoff [0.05]
         -a FLOAT comparison: FDR-adjusted p-value cutoff [1]
```


## Command orders

* Use the prebuilt database (UniRef90 and bacteria/archaea/fungi)
  1. FMAP_download.pl
  2. FMAP_mapping.pl
  3. FMAP_quantification.pl
  4. FMAP_table.pl
  5. FMAP_comparison.pl
  6. FMAP_pathway.pl
  7. FMAP_operon.pl

* Use a custom database (you can define UniRef and taxonomy.)
  1. FMAP_database.pl
  2. FMAP_prepare.pl
  3. FMAP_mapping.pl
  4. FMAP_quantification.pl
  5. FMAP_table.pl
  6. FMAP_comparison.pl
  7. FMAP_pathway.pl
  8. FMAP_operon.pl


## Citation

Kim J, Kim MS, Koh AY, Xie Y, Zhan X.
"FMAP: Functional Mapping and Analysis Pipeline for metagenomics and metatranscriptomics studies"
BMC Bioinformatics. 2016 Oct 10;17(1):420.
PMID: [27724866](https://www.ncbi.nlm.nih.gov/pubmed/27724866)
