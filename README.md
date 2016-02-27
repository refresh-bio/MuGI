##MuGI - **Mu**ltiple **G**enome **I**ndex


###MuGI—What is it?

Multiple Genome Index is a compressed index of collections of genomes of the same species. It allows to ask for both exact and approximate queries. The input data must be given as a reference sequence and Variant Call Format (VCF) file describing the differences between individual genomes and reference genome.


####The architecture of  MuGI

MuGI is composed of several programs that were used in our experiments on the genomes from the 1000 Genomes project. A description of the tool can be found in our paper (see below).


####How good is MuGI?

The collection of 1092 diploid human genomes (6.7TB of raw sequences) can be stored in 6.5–32GB of RAM. The software is configurable and different tradeoffs can be chosen depending on what is more important: fast queries or so small size that fits even a PC with 8GB of RAM. For the moderate settings the index occupies 11.4GB of RAM and for patterns of length 120–170 answers: exact queries in 43μs, 2-mismatch queries in 151μs, 5-mismatch queries in 4461μs. All matches up to given threshold in all individuals are reported. (Test machine: 3.4GHz i7 4770 CPU.)

For more details please consult our paper.

####Terms of use of MuGI

MuGI is in general a free compression program available in source code release. 

####Publications

+ Danek, A., Deorowicz, S., Grabowski, Sz., Indexes of Large Genome Collections on a PC, [PLOS ONE](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0109384), 2014; 9(10):e109384


#### Developers

The MuGI algorithm was invented by Agnieszka Danek, Sebastian Deorowicz, and Szymon Grabowski.
The implementation is by [Agnieszka Danek](https://github.com/agnieszkadanek).


###Usage

Here we describe the usage of the tools to build the index (MuGI_build), perform the search (MuGI_search), and parse the result to a readable form (MuGI_parse). We also provide a code generating example queries (MuGI_exGen). All these tools are implemented in C/C++.

## MuGI build
The program builds the index with set sparsity and length of k-mers, based on the input refer- ence sequence and Variant Call Format (VCF) file describing the differences between individual genomes and reference genome.

Usage:
```MuGI_build <input.fasta> <input.vcf> <k> <sp> <pl>```

Parameters:
+ ```input.fasta``` - name of the input FASTA file with the reference • input.vcf — name of the input VCF file
+ ```k``` — desired length of k-mers in the index
+ ```sp``` — desired sparsity of the index
+ ```pl``` — ploidity of the individuals in the collection (only 1 and 2 supported) 

The output files:
+ ```index-<k>-<sp>``` — built index 
+ ```index-dict``` — dictionary file

The following are required to build the sources:
- Intel Threading Building Blocks library (https://www.threadingbuildingblocks.org)
- Compiler with OpenMP support (http://openmp.org)

If the ploidity is set to 2, any haploid genotype call is treated as diploid. The diploid genotypes should be phased. The first word in successive chromosome headers should match successive chromosomes described in the VCF file (its first column). The positions of the variants in the VCF file should be within the current chromosome (first column). During the process they are adjusted to present the position within the whole reference based on the size of chromosomes sequences. The index-dict file is created along with the index. It contains data about sizes of all chromosomes and names of all individuals. It is used to parse the search results to a readable form.
 
For testing we implemented the possibility to build few indexes with the same k-mer length and different sparsities (from 1 up to < sp >) at the same run. To use this mode, CREATE_ALL_INDEXES must be defined in the defs.h file (see line 14 of defs.h).

##  MuGI search
The program performs the search of the input queries in the compressed collection (input index file) and reports to the stdout statistics about queries execution times. Asmlib library (http://www.agner.org/optimize/asmlib.zip) is required to build the sources.

Usage:
- ```MuGI_search <index-file> <queries-file> {<maxError>}```

Parameters:
+ ```index-file``` — path to file with the index
+ ```queries-file``` — path to file with the sequences to search (FASTA or FASTQ or plain sequences - one in each line)
+ ```maxError``` — optional, maximum number of mismatches in found sequence (default 0,
exact search)

The output file:
+ ```result.out``` — binary file with the result

## MuGI parse
The program parses the binary file with result (output of MuGI_search) and creates a text file with readable data.
For each queries, at first its name (for FASTQ and FASTA) or ordinal number (for simple list) is written. Then, all found positions are reported (chromosome number and position within it) along with list of all individuals in which the match was found. In a ”hex” mode (0), the list is represented by a bit vector written in a hexadecimal form, where each jth bit set corresponds to match found in jth haploid genome. In the ”full” mode (1) the list consist of names of all individuals with the match. For diploid individuals, to distinguish between two haploid sequences of each, a suffix ”-1” or ”-2” is added to the name. In the bit vector representation, each diploid individual is represented by two consecutive bits.

Usage:
```MuGI_parse <dict> <file_to_parse> <output_file> <mode>```

Parameters:
+ ```dict``` — dictionary file of the index
+ ```file_to_parse``` — binary file with the result
+ ```output_file``` — name of the output file
+ ```mode``` — optional, mode of the output: 0 (hexadecimal individuals representation) or 1 (full individual names) default: 0
  
## MuGI genEx
The program generates a file with a set number of queries. The file an be in FASTQ or FASTA format, or be a text file with list of sequences (three acceptable formats of queries for MuGI_search). It uses the created index file, to take an excerpts from the collection of genomes. As only Reference (REF), Variant Database (VD) and Bit Vectors (BV) are used, the k and sp parameters of the index building process does not matter. The final queries are made as described in the main paper.

Usage:
```MuGI_genEx <index-file> <out-file> <no_reads> {<mode>}```

Parameters:
+ ```index-file``` — path to file with the index
+ ```out-file``` — output file with example reads
+ ```no_reads``` — number of reads to generate
+ ```mode``` — optional, output format: "fastq" or "fasta" or "list" (default: "list")
