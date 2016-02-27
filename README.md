##MuGI - **Mu**ltiple **G**enome **I**ndex

###Getting Started 

```
# Clone repository
git clone https://github.com/refresh-bio/mugi.git
# Build 
cd tgc
cd src
make
cd ..
# Repeat experiments described in the publication
./run -abcdefgijklcd tgc
```


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
