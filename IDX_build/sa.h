/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */

#ifndef __build_INDEX__sa__
#define __build_INDEX__sa__

#include <iostream>
#include "fasta.h"
#include "coding.h"
#include "data.h"




void SA_IS(unsigned char *s, unsigned int *SA, unsigned int n, unsigned int K, unsigned int cs, unsigned int level);

class SA_index{
private:
	///unsigned int n;
    unsigned int n_sparse;
    unsigned int suffix_count;
    //unsigned char *s_ch;
    unsigned int *SA;
    unsigned int *LUT;
    
    int encodeSuffix(char * kmer);

    
public:
    SA_index();
    ~SA_index();
    SA_index(unsigned int sparse);
	int make_sa(CFastaFile * fasta);
    
    unsigned int get_count();
    void  write_sa(FILE * file);
    
    void  write_sa(FILE * file, int32_t sparse);
    
    void  write_sa_lut(FILE * file);
    
    
    /******/
    void readData(FILE *in);
    void calcLUT(CFastaFile * fasta); //assuming sa is read from input file
};



#endif /* defined(__build_INDEX__sa__) */
