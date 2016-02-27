/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */

#ifndef __search_INDEX__fastq__
#define __search_INDEX__fastq__

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "defs.h"
#include <string>

using namespace std;


typedef enum {none, fasta, fastq, simlple_list} file_mode;

class SeqFile
{
    ifstream file;
    file_mode mode;
    
    char * seq;
    size_t maxSeqLen;
    
    string currName;
    string currName2;
    string currSeq;
    string currQual;
    
public:
    SeqFile() {
        mode = none;
        maxSeqLen = 300;
        seq = (char *) malloc (maxSeqLen * sizeof(char));
        
	}
    
    ~SeqFile() {
        
    }
    
    bool Open(std::string file_name);
    bool Close();
    
    const char * getNextSeq();
    char * getCurrName();
};


#endif /* defined(__search_INDEX__fastq__) */