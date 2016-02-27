/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */

#ifndef __build_INDEX__vcf__
#define __build_INDEX__vcf__


#define BUFFER_SIZE 512000
#define INITIAL_VT_LINE_WIDTH 5000
#define NO_LINES_PER_THREAD 512

#include <iostream>
#include "defs.h"
#include "variant_list.h"
#include "bit_vectors.h"
//#include <omp.h>
#include "omp.h"
#include "data.h"
#include <errno.h>


using namespace std;




class VCF{
private:
    
    FILE * file;
    uint64_t file_size;
    
    
    
    uint32_t no_var;
    uint32_t no_gen;
    string * genomesNames;
    
    
    
public:
   
    
    
    VCF();
    
	~VCF() {
		if(file)
			fclose(file);
    }
        

	bool Open(string file_name);
	bool ProcessVCFData(VariantList * vl, BitVectors * bv);
    
    bool ProcessVCFDataPloidity1(VariantList * vl, BitVectors * bv);
    bool ProcessVCFDataPloidity2(VariantList * vl, BitVectors * bv);
    
    
    
    bool Close();
    
    
    uint32_t get_no_var();
    uint32_t get_no_gen();
    
    string * get_gen_names();
};



#endif /* defined(__build_INDEX__vcf__) */
