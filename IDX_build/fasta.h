/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */

#ifndef _FASTA_H
#define _FASTA_H

#include "defs.h"
#include "coding.h"
#include <iostream>
#include <vector>
using namespace std;

typedef enum {mode_fa_none, mode_fa_read, mode_fa_write, mode_fa_read_ra} t_fa_mode;

// ********************************************************************************************
class CFastaSequence {
public:
	uint32_t size;
	string seq_name;
    string chrom;
	uint32_t line_len;
	uint32_t start_pos;
	uchar eol_type;

public:
	CFastaSequence();
	CFastaSequence(uchar *_raw_data, string &_seq_name, string &_chrom, uint32 _line_len, uint32 _start_pos, uint32 _size, uchar _eol_type);
    CFastaSequence(const CFastaSequence &x);
	~CFastaSequence();
};

// ********************************************************************************************
class CFastaFile {
	FILE *file;
	char *io_buffer;
    
    uchar * compressedSeq;
    

public:
	int64 file_size;
	t_fa_mode mode;

	uchar *raw_data;
	int64 raw_data_pos;
	int64 trans_data_pos;
	uint32_t data_size;
	uchar *raw_data_ptr;
	uchar *trans_data_ptr;
	uchar *eof_data_ptr;

	vector<CFastaSequence> sequences;


	bool ReadSequence();
	//bool WriteSequenceToMemory(CFastaSequence &seq);
	//bool WriteSequenceToFile(CFastaSequence &seq);

public:
	CFastaFile() {
		file     = NULL;
		mode     = mode_fa_none;
		raw_data = NULL;
		io_buffer = NULL;
	};

	~CFastaFile() {
		if(file)
			fclose(file);
		if(raw_data)
			delete[] raw_data;
		if(io_buffer)
			delete[] io_buffer;
	}

	bool Open(string file_name);
	bool Create(string file_name);
	bool Close();
	bool Read(vector<CFastaSequence> &_sequences);
	int64 GetFileSize() {return file_size;}
    
    
    void compressSequence();
    void writeCompressedSequenceToFile(FILE * file);
    
    /******/
    void readData(FILE * in);
    uchar getRefFromCompressed(uint32_t n);
    
    uint32_t getDataSize(){return data_size; };
    
};

#endif
