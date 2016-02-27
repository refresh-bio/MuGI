/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */

#include "defs.h"
#include "fasta.h"
#include <algorithm>

using namespace std;

const uint32 io_buffer_size = 1 << 16;

// ********************************************************************************************
// CFastaSequence
// ********************************************************************************************
CFastaSequence::CFastaSequence()
{
	start_pos = 0;
	size = 0;
}

// ********************************************************************************************
CFastaSequence::CFastaSequence(uchar *_raw_data, string &_seq_name, string &_chrom, uint32 _line_len, uint32 _start_pos, uint32 _size, uchar _eol_type)
{
	seq_name  = _seq_name;
	line_len  = _line_len;
	size      = _size;
	start_pos = _start_pos;
	eol_type  = _eol_type;
    chrom = _chrom;
}

// ********************************************************************************************
CFastaSequence::CFastaSequence(const CFastaSequence &x)
{
	seq_name  = x.seq_name;
    chrom     = x.chrom;
	start_pos = x.start_pos;
	size      = x.size;
	line_len  = x.line_len;
	eol_type  = x.eol_type;
}

// ********************************************************************************************
CFastaSequence::~CFastaSequence()
{
}


// ********************************************************************************************
// CFastaFile
// ********************************************************************************************
bool CFastaFile::Open(string file_name)
{
	if(mode != mode_fa_none)
		return false;

	if((file = my_fopen(file_name.c_str(), "rb")) != NULL)
		mode = mode_fa_read;

	if(mode != mode_fa_read)
		return false;

	my_fseek(file, 0, SEEK_END);
	file_size = my_ftell(file);
	my_fseek(file, 0, SEEK_SET);

	return true;
}

// ********************************************************************************************
bool CFastaFile::Create(string file_name)
{
	if((file = my_fopen(file_name.c_str(), "wb")) != NULL)
		mode = mode_fa_write;
	file_size = 0;

	io_buffer = new char[io_buffer_size];
	setvbuf(file, io_buffer, _IOFBF, io_buffer_size);

	return mode == mode_fa_write;
}

// ********************************************************************************************
bool CFastaFile::Close()
{
	if(!file)
		return false;

	fclose(file);
	file = NULL;
	mode = mode_fa_none;

	if(raw_data)
		delete[] raw_data;
	raw_data = NULL;

	return true;
}

// ********************************************************************************************
bool CFastaFile::ReadSequence()
{
	string seq_name, chrom;
	int32 line_len;
	int size;
	int32 i;
	uchar c;
	uchar *start_ptr = trans_data_ptr;
	uchar eol_type;

	if(mode != mode_fa_read)
		return false;
	
	line_len = 0;

	c = *raw_data_ptr++;

	if(c != '>' && c != '@')
	{
		size = 0;
		return false;
	}
	if(raw_data_ptr >= eof_data_ptr)
		return false;

	seq_name = "";
	while(true)
	{
		c = *raw_data_ptr++;
		if(c == '\n' || c == '\r')
			break;
		seq_name.push_back(c);
	}

    
	eol_type = c;
	if(*raw_data_ptr == '\n' || *raw_data_ptr == '\r')
	{
		++raw_data_ptr;
		eol_type = 255;				// 0x0d + 0x0a
	}

    for(i = 0; ; ++i)
	{
		c = seq_name[i];
        if(c == ' ' || c == eol_type )
            break;
        
        chrom.push_back(c);

        
    }
    
    
	uchar *tmp_data_ptr = raw_data_ptr;
	for(i = 0; ; ++i)
	{
		c = *tmp_data_ptr++;
		if(c == '>' || c == '\n' || c == '\r' || c == '@')
			break;
	}
	line_len = i;

	while(true)
	{
		c = (uchar)*raw_data_ptr++;
		if(c >= 'A')
		{
			*trans_data_ptr++ = c;
		}
		else if(c == '>' || c == '@')
		{
			raw_data_ptr--;
			break;
		}
        
	}

	sequences.push_back(CFastaSequence(NULL, seq_name, chrom, line_len, (uint32) (start_ptr - raw_data), (uint32) (trans_data_ptr - start_ptr), eol_type));

	return true;
}

// ********************************************************************************************
bool CFastaFile::Read(vector<CFastaSequence> &_sequences)
{
	raw_data = new uchar[file_size+1];
	raw_data_pos   = 0;
	trans_data_pos = 0;

	raw_data_ptr = raw_data;
	trans_data_ptr = raw_data;
	eof_data_ptr = raw_data+file_size;

    
    
    
    int i = 0;
    
    
    //to work on all mac os -> bug in fread (some mac os), which cannot read more than 2GB at once
    uint64 temp_size = file_size ;
    while (temp_size > 2147483647)
    {
        fread(raw_data+i*2147483647, 1, 2147483647, file);
        temp_size = temp_size - 2147483647;
        i++;
	}
    fread(raw_data+i*2147483647, 1, temp_size, file);
          
          
    raw_data[file_size] = '>';

	sequences.clear();
	while(ReadSequence())
		;

    *trans_data_ptr = '\0';	
	

	_sequences = sequences;
    
  	data_size = sequences.back().start_pos + sequences.back().size;

    
   
	return true;
}



void CFastaFile::compressSequence()
{
    uint32_t cnt = 0,  i = 0;
    uchar code;
    compressedSeq = (uchar *) malloc((data_size/2+data_size%2) * sizeof(uchar));
    char c;
    while ((c = raw_data[i++]))
    {
        switch(c)
        {
            case '>':
                while ((c = raw_data[i++]) != '\n')
                    ;
                break;
            case '\n':
                break;
            default:
                code = Coding::encodeChar(c);
                if((cnt&1) == 0)
                {
                    code = code << 4;
                    compressedSeq[cnt/2] = code;
                }
                else
                    compressedSeq[cnt/2] = compressedSeq[cnt/2] | code;
                
                cnt++;
                
            break;
        }
    }
    if(cnt&1)
        compressedSeq[cnt/2] = compressedSeq[cnt/2] | 0xF;
       
}

void CFastaFile::writeCompressedSequenceToFile(FILE * file)
{
 
    uint32_t size = data_size/2 + data_size%2;
    fwrite(&size, sizeof(uint32_t), 1, file);
    fwrite(compressedSeq, 1, size, file);
    
}




/*************/
void CFastaFile::readData(FILE * in)
{
    uint32_t compRefSize;
    fread(&compRefSize,  sizeof(uint32_t), 1, in);
    compressedSeq = (uchar *)malloc((uint32_t)(compRefSize+1) * sizeof(uchar));
    fread(compressedSeq, 1, compRefSize, in);
    data_size = compRefSize * 2;
    if((compressedSeq[compRefSize-1] & 0xF) == 0xF)
        data_size--;//cout << "last EMPTY"; //last empty -ok
    else
        compressedSeq[compRefSize] = 0xFF; //add sentinel
    
    
    raw_data = (uchar *)malloc((uint32_t)(data_size) * sizeof(uchar));
    for(uint32_t i = 0; i < data_size; i++)
        raw_data[i] = getRefFromCompressed(i);

    
}


uchar CFastaFile::getRefFromCompressed(uint32_t n)
{
    unsigned char c;
    //if(n%2)
    if(n&1)
        c = compressedSeq[n/2] & 0xF;
    else
        c = (compressedSeq[n/2] >> 4) & 0x0F ;
    
    
    switch(c)
    {
        case 0:
            c = 'A';
            break;
        case 1:
            c = 'C';
            break;
        case 2:
            c = 'G';
            break;
        case 8:
            c = 'T';
            break;
        case 5:
            c = 'N';
            break;
        case 3:
            c = 'K';
            break;
        case 4:
            c = 'M';
            break;
        case 6:
            c = 'R';
            break;
        case 7:
            c = 'S';
            break;
        case 9:
            c = 'U';
            break;
        case 10:
            c = 'W';
            break;
        case 11:
            c = 'Y';
            break;
        default:
            printf(" Unknown nucleotide: %d\n", c);
            c = 0XF;
            break;
            
    }
    return c;
}


