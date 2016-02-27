/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */
#include "seq_file.h"

bool SeqFile::Open(string file_name)
{
	
    
	file.open(file_name.c_str());
    
    if(!file.good())
    {
        return false;
    }
    
    int c = file.get();
    
    
    if(c == '@')
        mode = fastq;
    else if(c == '>')
        mode = fasta;
    else
        mode = simlple_list;
    
    
    file.unget();
    
	return true;
}

const char * SeqFile::getNextSeq()
{
    
    
    
    
    switch (mode)
    {
        case fastq:
            if(!getline(file, currName) || currName.length() == 0)
                return NULL;
            getline(file, currSeq);
            getline(file, currName2);
            getline(file, currQual);
            break;
        case fasta:
            if(!getline(file, currName) || currName.length() == 0)
                return NULL;
            getline(file, currSeq);
            break;
        case simlple_list:
            if(!getline(file, currSeq) || currSeq.length() == 0)
                return NULL;
            break;
        case none:
            return NULL;
            break;
    }
    
    
    return currSeq.c_str();
    
}

char * SeqFile::getCurrName()
{
    
    
    if(currName.size() >= maxSeqLen)
    {
        char * tmp;
        maxSeqLen = maxSeqLen * 2 > currName.size() ?  maxSeqLen * 2 : currName.size();
        tmp = (char *) realloc(seq, maxSeqLen*sizeof(char));
        if(tmp)
            seq = tmp;
        else
            return NULL;
    }
    seq[0] = '\0';
    memcpy(seq,currName.c_str(),currName.size());
    
    return seq;
}


bool SeqFile::Close()
{
	if(!file)
		file.close();
    mode = none;
    
    return true;
}
