/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */

#include "coding.h"



uchar Coding::encodeChar(char x)
{
    uchar code = 0xFF;
    switch(x)
    {
        case 'N':
            code = 5;
            break;
        case 'A':
            code = 0;
            break;
        case 'C':
            code = 1;
            break;
        case 'G':
            code = 2;
            break;
        case 'T':
            code = 8;
            break;
        case 'K':
            code = 3;
            break;
        case 'M':
            code = 4;
            break;
        case 'R':
            code = 6;
            break;
        case 'S':
            code = 7;
            break;
        case 'U':
            code = 9;
            break;
        case 'W':
            code = 10;
            break;
        case 'Y':
            code = 11;
            break;
        default:
            printf("Unknown nucleotide: %c, assigning 'N'\n", x);
            code = 5;
            break;
    }
    
    return code;
}


char Coding::decodeCode(uchar x)
{
    char c;
    switch(x)
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
            printf("decodeCode() Unknown nucleotide: %d\n", x);
            c = 0XF;
            break;
            
    }
    
    return c;
    
}


int Coding::encodeSuffix(char * kmer)
{
    int encoded_kmer = 0;
    int kmer_pos = 0;
    int charToEncode = LUT_ENCODED_CHARS;
    int tempCode = 0;
    char tempChar;
    
    if((int)strlen(kmer) < charToEncode)
        charToEncode = (int)strlen(kmer);
    
    int32_t diff = LUT_ENCODED_CHARS - charToEncode;

    
    
    
    while (charToEncode--)
    {
        tempChar = kmer[kmer_pos++];
        if(tempChar == 'A' || tempChar == 'a')
            tempCode = 0;
        else if (tempChar == 'C' || tempChar == 'c')
            tempCode = 1;
        else if (tempChar == 'G' || tempChar == 'g')
            tempCode = 2;
        else if (tempChar == 'T' || tempChar == 't')
            tempCode = 3;
        else
        {
            return -1;
        }
        
        encoded_kmer = encoded_kmer | (tempCode << (2*(charToEncode+diff)));
        
    }
    
    return encoded_kmer;
}




void Coding::StoreUInt(unsigned char *dest, uint32_t  data)
{
  
    *dest++ = (data      ) & 0xFF;
    *dest++ = (data >>  8) & 0xFF;
    *dest++ = (data >> 16) & 0xFF;
    *dest++ = (data >> 24) & 0xFF;
    
    
    
}


void Coding::LoadUInt(unsigned char *dest, uint32_t  &data)
{

    
    
    data = 0;
    data = *dest++;
    data = data + (*dest++ << 8);
    data = data + (*dest++ << 16);
    data = data + (*dest++ << 24);
    
}


