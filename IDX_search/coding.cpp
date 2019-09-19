/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */

#include "coding.h"


void StoreUInt(unsigned char *dest, uint32_t  data)
{
    *dest++ = (data      ) & 0xFF;
    *dest++ = (data >>  8) & 0xFF;
    *dest++ = (data >> 16) & 0xFF;
    *dest++ = (data >> 24) & 0xFF;
    
}

void LoadUInt(unsigned char *dest, uint32_t  &data)
{
    
  // ver1
    data = 0;
    data = *dest++;
    data = data + ((*dest++) << 8);
    data = data + ((*dest++) << 16);
    data = data + ((*dest  ) << 24);
   
    
    //ver2
    //A_memcpy(&data, dest, sizeof(uint32_t));
    
    
    //ver3
    //memcpy(&data, dest, sizeof(uint32_t));

}




 int32_t  encodeSuffix(const char * kmer, int32_t len)
{
   int32_t encoded_kmer = 0;
   int32_t kmer_pos = 0;
   int32_t charToEncode = LUT_ENCODED_CHARS < len? LUT_ENCODED_CHARS : len ;
   if(charToEncode > kl)
        charToEncode = kl;
   int32_t diff = LUT_ENCODED_CHARS - charToEncode;
   int32_t tempCode = 0;

    
    while (charToEncode--)
    {
        //tempChar = kmer[kmer_pos++];
        
        switch (kmer[kmer_pos++]) {
            case 'A':
                tempCode = 0;
                break;
            case 'a':
                tempCode = 0;
                break;
            case 'C':
                tempCode = 1;
                break;
            case 'c':
                tempCode = 1;
                break;
            case 'G':
                tempCode = 2;
                break;
            case 'g':
                tempCode = 2;
                break;
            case 'T':
                tempCode = 3;
                break;
            case 't':
                tempCode = 3;
                break;
                
            default:
                cout << "There shouldn't be other characters than AaCcTtGg - error" << endl;
                break;
        }
        
        encoded_kmer = encoded_kmer | (tempCode << ((charToEncode+diff)<<1));
        
    }
    
    return encoded_kmer;
}

unsigned char encodeChar(char x)
{
    unsigned char code = 0xFF;
    switch(x)
    {
            
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
        case 'N':
            code = 5;
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
            printf("Unknown nucleotide: %c\n", x);
            break;
    }
    
    return code;
}

unsigned char getRefFinal(uint32_t n)
{
   
    if(n&1)
        return reference[n>>1] & 0xF;
    else
        return (reference[n>>1] >> 4) & 0x0F ;
}

unsigned char getRef(uint32_t n)
{
    unsigned char c;
    if(n&1)
        c = reference[n/2] & 0xF;
    else
        c = (reference[n/2] >> 4) & 0x0F ;
    
    
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

char decodeCode(unsigned char x)
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

void getReadCode(const char * read, unsigned char * readCode, int64_t len)
{
   int32_t i = 0;

    
    while(i < len)
    {
        switch(read[i])
        {
            case 'A':
                readCode[i] = 0;
                break;
            case 'C':
                readCode[i] = 1;
                break;
            case 'G':
                readCode[i] = 2;
                break;
            case 'T':
                readCode[i] = 8;
                break;
            case 'N':  //different than in REF
                readCode[i] = 12;//5;
                break;
            case 'K':
                readCode[i] = 3;
                break;
            case 'M':
                readCode[i] = 4;
                break;
            case 'R':
                readCode[i] = 6;
                break;
            case 'S':
                readCode[i] = 7;
                break;
            case 'U':
                readCode[i] = 9;
                break;
            case 'W':
                readCode[i] = 10;
                break;
            case 'Y':
                readCode[i] = 11;
                break;
            default:
                printf("Unknown nucleotide: %c\n", read[i]);
                break;
        }
        
        ++i;
    }
    readCode[i] = 0xFF;
}

