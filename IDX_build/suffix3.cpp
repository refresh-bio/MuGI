/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */

#include "suffix3.h"


char *  getSuffix(suffix3 * suffix, unsigned int * nextPos)
{
    
    char * result = (char *) malloc ((Data::kl+2) * sizeof( char));
    int pos = 0;
    int max = Data::kl+1;
    bool firstVT = true;
    unsigned int chCovered = 0;
    unsigned int refPos = suffix->prevPos+1;
    int back, prevIns;
    
    ///BINARY SEARCH FOR FIRST VT
    unsigned int vtNo = 0, l = 0, r = Data::vl->get_no_var()-1, m;
    while (l < r)
    {
        m = (l+r) / 2;
        if ((Data::vl->vt[m].prePos < suffix->prevPos || Data::vl->vt[m].prePos == 0xFFFFFFFF) && suffix->prevPos != 0xFFFFFFFF)
        {
            l = m + 1;
        }
        else
        {
            r = m;
        }
    }
    vtNo = l;
    
    unsigned int vtCovered = 0, newRefPos;
    
    while (chCovered < Data::kl)
    {
        
        while(Data::vl->vt[vtNo].prePos+1 < refPos && chCovered < Data::kl && vtNo < Data::vl->get_no_var())
        {
            if(Data::vl->getType(vtNo) == DEL || Data::vl->getType(vtNo) == SV)
            {
                
                newRefPos = Data::vl->vt[vtNo].prePos+1+Data::vl->getDelLen(vtNo) > refPos ? Data::vl->vt[vtNo].prePos+1+Data::vl->getDelLen(vtNo) : refPos;
                if(newRefPos > refPos)
                {
                    if(suffix->vtList[vtCovered/NO_BITS_UINT] & (1 << (NO_BITS_UINT - vtCovered%NO_BITS_UINT - 1))) //variant found/introduced
                    {
                        refPos = newRefPos;
                        if(Data::vl->getInsLen(vtNo) != 0)
                        {
                            back = 1;
                            prevIns = 0;
                            while(Data::vl->vt[vtNo].prePos == Data::vl->vt[vtNo-back].prePos && (Data::vl->getType(vtNo-back) == INS || (Data::vl->getType(vtNo-back) == SV && Data::vl->getInsLen(vtNo-back) > 0)))
                            {
                                prevIns = (unsigned int)Data::vl->getInsLen(vtNo-back) > prevIns? (unsigned int)Data::vl->getInsLen(vtNo-back) : prevIns;
                                back++;
                            }
                            if(prevIns == 0 || Data::vl->getInsLen(vtNo) > prevIns)
                            {
                                if(firstVT)
                                {
                                    for (int c = suffix->gap; c < Data::vl->getInsLen(vtNo); c++)
                                    {
                                        result[pos++] = Coding::decodeCode(Data::vl->getPtrToFirstInsChar(vtNo)[c]);
                                        if (pos == max)
                                        {
                                            max = max * 2;
                                            char * tempResult = ( char *) realloc (result, (max+1) * sizeof( char));
                                            if (NULL == tempResult)
                                            {
                                                printf("Error in realloc\n");
                                                exit(8);
                                            }
                                            else
                                            {
                                                result = tempResult;
                                            }
                                        }
                                    }
                                    //cout << vt[vtNo].alt+suffix->gap;
                                    chCovered = chCovered + (unsigned int)Data::vl->getInsLen(vtNo)-suffix->gap;
                                    firstVT = false;
                                }
                                else
                                {
                                    for (int c = 0; c < Data::vl->getInsLen(vtNo); c++)
                                    {
                                        result[pos++] = Coding::decodeCode(Data::vl->getPtrToFirstInsChar(vtNo)[c]);
                                        if (pos == max)
                                        {
                                            max = max * 2;
                                            char * tempResult = ( char *) realloc (result, (max+1) * sizeof( char));
                                            if (NULL == tempResult)
                                            {
                                                printf("Error in realloc\n");
                                                exit(8);
                                            }
                                            else
                                            {
                                                result = tempResult;
                                            }
                                        }
                                    }
                                    //cout << vt[vtNo].alt+suffix->gap;
                                    chCovered = chCovered + Data::vl->getInsLen(vtNo);
                                }
                            }
                            
                            
                        }
                        
                    }
                    vtCovered++;
                }
                
                
            }
            vtNo++;
        }
        
        if(chCovered < Data::kl)
        {
            if(vtNo >= Data::vl->get_no_var() || Data::vl->vt[vtNo].prePos+1 > refPos)  //no variant at currPos
            {
                if(refPos < Data::fasta->data_size)
                {
                    result[pos++] = Data::fasta->raw_data[refPos++];
                    if (pos == max)
                    {
                        max = max * 2;
                        result = ( char *) realloc (result, (max+1) * sizeof( char ));
                    }
                    
                    //cout << data[0].raw_data[refPos++];
                    chCovered++;
                }
                else
                {
                    break;
                }
            }
            else
            {
                if(suffix->vtList[vtCovered/NO_BITS_UINT] & (1 << (NO_BITS_UINT - vtCovered%NO_BITS_UINT - 1))) //variant found/introduced
                {
                    switch (Data::vl->getType(vtNo))
                    {
                        case SNP:
                            //introduce variant
                            result[pos++] = Coding::decodeCode(Data::vl->vt[vtNo].alt);
                            if (pos == max)
                            {
                                max = max * 2;
                                result = ( char *) realloc (result, (max+1) * sizeof( char));
                            }
                            
                            chCovered++;
                            refPos++;
                            vtCovered++;
                            vtNo++;
                            firstVT = false;
                            break;
                        case DEL:
                            //introduce variant
                            refPos = refPos + Data::vl->getDelLen(vtNo);
                            vtCovered++;
                            vtNo++;
                            firstVT = false;
                            break;
                        case INS:
                            //introduce variant
                            if(firstVT)
                            {
                                for (int c = suffix->gap; c < Data::vl->getInsLen(vtNo); c++)
                                {
                                    result[pos++] = Coding::decodeCode(Data::vl->getPtrToFirstInsChar(vtNo)[c]);
                                    if (pos == max)
                                    {
                                        max = max * 2;
                                        result = ( char *) realloc (result, (max+1) * sizeof( char));
                                    }
                                }
                                //cout << vt[vtNo].alt+suffix->gap;
                                chCovered = chCovered + Data::vl->getInsLen(vtNo)-suffix->gap;
                                firstVT = false;
                            }
                            else
                            {
                                for (int c = 0; c < Data::vl->getInsLen(vtNo); c++)
                                {
                                    result[pos++] = Coding::decodeCode(Data::vl->getPtrToFirstInsChar(vtNo)[c]);
                                    if (pos == max)
                                    {
                                        max = max * 2;
                                        result = ( char *) realloc (result, (max+1) * sizeof( char));
                                    }
                                }
                                //cout << vt[vtNo].alt+suffix->gap;
                                chCovered = chCovered + Data::vl->getInsLen(vtNo);
                            }
                            vtCovered++;
                            vtNo++;
                            
                            break;
                        case SV:
                            //introduce variant
                            if(Data::vl->getInsLen(vtNo) == 0)
                            {
                                refPos = refPos + Data::vl->getDelLen(vtNo);
                                vtCovered++;
                                vtNo++;
                                firstVT = false;
                            }
                            else
                            {
                                if(firstVT)
                                {
                                    for (int c = suffix->gap; c < Data::vl->getInsLen(vtNo); c++)
                                    {
                                        result[pos++] = Coding::decodeCode(Data::vl->getPtrToFirstInsChar(vtNo)[c]);
                                        if (pos == max)
                                        {
                                            max = max * 2;
                                            result = ( char *) realloc (result, (max+1) * sizeof( char));
                                        }
                                    }
                                    
                                    chCovered = chCovered + Data::vl->getInsLen(vtNo)-suffix->gap;
                                    firstVT = false;
                                }
                                else
                                {
                                    for (int c = 0; c < Data::vl->getInsLen(vtNo); c++)
                                    {
                                        result[pos++] = Coding::decodeCode(Data::vl->getPtrToFirstInsChar(vtNo)[c]);
                                        if (pos == max)
                                        {
                                            max = max * 2;
                                            result = ( char *) realloc (result, (max+1) * sizeof( char));
                                        }
                                    }
                                    //cout << vt[vtNo].alt+suffix->gap;
                                    chCovered = chCovered + Data::vl->getInsLen(vtNo);
                                }
                                refPos = refPos + Data::vl->getDelLen(vtNo);
                                vtCovered++;
                                vtNo++;
                            }
                            break;
                    }
                    
                }
                else //variant not found/not introduced
                {
                    vtCovered++;
                    vtNo++;
                }
            }
        }
    }
    
    result[pos++] = '\0';
    if (pos == max)
    {
        max = max * 2;
        result = ( char *) realloc (result, (max+1) * sizeof( char));
    }
    //cout << endl;
    
    
    *nextPos = refPos;
    return result;
}

char  getFirstChar(suffix3 * suffix)
{
    
    char ch = 0;
    
    bool firstVT = true;
    unsigned int chCovered = 0;
    unsigned int refPos = suffix->prevPos+1;
    int back, prevIns;
    
    ///BINARY SEARCH FOR FIRST VT
    unsigned int vtNo = 0, l = 0, r = Data::vl->get_no_var()-1, m;
    while (l < r)
    {
        m = (l+r) / 2;
        if ((Data::vl->vt[m].prePos < suffix->prevPos || Data::vl->vt[m].prePos == 0xFFFFFFFF) && suffix->prevPos != 0xFFFFFFFF)
        {
            l = m + 1;
        }
        else
        {
            r = m;
        }
    }
    vtNo = l;
    
    unsigned int vtCovered = 0, newRefPos;
    
    while (!ch)
    {
        
        while(Data::vl->vt[vtNo].prePos+1 < refPos && chCovered < Data::kl && vtNo < Data::vl->get_no_var())
        {
            if(Data::vl->getType(vtNo) == DEL || Data::vl->getType(vtNo) == SV)
            {
                newRefPos = Data::vl->vt[vtNo].prePos+1+Data::vl->getDelLen(vtNo) > refPos ? Data::vl->vt[vtNo].prePos+1+Data::vl->getDelLen(vtNo) : refPos;
                if(newRefPos > refPos)
                {
                    if(suffix->vtList[vtCovered/NO_BITS_UINT] & (1 << (NO_BITS_UINT - vtCovered%NO_BITS_UINT - 1))) //variant found/introduced
                    {
                        refPos = newRefPos;
                        if(Data::vl->getInsLen(vtNo) != 0)
                        {
                            back = 1;
                            prevIns = 0;
                            while(Data::vl->vt[vtNo].prePos == Data::vl->vt[vtNo-back].prePos && (Data::vl->getType(vtNo-back) == INS || (Data::vl->getType(vtNo-back) == SV && Data::vl->getInsLen(vtNo-back) > 0)))
                            {
                                prevIns = (unsigned int)Data::vl->getInsLen(vtNo-back) > prevIns? (unsigned int)Data::vl->getInsLen(vtNo-back) : prevIns;
                                back++;
                            }
                            if(prevIns == 0 || Data::vl->getInsLen(vtNo) > prevIns)
                            {
                                if(firstVT)
                                {
                                    if ( suffix->gap < Data::vl->getInsLen(vtNo))
                                        return Coding::decodeCode(Data::vl->getPtrToFirstInsChar(vtNo)[suffix->gap]) ;//
                                    
                                    firstVT = false;
                                    vtCovered++;
                                    vtNo++;
                                    
                                }
                                else
                                {
                                    return  Coding::decodeCode(Data::vl->getPtrToFirstInsChar(vtNo)[0]);
                                }
                            }
                            
                            
                        }
                        
                    }
                    vtCovered++;
                }
                
                
            }
            vtNo++;
        }
        
        if(!ch)
        {
            if(vtNo >= Data::vl->get_no_var() || Data::vl->vt[vtNo].prePos+1 > refPos)  //no variant at currPos
            {
                if(refPos < Data::fasta->data_size)
                {
                    if(refPos >= Data::fasta->data_size)
                        return '\0';
                    else
                        return Data::fasta->raw_data[refPos];
                }
            }
            else
            {
                if(suffix->vtList[vtCovered/NO_BITS_UINT] & (1 << (NO_BITS_UINT - vtCovered%NO_BITS_UINT - 1))) //variant found/introduced
                {
                    switch (Data::vl->getType(vtNo))
                    {
                        case SNP:
                            //introduce variant
                            return Coding::decodeCode(Data::vl->vt[vtNo].alt);
                            
                            break;
                        case DEL:
                            //introduce variant
                            refPos = refPos + Data::vl->getDelLen(vtNo);
                            vtCovered++;
                            vtNo++;
                            firstVT = false;
                            break;
                        case INS:
                            //introduce variant
                            if(firstVT)
                            {
                                if ( suffix->gap < Data::vl->getInsLen(vtNo))
                                    return Coding::decodeCode(Data::vl->getPtrToFirstInsChar(vtNo)[suffix->gap]) ;//
                                
                                firstVT = false;
                                vtCovered++;
                                vtNo++;
                            }
                            else
                            {
                                
                                return Coding::decodeCode(Data::vl->getPtrToFirstInsChar(vtNo)[0]);
                            }
                            
                            break;
                        case SV:
                            //introduce variant
                            if(Data::vl->getInsLen(vtNo) == 0)
                            {
                                refPos = refPos + Data::vl->getDelLen(vtNo);
                                vtCovered++;
                                vtNo++;
                                firstVT = false;
                            }
                            else
                            {
                                if(firstVT)
                                {
                                    if ( suffix->gap < Data::vl->getInsLen(vtNo))
                                        return Coding::decodeCode(Data::vl->getPtrToFirstInsChar(vtNo)[suffix->gap]) ;//
                                    
                                    firstVT = false;
                                    vtCovered++;
                                    vtNo++;
                                    
                                }
                                else
                                {
                                    
                                    return Coding::decodeCode(Data::vl->getPtrToFirstInsChar(vtNo)[0]);
                                }
                              
                            }
                            break;
                    }
                    
                }
                else //variant not found/not introduced
                {
                    vtCovered++;
                    vtNo++;
                }
            }
        }
    }
    
    return ch;
}

bool compareSuffix3(const suffix3& s1, const suffix3& s2)
{
    char * kmer1, *kmer2;
    
    
    unsigned int nextPos1, nextPos2;
    bool result = true;
    
    kmer1 = getSuffix((suffix3*) &s1, &nextPos1);
    kmer2 = getSuffix((suffix3*) &s2, &nextPos2);
    
    int kmerComp = strncmp((const char*) kmer1, (const char*) kmer2, Data::kl);
    
   

    
    
    if (kmerComp < 0)
        result = true;
    else //if (kmerComp > 0)
        result = false;
    
    
    free(kmer1);
    free(kmer2);
    return result;
}


