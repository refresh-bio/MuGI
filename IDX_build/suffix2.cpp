/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */

#include "suffix2.h"


char *  getSuffix(suffix2 * suffix, unsigned int * nextPos)
{
    
    char * result = ( char *) malloc ((Data::kl+2)* sizeof( char));
    int pos = 0;
    int max = (Data::kl+1);
    bool firstVT = true;
    unsigned int chCovered = 0;
    unsigned int refPos = suffix->prevPos+1;
    
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
    
    
    
    unsigned int vtCovered = 0;
    
    
    while (chCovered < Data::kl)
    {
        while(vtNo < Data::vl->get_no_var() && Data::vl->vt[vtNo].prePos+1 < refPos)
        {
           
            
            cout << "shouldn't end up here with suffix2 (no var in del)" << suffix->prevPos << " " << suffix->vtList << " " << suffix->gap << " " << refPos << " " <<  vtNo << " " << Data::vl->vt[vtNo].prePos << endl;
            
             vtNo++;
        }
        
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
            if(suffix->vtList & (1 << (NO_BITS_UINT - vtCovered - 1))) //variant found/introduced
            {
                switch (Data::vl->getType(vtNo))
                {
                    case SNP:
                        //introduce variant
                        result[pos++] =  Coding::decodeCode(Data::vl->vt[vtNo].alt);
                        if (pos == max)
                        {
                            max = max * 2;
                            result = ( char *) realloc (result, (max+1) * sizeof( char));
                        }
                        
                        chCovered++;
                        refPos++;
                        vtCovered++;
                        vtNo++;
                        
                        break;
                    case DEL:
                        //introduce variant
                        refPos = refPos + Data::vl->getDelLen(vtNo);
                        vtCovered++;
                        vtNo++;
                        
                        break;
                    case INS:
                        //introduce variant
                        if(firstVT)
                        {
                            for (int c = suffix->gap; c < Data::vl->getInsLen(vtNo); c++)
                            {
                                result[pos++] =  Coding::decodeCode(Data::vl->getPtrToFirstInsChar(vtNo)[c]);
                                if (pos == max)
                                {
                                    max = max * 2;
                                    result = ( char *) realloc (result, (max+1) * sizeof( char));
                                }
                            }
                            chCovered = chCovered + Data::vl->getInsLen(vtNo) - suffix->gap;
                            firstVT = false;
                        }
                        else
                        {
                            for (int c = 0; c < Data::vl->getInsLen(vtNo); c++)
                            {
                                result[pos++] =  Coding::decodeCode(Data::vl->getPtrToFirstInsChar(vtNo)[c]);
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
    
    result[pos++] = '\0';
    
    
    *nextPos = refPos;
    return result;
}



char  getFirstChar(suffix2 * suffix)
{
    
    char ch = 0;
    
    bool firstVT = true;
    unsigned int refPos = suffix->prevPos+1;
    
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
   
    
    
    unsigned int vtCovered = 0;
    
    
    while (!ch)
    {
        while(vtNo < Data::vl->get_no_var() && Data::vl->vt[vtNo].prePos+1 < refPos)
        {
            vtNo++;
            
            cout << "[getfirstchar] shouldn't end up here with suffix2 (no var in del)" << suffix->prevPos << " " << suffix->vtList << " " << suffix->gap << " " << refPos << " " <<  vtNo << " " << Data::vl->vt[vtNo].prePos << endl;
        }
        
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
            if(suffix->vtList & (1 << (NO_BITS_UINT - vtCovered - 1))) //variant found/introduced
            {
                switch (Data::vl->getType(vtNo))
                {
                    case SNP:
                         return  Coding::decodeCode(Data::vl->vt[vtNo].alt);
                        break;
                    case DEL:
                        //introduce variant
                        refPos = refPos + Data::vl->getDelLen(vtNo);
                        vtCovered++;
                        vtNo++;
                        
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
                            return   Coding::decodeCode(Data::vl->getPtrToFirstInsChar(vtNo)[0]);
                        }
                        
                        break;
                    case SV:
                        //introduce variant
                        if(Data::vl->getInsLen(vtNo) == 0)
                        {
                            refPos = refPos + Data::vl->getDelLen(vtNo);
                            vtCovered++;
                            vtNo++;
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
                                return   Coding::decodeCode(Data::vl->getPtrToFirstInsChar(vtNo)[0]);
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
    
      return ch;
}


bool compareSuffix2(const suffix2& s1, const suffix2& s2)
{
    char * kmer1, *kmer2;
    unsigned int nextPos1, nextPos2;
    bool result = true;
    
    kmer1 = getSuffix((suffix2*) &s1, &nextPos1);
    kmer2 = getSuffix((suffix2*) &s2, &nextPos2);
    
    int kmerComp = strncmp((const char*) kmer1, (const char*) kmer2, Data::kl);
    if (kmerComp < 0)
        result = true;
    else //if (kmerComp > 0)
        result = false;
       
    free(kmer1);
    free(kmer2);
    return result;
}