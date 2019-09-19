/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */


#include <cstdio>
#include <math.h>
#include <stdlib.h>
#include <fcntl.h>
#include <algorithm>
#include <cstring>
#include <string>

#include <iostream>


#include "stat.h"
#include "timer.h"

#include "coding.h"
#include "SuffixArray.h"
#include "SuffixArray1.h"
#include "SuffixArray2.h"
#include "SuffixArray3.h"

#include "seq_file.h"


uint32_t ileSame = 0;

uint64_t *resV;
uint32_t * vtList;
uint32_t max_vtToCheck, vtListSize;
uint32_t * vtToCheck;
uint64_t resSize;

char * altVt;
using namespace std;
unsigned char * ins_sv_arr;

#define NO_BITS_UINT 32U

unsigned char* reference;
uint32_t refSize;

unsigned char *s_ch;

uint32_t kl; //k-mer length
uint32_t n_sparse;
FILE *res, *SA, *SA1, *SA2, *SA3;
uint32_t *LUT, *LUT1, *LUT2, *LUT3, *VT_LUT;


singleDelFromRight * delFromRight;
variant * vt;
uchar * vtType;
uint32_t lastKmerNoVar = 0;
uint32_t noVar=0, filtered = 0;
int32_t inBitV, ploidity;
char  **inputBuffer;
char *inputBuffer_t;
char **indNames;
uint64_t **vtVector;



uint32_t parts, blocking, cnt_uniq_parts, cnt_uniq_var, cnt_var;
uint32_t *all_var;
uint32_t *uniq_var;
uint64_t *uniq_parts;


suffix1 * suffixArray1;
suffix2 * suffixArray2;
suffix3 * suffixArray3;
uint32_t * suffixArray;

uint32_t suffixArrayCount = 0;

uint32_t suffixArray1Count = 0;
uint32_t suffixArray1Max = 0;
uint32_t suffixArray2Count = 0;
uint32_t suffixArray2Max = 0;
uint32_t suffixArray3Count = 0;
uint32_t suffixArray3MaxUInt = 0;

uint32_t maxDel = 0, maxSVdel = 0;
uint32_t * listSV, noSV = 0;
uint32_t * prevSV, * prevDeletion;


size_t my_fread(void *dest, size_t size, size_t count, FILE *file)
{
    size_t to_read = size * count;
    char *dest_char = (char *) dest;
    size_t readed = 0;
    
    while(to_read)
    {
        size_t block_size = 1 << 28;
        if(block_size > to_read)
            block_size = to_read;
        size_t block_readed = fread(dest_char, 1, block_size, file);
        if(block_size != block_readed)
            break;
        
        dest_char += block_readed;
        to_read -= block_readed;
        readed += block_readed;
    }
    
    return readed;
}


          

bool compareDelFromRight(singleDelFromRight a, singleDelFromRight b)
{
    if(a.nextRefPos < b.nextRefPos)
        return true;
    else if (a.nextRefPos > b.nextRefPos)
        return false;
    else
    {
        if(a.vtNo <= b.vtNo)
            return true;
        else
            return false;
    }
}




int32_t main(int32_t argc, const char * argv[])
{
    
    if(!argv[1] || !argv[2])
    {
        cout << "Wrong number of parameters." << endl;
        cout << "Required parameters:\n[index-file] [seq-file] {[maxError]}" << endl;
        cout << "where:\n[index-file] path to file with the index," << endl;
        cout << "[seq-file] path to file with the sequences to search (FASTA or FASTQ or plain sequences - one in each line)," << endl;
        cout << "[maxError] optional, maximum number of mismatches in found sequence (default 0, exact search)." << endl;
        exit(1);
    }
    
    
    // Get maximal number of error/mismatches in searched sequence
    int32_t maxError;
    if(argv[3])
    {
        for(int c = 0; c < strlen(argv[3]); c++)
        {
            if(argv[3][c] < '0' || argv[3][c] > '9')
            {
                cout << "Third parameter: \"" << argv[3] << "\" is not a number" << endl;
                cout << "Required parameters:\n[index-file] [seq-file] {[maxError]}" << endl;
                cout << "where:\n[index-file] path to file with the index," << endl;
                cout << "[seq-file] path to file with the sequences to search (FASTA or FASTQ or plain sequences - one in each line)," << endl;
                cout << "[maxError] optional, maximum number of mismatches in found sequence (default 0, exact search)." << endl;
                exit(1);
            }
        }
        maxError = atoi(argv[3]);
    }
    else
        maxError = DEFAULT_MAX_ERROR;
    cout << "maxError: " << maxError << endl;
    
    cout << "Reading data from: " << argv[1] << endl;
    FILE * in = fopen(argv[1], "rb");
    if (in == NULL) {
        printf("Cannot open %s\n", argv[1]);
        printf("The message is - %s\n", strerror(errno));
        exit(8);
    }
    
    uint32_t temp ; //for LUT_ENCODED_CHARS
    fread(&kl, 1, sizeof(uint32_t), in);
    fread(&n_sparse, 1, sizeof(uint32_t), in);
    fread(&temp, 1, sizeof(uint32_t), in);
    
    uint32_t compRefSize;
    fread(&compRefSize,  sizeof(uint32_t), 1, in);
    reference = (uchar *)malloc((uint32_t)(compRefSize+1) * sizeof(uchar));
    my_fread(reference, 1, compRefSize, in);
    refSize = compRefSize * 2;
    if((reference[compRefSize-1] & 0xF) == 0xF)
        refSize--;//cout << "last EMPTY"; //last empty -ok
    else
        reference[compRefSize] = 0xFF; //add sentinel
    
    
    uint32_t allChars = 0;
    uint32_t noIns = 0;
    uint32_t noDel = 0;
    
    fread(&noVar, 1, sizeof(uint32_t), in);
    fread(&noSV, 1, sizeof(uint32_t), in);
    fread(&noIns, 1, sizeof(uint32_t), in);
    fread(&noDel, 1, sizeof(uint32_t), in);
    fread(&allChars, 1, sizeof(uint32_t), in);
    
    
    fread(&maxSVdel, 1, sizeof(uint32_t), in);
    fread(&maxDel, 1, sizeof(uint32_t), in);
    
    vtType = (uchar *)malloc((noVar+1)*sizeof(uchar)); //add sentinel
    vt = (variant*) malloc((noVar+1)*sizeof(variant)); //add sentinel
    uint32_t size_of_ins_sv_arr = 2*sizeof(uint32_t)*noSV + sizeof(uint32_t)*noIns + allChars;
    ins_sv_arr = (unsigned char*) malloc((size_of_ins_sv_arr)*sizeof(unsigned char));
    
    
    fread(vtType, sizeof(uchar), noVar+1, in);
    fread(vt, sizeof(variant), noVar+1, in);
    fread(ins_sv_arr, sizeof(uchar), (2*sizeof(uint32_t)*noSV + sizeof(uint32_t)*noIns + allChars), in);
    
    
    prevSV = (uint32_t * ) malloc((noVar)*sizeof(uint32_t));
    prevDeletion = (uint32_t * ) malloc((noVar)*sizeof(uint32_t));
    listSV = (uint32_t *) malloc (noSV * sizeof(uint32_t));
    delFromRight = (struct singleDelFromRight * ) malloc (((noDel+noSV) + 1) * sizeof(struct singleDelFromRight));
    
    
    fread(listSV, sizeof(uint32_t), noSV, in);
    fread(prevSV, sizeof(uint32_t), noVar, in);
    fread(prevDeletion, sizeof(uint32_t), noVar, in);
    fread(delFromRight, sizeof(struct singleDelFromRight), (noDel+noSV) + 1, in);
    
 
    fread(&temp, sizeof(uint32_t), 1, in);
    VT_LUT = (uint32_t * ) malloc ((uint32_t)(refSize/VT_LUT_STEP+1) * sizeof(uint32_t));
    fread(VT_LUT, sizeof(uint32_t), refSize/VT_LUT_STEP, in);
    VT_LUT[refSize/VT_LUT_STEP] = noVar; //guard
    
    fread(&inBitV, sizeof(uint32_t), 1, in);
    fread(&ploidity, sizeof(uint32_t), 1, in);
    inBitV = inBitV * ploidity;
    
    fread(&blocking, sizeof(uint32_t), 1, in);
    fread(&parts, sizeof(uint32_t), 1, in);
    fread(&cnt_var, sizeof(uint32_t), 1, in);
    fread(&cnt_uniq_var, sizeof(uint32_t), 1, in);
    fread(&cnt_uniq_parts, sizeof(uint32_t), 1, in);
    
    all_var = new uint32_t[cnt_var];
    uniq_var = new uint32_t[cnt_uniq_var * parts];
    uniq_parts= new uint64_t[cnt_uniq_parts * blocking];
    
    
    fread(all_var, sizeof(uint32_t), cnt_var, in);
    my_fread(uniq_var, sizeof(uint32_t), cnt_uniq_var * parts, in);
    my_fread(uniq_parts, sizeof(unsigned long long), cnt_uniq_parts * blocking, in);
    
    
    fread(&suffixArrayCount,  sizeof(uint32_t),1, in);
  //  cout << "sa " << suffixArrayCount << endl;
    suffixArray = (uint32_t *) malloc ((uint64_t)suffixArrayCount * sizeof(uint32_t));
    my_fread(suffixArray, sizeof(uint32_t), suffixArrayCount, in);
    
    LUT = (uint32_t * ) malloc ((uint32_t)(LUT_MAX+1) * sizeof(uint32_t));
    fread(LUT, sizeof(unsigned int), LUT_MAX+1, in);
    
    fread(&suffixArray1Count,  sizeof(uint32_t),1, in);
   // cout <<"sa1 " <<  suffixArray1Count << endl;
    suffixArray1 = (suffix1 *) malloc (suffixArray1Count * sizeof(suffix1));
    my_fread(suffixArray1, sizeof(uint32_t), suffixArray1Count*2, in);
    
   
    LUT1 = (uint32_t * ) malloc ((uint32_t)(LUT_MAX+1) * sizeof(uint32_t));
    fread(LUT1, sizeof(uint32_t), LUT_MAX+1, in);
    
    
    
    
    fread(&suffixArray2Count, 1, sizeof(uint32_t), in);
  //  cout << "sa2 " <<  suffixArray2Count << endl;
    suffixArray2 = (suffix2 *) malloc (suffixArray2Count * sizeof(suffix2));
    fread(suffixArray2, sizeof(unsigned int), suffixArray2Count*3, in);
    
    
    
    
    
    LUT2 = (uint32_t * ) malloc ((uint32_t)(LUT_MAX+1) * sizeof(uint32_t));
    fread(LUT2, sizeof(uint32_t), LUT_MAX+1, in);
    
    fread(&suffixArray3Count,  sizeof(uint32_t),1, in);
  //  cout << "sa3 " <<  suffixArray3Count << endl;
    suffixArray3 = (suffix3 *) malloc (suffixArray3Count * sizeof(suffix3));
    for(uint32_t i = 0; i < suffixArray3Count; i++)
    {
        
        
        fread(&suffixArray3[i].prevPos, sizeof(unsigned int), 1, in);
        fread(&suffixArray3[i].gap, sizeof(unsigned int), 1, in);
        fread(&suffixArray3[i].vtListSize, sizeof(unsigned int), 1, in);
        suffixArray3[i].vtList = (uint32_t * ) malloc (suffixArray3[i].vtListSize * sizeof(uint32_t));
        for(int j = 0; j < suffixArray3[i].vtListSize; j++)
            fread(&suffixArray3[i].vtList[j], sizeof(unsigned int), 1, in);
        

        
    }
    LUT3 = (uint32_t * ) malloc ((uint32_t)(LUT_MAX+1) * sizeof(uint32_t));
    fread(LUT3, sizeof(uint32_t), LUT_MAX+1, in);
    
    
    cout << "total size: " << my_ftell(in) << endl;
    
    /************** SEARCH READS FROM EX-FILE *****************************************/
    
    uint32_t pos = 0, rememberPos;
  
    
    unsigned char read1[MAX_READ_SIZE];
    char subRead[MAX_READ_SIZE];
    
    // fflush(stdout);
    SDStat<double> stat_time;
    stat_time.reset();
    
    TimeMeasure timeMeasure;
    resSize = inBitV/64+(inBitV%64>0);
    
    
    uint64_t * positions, *resultAll, * positionsInResAll, * tempPos;
    uint32_t *encodedPattern;
    int32_t counter;
    //uint64_t foundPos;
    bool newPosition;
    
    encodedPattern = (uint32_t*) malloc(n_sparse * sizeof(uint32_t));
    if(!encodedPattern)
    {
        printf("Cannot allocate memory for encodedPattern, n_sparse: %d\n", n_sparse);
        printf("The message is - %s\n", strerror(errno));
        exit(8);
    }
    
    
    int32_t delta = 0;
    FILE *RESULT;
    RESULT = fopen("result.out", "wb");
    if(!RESULT)
    {
        printf("Cannot open result.out\n");
        printf("The message is - %s\n", strerror(errno));
        exit(8);
    }
    setvbuf(RESULT, NULL, _IOFBF, 1 << 28);
    
    unsigned char subRead1[MAX_READ_SIZE];
  
    
    int32_t readLen;
    int32_t a, r;
    
    
    
    uint32_t  searchLen;
    uint64_t posRes, maxRes;
    
    
    maxRes = 10000;
    uint64_t maxPos = 10000; 
    
    max_vtToCheck = 32;
    vtListSize = 1;
    resV = (uint64_t *)malloc(resSize*sizeof(uint64_t ));
    vtList = (uint32_t *) malloc (vtListSize * sizeof(uint32_t));
    vtToCheck = (uint32_t *) malloc ((max_vtToCheck)*sizeof(uint32_t));
    
    
    
    
    resultAll = (uint64_t*)malloc ((maxRes)*sizeof(uint64_t));
    positions = (uint64_t*) malloc ((maxPos)*sizeof(uint64_t));
    positionsInResAll = (uint64_t*)malloc ((maxPos)*sizeof(uint64_t));
    
    
    SeqFile seqFile;
    if(!seqFile.Open(argv[2]))
    {
        printf("Cannot open %s\n", argv[2]);
        printf("The message is - %s\n", strerror(errno));
        exit(8);
    }
    cout << "File with sequences to search: " << argv[2] << endl;

    
    const char *read;
    int32_t lenToCheck;
    uint32_t rightBorder, nameLen;
    bool all;
    uint32_t all_ones = 0xFFFFFFFF, position;
    int i = 0;
    if(maxError > 0)
    {
        read = seqFile.getNextSeq();
        
        while(read)
        {
            
            
            posRes = 0;
            
            
            
            readLen = (uint32_t)strlen(read);
            if( readLen <  n_sparse*(maxError+1))
            {
                cout << "Problem with query: " << read << "\nMinimum query length is equal to sparsity*(maxError+1) = " << n_sparse << "*(" << maxError<<"+1) = " << n_sparse*(maxError+1) << ", while the query is of length " << readLen << endl;;
                return(7);
            }

            
            getReadCode(read, read1, readLen);
            
            
            
            timeMeasure.reset_clock();
            timeMeasure.start_clock();
            
            
            
            searchLen = readLen/(maxError+1);
            
            
            for(int32_t e = 0; e <= maxError; ++e)
            {
                
                
                A_memcpy(subRead, (read+e*searchLen), (size_t)searchLen);
                subRead[searchLen] = '\0';
                
                A_memcpy(subRead1, (read1+e*searchLen), (size_t)searchLen);
                
                
                delta = 0;
                
                
                
                encodedPattern[delta] = encodeSuffix((char*)subRead+delta, searchLen-delta);
                lenToCheck = searchLen-delta < kl ? searchLen-delta : kl;
                rightBorder = LUT_ENCODED_CHARS > lenToCheck ? encodedPattern[delta] + (1<<((LUT_ENCODED_CHARS-lenToCheck)<<1)) : encodedPattern[delta]+1;
                
  
                
                searchString(subRead1, searchLen,  encodedPattern[delta], delta, &resultAll, &posRes, &maxRes, read1, readLen, e*searchLen-1, (e+1)*searchLen, maxError,  lenToCheck, rightBorder);
                
                searchSA1(subRead1, searchLen, encodedPattern[delta], delta, &resultAll, &posRes, &maxRes, read1, readLen, e*searchLen-1, (e+1)*searchLen, maxError,  lenToCheck, rightBorder);
                
                searchSA2(subRead1, searchLen, encodedPattern[delta], delta, &resultAll, &posRes, &maxRes, read1, readLen, e*searchLen-1, (e+1)*searchLen, maxError,  lenToCheck, rightBorder);
                
                searchSA3(subRead1, searchLen, encodedPattern[delta], delta, &resultAll, &posRes, &maxRes, read1, readLen, e*searchLen-1, (e+1)*searchLen, maxError,  lenToCheck, rightBorder);
                
                
                ++delta;
                
                while(delta < n_sparse)
                {
                    
                    encodedPattern[delta] = encodeSuffix((char*)subRead+delta, searchLen-delta);
                    lenToCheck = searchLen-delta < kl ? searchLen-delta : kl;
                    rightBorder = LUT_ENCODED_CHARS > lenToCheck ? encodedPattern[delta] + (1<<((LUT_ENCODED_CHARS-lenToCheck)<<1)) : encodedPattern[delta]+1;
           
                    searchString(subRead1, searchLen,  encodedPattern[delta], delta, &resultAll, &posRes, &maxRes, read1, readLen, e*searchLen-1, (e+1)*searchLen, maxError, lenToCheck, rightBorder);
                    searchSA1(subRead1, searchLen, encodedPattern[delta], delta, &resultAll, &posRes, &maxRes, read1, readLen, e*searchLen-1, (e+1)*searchLen, maxError,  lenToCheck, rightBorder);
                    searchSA3(subRead1, searchLen, encodedPattern[delta], delta, &resultAll, &posRes, &maxRes, read1, readLen, e*searchLen-1, (e+1)*searchLen, maxError,  lenToCheck, rightBorder);
                    ++delta;
                    
                }
                
            }
            
            
            if(posRes >= maxPos)
            {
                
                maxPos = 2 * maxPos > posRes ? 2 * maxPos : posRes+1 ;
                
                tempPos = (uint64_t *) realloc (positions, maxPos * sizeof( uint64_t));
                if (NULL == tempPos)
                {
                    printf("Error in realloc, the message is - %s\n max = %llu\n", strerror(errno), maxPos);
                    fflush(stdout);
                    //free(result);
                    exit(8);
                }
                else
                {
                    positions = tempPos;
                }
                
                
                tempPos = ( uint64_t *) realloc (positionsInResAll, maxPos * sizeof( uint64_t));
                if (NULL == tempPos)
                {
                    printf("Error in realloc, the message is - %s\n max = %llu\n", strerror(errno), maxPos);
                    fflush(stdout);
                    //free(result);
                    exit(8);
                }
                else
                {
                    positionsInResAll = tempPos;
                }
                
            }
            
            
            /* TO FILE */
            
            
            counter = 0;
            pos = 0;
            while(resultAll[pos++]>=0xFFFFFFFFFFFFFFFE)
            {
                
                newPosition = true;
                
                for(a = 0 ; a < counter; a++)
                    if(positions[a]==resultAll[pos])
                    {
                        newPosition = false;
                        break;
                    }
                if(newPosition)
                {
                    positions[counter] = resultAll[pos++];
                    positionsInResAll[counter++] = pos;
                    if(resultAll[pos-2] == 0xFFFFFFFFFFFFFFFF)
                        pos+=resSize;
                }
                else
                {
                    pos++;
                    if(resultAll[pos-2] == 0xFFFFFFFFFFFFFFFF && resultAll[positionsInResAll[a]-2] == 0xFFFFFFFFFFFFFFFF)
                        for( r=0; r < resSize; r++)
                            resultAll[positionsInResAll[a]+r]= resultAll[positionsInResAll[a]+r]|resultAll[pos++];
                    else
                    {
                        if(resultAll[positionsInResAll[a]-2] == 0xFFFFFFFFFFFFFFFF)
                            for( r=0; r < resSize; r++)
                                resultAll[positionsInResAll[a]+r]= 0xFFFFFFFFFFFFFFFF;
                    }
                }
            }
            
            
            timeMeasure.stop_clock();
            stat_time.add(timeMeasure.get_time());
            
            
            sort(positions, positions+counter);
            //move prePos-1 the beginning
            while(positions[counter-1] == 0xFFFFFFFF)
            {
                for(uint32 i = counter-1; i > 0; i--)
                    positions[i] = positions[i-1];
                positions[0] = 0xFFFFFFFE; //special to indicate -1 pos (other than "not found = 0xFFFFFFFF"
            }
            
            fwrite(&all_ones, sizeof(uint32_t), 1, RESULT);
            
            
            nameLen = (uint32_t)strlen((seqFile.getCurrName()));
            if(nameLen)
            {
                fwrite(&nameLen, sizeof(uint32_t), 1, RESULT);
                fwrite(seqFile.getCurrName(), sizeof(char), nameLen, RESULT);
            }
            else
            {
                fwrite(&all_ones, sizeof(uint32_t), 1, RESULT);
            }
            
            
            if(counter == 0)
                fwrite(&all_ones, sizeof(uint32_t), 1, RESULT);
            else
            {
                for( a=0; a<counter; a++)
                {
                    pos = 0;
                    while(resultAll[pos++]>=0xFFFFFFFFFFFFFFFE)
                    {
                        if(resultAll[pos]==positions[a] || (resultAll[pos]==0xFFFFFFFF && positions[a]==0xFFFFFFFE))
                        {
                            position = (uint32_t)positions[a];
                            pos++;
                            fwrite(&position, sizeof(uint32_t), 1, RESULT);
                            
                            if(resultAll[pos-2] == 0xFFFFFFFFFFFFFFFF)
                            {
                                rememberPos = pos;
                                all = true;
                                for( r=0; r < resSize; r++)
                                    if(resultAll[pos++] != 0xFFFFFFFFFFFFFFFF)
                                    {
                                        all = false;
                                        pos = rememberPos;
                                        break;
                                    }
                                
                                if(all)
                                    fwrite(&all_ones, sizeof(uint32_t), 1, RESULT);
                                else
                                {
                                    
                                    position = 0;
                                    fwrite(&position, sizeof(uint32_t), 1, RESULT);
                                    for( r=0; r < resSize; r++)
                                        fwrite(&resultAll[pos++], sizeof(uint64_t), 1, RESULT);
                                }
                            }
                            else
                            {
                                fwrite(&all_ones, sizeof(uint32_t), 1, RESULT);
                               
                            }
                            
                            break;
                        }
                        else
                        {
                            pos++;
                            if(resultAll[pos-2] == 0xFFFFFFFFFFFFFFFF)
                                pos = pos + (uint32_t)resSize;
                        }
                    }
                }
            }
            i++;
           read = seqFile.getNextSeq();
        }
    }
    else if(maxError == 0)
    {
        read = seqFile.getNextSeq();
        i = 0;
        while(read)
        {
            posRes = 0;
            delta = 0;
            
           
            
            readLen = (uint32_t)strlen(read);
            
            if( readLen < n_sparse*(maxError+1))
            {
                cout << "Problem with query: " << read << "\nMinimum query length is equal to sparsity*(maxError+1) = " << n_sparse << "*(" << maxError<<"+1) = " << n_sparse*(maxError+1) << ", while the query is of length " << readLen << endl;;
                return(7);
            }

            
            getReadCode(read, read1, readLen);
            
            
            
            timeMeasure.reset_clock();
            timeMeasure.start_clock();
            
            
            
            encodedPattern[delta] = encodeSuffix((char*)read+delta, readLen);
            lenToCheck = readLen-delta < kl ? readLen-delta : kl;
            rightBorder = LUT_ENCODED_CHARS > lenToCheck ? encodedPattern[delta] + (1<<((LUT_ENCODED_CHARS-lenToCheck)<<1)) : encodedPattern[delta]+1;
            
            searchStringExact(read1, readLen,  encodedPattern[delta], delta, &resultAll, &posRes, &maxRes, lenToCheck, rightBorder);
            
            
            searchSA1Exact(read1, readLen,  encodedPattern[delta], delta, &resultAll, &posRes, &maxRes, lenToCheck, rightBorder);
            
            searchSA2Exact(read1, readLen, encodedPattern[delta], delta, &resultAll, &posRes, &maxRes, lenToCheck, rightBorder);
            
            searchSA3Exact(read1, readLen,  encodedPattern[delta], delta, &resultAll, &posRes, &maxRes, lenToCheck, rightBorder);
            
            
            ++delta;
            
            while(delta < n_sparse)
            {
                
                encodedPattern[delta] = encodeSuffix((char*)read+delta, readLen);
                lenToCheck = readLen-delta < kl ? readLen-delta : kl;
                rightBorder = LUT_ENCODED_CHARS > lenToCheck ? encodedPattern[delta] + (1<<((LUT_ENCODED_CHARS-lenToCheck)<<1)) : encodedPattern[delta]+1;
            
                
                 searchStringExact(read1, readLen,  encodedPattern[delta], delta, &resultAll, &posRes, &maxRes, lenToCheck, rightBorder);
                
                
                searchSA1Exact(read1, readLen,  encodedPattern[delta], delta, &resultAll, &posRes, &maxRes, lenToCheck, rightBorder);
                
                //searchSA2Exact(read1, readLen, encodedPattern[delta], delta, &resultAll, &posRes, &maxRes, lenToCheck, rightBorder);
                
                searchSA3Exact(read1, readLen,  encodedPattern[delta], delta, &resultAll, &posRes, &maxRes, lenToCheck, rightBorder);
                
                ++delta;
                
            }
            
            
            
            
            if(posRes >= maxPos)
            {
                
                maxPos = 2 * maxPos > posRes ? 2 * maxPos : posRes+1 ;
                
                tempPos = ( uint64_t *) realloc (positions, maxPos * sizeof( uint64_t));
                if (NULL == tempPos)
                {
                    printf("Error in realloc, the message is - %s\n max = %llu\n", strerror(errno), maxPos);
                    fflush(stdout);
                    //free(result);
                    exit(8);
                }
                else
                {
                    positions = tempPos;
                }
                
                
                tempPos = ( uint64_t *) realloc (positionsInResAll, maxPos * sizeof( uint64_t));
                if (NULL == tempPos)
                {
                    printf("Error in realloc, the message is - %s\n max = %llu\n", strerror(errno), maxPos);
                    fflush(stdout);
                    //free(result);
                    exit(8);
                }
                else
                {
                    positionsInResAll = tempPos;
                }
                
            }
            
            
            /* TO FILE */
            
            
            counter = 0;
            pos = 0;
            while(resultAll[pos++]>=0xFFFFFFFFFFFFFFFE)
            {
                
                newPosition = true;
                
                for(a = 0 ; a < counter; a++)
                    if(positions[a]==resultAll[pos])
                    {
                        newPosition = false;
                        break;
                    }
                if(newPosition)
                {
                    positions[counter] = resultAll[pos++];
                    positionsInResAll[counter++] = pos;
                    if(resultAll[pos-2] == 0xFFFFFFFFFFFFFFFF)
                        pos+=resSize;//for(r=0; r < resSize; r++)
                    //    pos++; // fprintf(RESULT, "%llu ", resultSA[pos++]);
                }
                else
                {
                    pos++;
                    //resultAll[pos++] = 0xFFFFFFFFFFFFFFFF;
                    if(resultAll[pos-2] == 0xFFFFFFFFFFFFFFFF && resultAll[positionsInResAll[a]-2] == 0xFFFFFFFFFFFFFFFF)
                        for( r=0; r < resSize; r++)
                            resultAll[positionsInResAll[a]+r]= resultAll[positionsInResAll[a]+r]|resultAll[pos++];
                    else
                    {
                        if(resultAll[positionsInResAll[a]-2] == 0xFFFFFFFFFFFFFFFF)
                            for( r=0; r < resSize; r++)
                                resultAll[positionsInResAll[a]+r]= 0xFFFFFFFFFFFFFFFF;
                    }
                }
            }
            
            
            timeMeasure.stop_clock();
            stat_time.add(timeMeasure.get_time());
            
            
            sort(positions, positions+counter);
            //move prePos-1 the beginning
            while(positions[counter-1] == 0xFFFFFFFF)
            {
                for(uint32 i = counter-1; i > 0; i--)
                    positions[i] = positions[i-1];
                positions[0] = 0xFFFFFFFE;
            }
            
            fwrite(&all_ones, sizeof(uint32_t), 1, RESULT);
            
            nameLen = (uint32_t)strlen((seqFile.getCurrName()));
            if(nameLen)
            {
                fwrite(&nameLen, sizeof(uint32_t), 1, RESULT);
                fwrite(seqFile.getCurrName(), sizeof(char), nameLen, RESULT);
            }
            else
            {
                fwrite(&all_ones, sizeof(uint32_t), 1, RESULT);
            }
            
            if(counter == 0)
                fwrite(&all_ones, sizeof(uint32_t), 1, RESULT);
            else
            {
                for( a=0; a<counter; a++)
                {
                    pos = 0;
                    while(resultAll[pos++]>=0xFFFFFFFFFFFFFFFE)
                    {
                        if(resultAll[pos]==positions[a] || (resultAll[pos]==0xFFFFFFFF && positions[a]==0xFFFFFFFE))
                        {
                            position = (uint32_t)positions[a];
                            pos++;
                            fwrite(&position, sizeof(uint32_t), 1, RESULT);
                            
                            if(resultAll[pos-2] == 0xFFFFFFFFFFFFFFFF)
                            {
                                rememberPos = pos;
                                all = true;
                                for( r=0; r < resSize; r++)
                                    if(resultAll[pos++] != 0xFFFFFFFFFFFFFFFF)
                                    {
                                        all = false;
                                        pos = rememberPos;
                                        break;
                                    }
                                
                                if(all)
                                    fwrite(&all_ones, sizeof(uint32_t), 1, RESULT);
                                else
                                {
                
                                    position = 0;
                                    fwrite(&position, sizeof(uint32_t), 1, RESULT);
                                    for( r=0; r < resSize; r++)
                                        fwrite(&resultAll[pos++], sizeof(uint64_t), 1, RESULT);
                                }
                            }
                            else
                            {
                                fwrite(&all_ones, sizeof(uint32_t), 1, RESULT);
                              //  for(int32_t r=0; r < resSize; r++)
                              //      fwrite(&all_ones, sizeof(uint64_t), 1, RESULT);
                            }
                            
                            break;
                        }
                        else
                        {
                            pos++;
                            if(resultAll[pos-2] == 0xFFFFFFFFFFFFFFFF)
                                pos = pos + (uint32_t)resSize;
                        }
                    }
                }
            }
            read = seqFile.getNextSeq();
            i++;
            
        }
    }
    
    cout << i << " searches"<< endl;
    
    
    free(positions);
    free(positionsInResAll);
    free(resultAll);
    free(resV);
    
    stat_time.sort();
    printf("Average:%10.6f ms \nMedian: %10.6f ms  \n",   stat_time.get_average(), stat_time.get_median());//,
    
    float percent = 0;
    while (percent < 80)
    {
        printf("%6.2f%%\t%f\n", percent, stat_time.get_time_for_percent(percent));
        percent += 1;
    }
    
    while (percent <= 100)
    {
        printf("%6.2f%%\t%f\n", percent, stat_time.get_time_for_percent(percent));
        percent += 0.25;
    }
    
    
  //  stat_time.get_std_dev(stat_time.get_average());
    fflush(stdout);
    

    
    free(suffixArray);
    free(suffixArray1);
    free(suffixArray2);
    for( i = 0; i< suffixArray3Count; ++i)
        free(suffixArray3[i].vtList);
    free(suffixArray3);
    free(listSV);
    
    
    
    free(vt);
    

    
    
    free(VT_LUT);
    free(LUT);
    free(LUT1);
    free(LUT2);
    free(LUT3);
    
    
    fclose(RESULT);
    
    delete all_var ;
    delete uniq_var;
    delete uniq_parts;
    
    return 0;
}


uint64_t checkVTconfigurationExact(uint32_t prePos,  uint32_t noVt,  uint32_t firstVar, uint32_t currVt,  uint32_t lastVt,  uint32_t firstVarRead1, uint32_t refPos, uint32_t chCovered, unsigned char * read1, int32_t readLen, uint32_t delta,  int32_t preChToCheck,  uint64_t ** passed_result, uint64_t pos, uint64_t * max, uint32_t prevDelCount, uint32_t firstSV, uint32_t lastSV,  uint32_t  pos_vtToCheck,  uint32_t nextRefPos, uint32_t gap, uint32_t remainedToInsert, uint32_t prevDelVar)

{

    
    
    
    uint64_t * result = *passed_result;
    
    
    int32_t ch;
    
    
    uint32_t tempLen;
    
    uint32_t * temp_vtToCheck;
    bool goodPath, currPathValid = true;;
    uint32_t last = 0, len = 0, delLen = 0, insLen = 0;
    uint32_t newRefPos;
    
    uint32_t patternLen = readLen-delta;
    
    
    //
    if(remainedToInsert && chCovered < patternLen)
    {
        if(vtType[lastVt] == INS)
        {
            LoadUInt(ins_sv_arr+vt[lastVt].alt, insLen);
            tempLen = vt[lastVt].alt+sizeof(uint32_t);
        }
        else if(vtType[lastVt] == SV)//SV
        {
            LoadUInt(ins_sv_arr+vt[lastVt].alt+sizeof(uint32_t), insLen);
            tempLen = vt[lastVt].alt+(sizeof(uint32_t)<<1);
        }
        else
        {
            cout << "Should be SV or INS! " <<  lastVt << " " << prePos <<" " << remainedToInsert << "\n";
            exit(9);
        }
        
        
        
        while (chCovered < patternLen && remainedToInsert > 0)
        {
            if ((ins_sv_arr[tempLen+insLen-remainedToInsert]) == read1[chCovered+delta])
            {
                chCovered++; remainedToInsert--;
            }
            else
            {
                currPathValid = false;
                break;
            }
        }
        if(!remainedToInsert)
        {
            lastVt++;
        }
    }
    
    while (currPathValid && chCovered < patternLen)
    {
        while(vt[currVt].prePos+1 < refPos && chCovered < patternLen) //currVt < noVar &&
        {
            
            if(vtType[currVt] == DEL || vtType[currVt] == SV)
            {
                if(vtType[currVt] == DEL)
                {
                    delLen = vt[currVt].alt;
                    insLen = 0;
                }
                else
                {
                    LoadUInt(ins_sv_arr+vt[currVt].alt, delLen);
                    
                    
                    LoadUInt(ins_sv_arr+vt[currVt].alt+sizeof(uint32_t), insLen);
                }
                
                
                newRefPos = vt[currVt].prePos+1+delLen > refPos ? vt[currVt].prePos+1+delLen : refPos;
                
                if(newRefPos > refPos)
                {
                    ///
                    (vtToCheck)[pos_vtToCheck++] = currVt;
                    if(pos_vtToCheck == max_vtToCheck)
                    {
                        max_vtToCheck = max_vtToCheck * 2;
                        
                        temp_vtToCheck = ( uint32_t *) realloc (vtToCheck, max_vtToCheck * sizeof( uint32_t));
                        if (NULL == temp_vtToCheck)
                        {
                            printf("SA: Error in realloc, the message is - %s\n max = %d\n", strerror(errno), max_vtToCheck);
                            fflush(stdout);
                            //free(result);
                            exit(8);
                        }
                        else
                        {
                            vtToCheck = temp_vtToCheck;
                        }
                    }
                    ///
                    
                    
                    last = vtList[noVt/NO_BITS_UINT];
                    vtList[noVt/NO_BITS_UINT] = vtList[noVt/NO_BITS_UINT] | (1U << (NO_BITS_UINT -  (noVt&NO_BITS_UINT_MINUS_ONE) -  1U));
                    ++noVt;
                    if(noVt == vtListSize * NO_BITS_UINT)
                    {
                        vtListSize++;
                        uint32_t * tempList = (uint32_t *) realloc(vtList, vtListSize * sizeof(uint32_t));
                        if (NULL ==tempList) {
                            printf("Error in realloc, the message is - %s\n", strerror(errno));
                            exit(8);
                        }
                        else
                        {
                            vtList = tempList;
                        }
                        vtList[vtListSize-1]=0;
                    }
                    
                    
                    
                    if(insLen == 0)
                    {
                        pos = checkVTconfigurationExact(prePos,  noVt, firstVar, currVt+1, currVt+1, firstVarRead1, newRefPos, chCovered, read1, readLen, delta, preChToCheck,  &result, pos, max, prevDelCount, firstSV, lastSV, pos_vtToCheck, newRefPos, gap, remainedToInsert, prevDelVar);
                    }
                    else
                    {
                        
                        ///
                        //int32_t ch;
                        goodPath=true;
                        for(ch = 0; ch < insLen && chCovered+ch < patternLen; ch++)
                        {
                            if(read1[chCovered+ch+delta] != ins_sv_arr[vt[currVt].alt+(sizeof(uint32_t)<<1)+ch])
                            {
                                goodPath=false;
                                break;
                            }
                        }
                        if(goodPath)
                        {
                            if(insLen-ch == 0)
                                pos = checkVTconfigurationExact(prePos, noVt, firstVar, currVt+1,  currVt+1, firstVarRead1, newRefPos, chCovered+ch, read1, readLen, delta, preChToCheck, &result, pos, max, prevDelCount, firstSV, lastSV, pos_vtToCheck,   newRefPos, gap, insLen-ch,  prevDelVar);
                            else
                                pos = checkVTconfigurationExact(prePos,  noVt, firstVar, currVt+1,  currVt, firstVarRead1, newRefPos, chCovered+ch, read1, readLen, delta, preChToCheck, &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  newRefPos, gap, insLen-ch, prevDelVar);
                        }
                        
                        
                        ///
                    }
                    
                    vtList[(noVt-1)/NO_BITS_UINT] = last;
                    
                }
                
                
                
            }
            
            
            currVt++;
            lastVt++;
            
        }
        
        if( chCovered < patternLen)
        {
            if( vt[currVt].prePos+1 == refPos && refPos < refSize) //currVt < noVar &&
            {
                (vtToCheck)[(pos_vtToCheck)++] = currVt;
                if(pos_vtToCheck == max_vtToCheck)
                {
                    max_vtToCheck = max_vtToCheck * 2;
                    
                    temp_vtToCheck = ( uint32_t *) realloc (vtToCheck, max_vtToCheck * sizeof( uint32_t));
                    if (NULL == temp_vtToCheck)
                    {
                        printf("SA: Error in realloc, the message is - %s\n max = %d\n", strerror(errno),max_vtToCheck);
                        fflush(stdout);
                        //free(result);
                        exit(8);
                    }
                    else
                    {
                        vtToCheck = temp_vtToCheck;
                    }
                }
                ///
                
                
                
                last = vtList[noVt/NO_BITS_UINT];
                vtList[noVt/NO_BITS_UINT] = vtList[noVt/NO_BITS_UINT] | (1U << (NO_BITS_UINT -  (noVt&NO_BITS_UINT_MINUS_ONE) -  1U));
                ++noVt;
                
                if(noVt == vtListSize * NO_BITS_UINT)
                {
                    vtListSize++;
                    uint32_t * tempList = (uint32_t *) realloc(vtList, vtListSize * sizeof(uint32_t));
                    if (NULL ==tempList) {
                        printf("Error in realloc, the message is - %s\n", strerror(errno));
                        exit(8);
                    }
                    else
                    {
                        vtList = tempList;
                    }
                    vtList[vtListSize-1]=0;
                }
                
                
                
                switch (vtType[currVt])
                {
                    case SNP:
                        //introduce variant
                        
                        if(read1[chCovered+delta] == vt[currVt].alt)
                            pos = checkVTconfigurationExact(prePos,  noVt, firstVar, currVt+1,  currVt+1, firstVarRead1, refPos+1, chCovered+1, read1, readLen,delta, preChToCheck, &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck, refPos+1, gap, remainedToInsert, prevDelVar);
                        
                        break;
                    case DEL:
                        //introduce variant
                        
                        
                        
                        pos = checkVTconfigurationExact(prePos, noVt, firstVar, currVt+1,  currVt+1, firstVarRead1, refPos+vt[currVt].alt, chCovered, read1, readLen,delta, preChToCheck, &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  refPos+vt[currVt].alt, gap, remainedToInsert,  prevDelVar);
                        
                        
                        break;
                    case INS:
                        //introduce variant
                        
                       // int32_t ch;
                        goodPath=true;
                        LoadUInt(ins_sv_arr+vt[currVt].alt, insLen);
                        //A_memcpy(&insLen, ins_sv_arr+vt[currVt].alt, sizeof(uint32_t));
                        for(ch = 0; ch <  insLen && chCovered+ch < patternLen; ch++)
                        {
                            if(read1[chCovered+ch+delta] != ins_sv_arr[vt[currVt].alt+sizeof(uint32_t)+ch])
                            {
                                goodPath=false;
                                break;
                            }
                        }
                        if(goodPath)
                        {  if(insLen-ch == 0)
                            pos = checkVTconfigurationExact(prePos,  noVt, firstVar, currVt+1, currVt+1,  firstVarRead1, refPos, chCovered+ch, read1, readLen, delta, preChToCheck, &result, pos, max, prevDelCount, firstSV, lastSV, pos_vtToCheck,  refPos, gap, insLen-ch, prevDelVar);
                        else
                            pos = checkVTconfigurationExact(prePos,  noVt, firstVar, currVt+1, currVt, firstVarRead1,  refPos, chCovered+ch, read1, readLen, delta, preChToCheck, &result, pos, max, prevDelCount, firstSV, lastSV, pos_vtToCheck, refPos, gap, insLen-ch, prevDelVar);
                        }
                        
                        break;
                    case SV:
                        //introduce variant
                        //no insertion -> treat like simple deletion
                        LoadUInt(ins_sv_arr+vt[currVt].alt+sizeof(uint32_t), insLen);
                        //A_memcpy(&insLen, ins_sv_arr+vt[currVt].alt+sizeof(uint32_t), sizeof(uint32_t));
                        
                        if(insLen  == 0)
                        {
                            //  if(chCovered > 0)
                            {
                                LoadUInt(ins_sv_arr+vt[currVt].alt, delLen);
                                //A_memcpy(&delLen , ins_sv_arr+vt[currVt].alt, sizeof(uint32_t));
                                
                                
                                pos = checkVTconfigurationExact(prePos, noVt, firstVar, currVt+1,  currVt+1, firstVarRead1, refPos+delLen , chCovered, read1, readLen, delta,preChToCheck, &result, pos,max, prevDelCount, firstSV, lastSV, pos_vtToCheck, refPos+delLen, gap, remainedToInsert, prevDelVar);
                                
                            }
                        }
                        else    // deletion+insertion
                        {
                            
                           // int32_t ch;
                            goodPath=true;
                            for(ch = 0; ch < insLen && chCovered+ch < patternLen; ch++)
                            {
                                if(read1[chCovered+ch+delta] != ins_sv_arr[vt[currVt].alt+(sizeof(uint32_t)<<1)+ch])
                                {
                                    goodPath=false;
                                    break;
                                }
                            }
                            if(goodPath)
                            {
                                LoadUInt(ins_sv_arr+vt[currVt].alt, delLen);
                                //A_memcpy(&delLen, ins_sv_arr+vt[currVt].alt, sizeof(uint32_t));
                                
                                if(insLen-ch==0)
                                    pos = checkVTconfigurationExact(prePos,  noVt, firstVar, currVt+1,  currVt+1, firstVarRead1,  refPos+delLen, chCovered+ch,read1, readLen,delta, preChToCheck, &result, pos, max, prevDelCount, firstSV, lastSV, pos_vtToCheck,  refPos+delLen, gap, insLen-ch,  prevDelVar);
                                else
                                    pos = checkVTconfigurationExact(prePos,  noVt, firstVar, currVt+1,  currVt, firstVarRead1,  refPos+delLen, chCovered+ch,read1, readLen,delta, preChToCheck, &result, pos, max, prevDelCount, firstSV, lastSV, pos_vtToCheck,  refPos+delLen, gap, insLen-ch,prevDelVar);
                                
                            }
                        }
                        
                        break;
                }
                
                // free(newVtList);
                // last = vtList[noVt/NO_BITS_UINT];
                
                
                //don't introduce variant
                //   newVtList[noVt/NO_BITS_UINT] = last;
                currVt++;
                lastVt++;
                //   ++noVt;
                vtList[(noVt-1)/NO_BITS_UINT] = last;
                
            }
            else
            {
               if(refPos >= refSize)
               {
                   currPathValid = false;
                   break;
               }
            len = vt[currVt].prePos - refPos + 1 < patternLen - chCovered ? vt[currVt].prePos - refPos + 1 : patternLen - chCovered ;
                // while (refPos <= vt[currVt].prePos && chCovered < patternLen)
                {
                    //if (refPos != data[0].size)
                    {
                        for(int32_t i = 0; i < len; ++i)
                            if(read1[chCovered+delta+i] != getRefFinal(refPos+i))
                            {
                                currPathValid = false;
                                break;
                            }
                        
                        if(!currPathValid)
                            break;
                        
                        chCovered += len;
                        refPos += len;
                        nextRefPos += len;
                        
                        
                        
                    }
                    
                }
                
            }
        }
    }
    int32_t d = 0;
    
    
    
    //delta region...
    if(preChToCheck && currPathValid)// && )
    {
        //if(prePos > 0)
        {
                //Check beginning!!
            currVt = firstVar > 0 ? firstVar - 1 : noVar;
            while (vtType[currVt] == SV || vtType[currVt] == DEL)
                currVt--;
            
            
            
            
            
            if(vt[currVt].prePos != prePos)
            {
                refPos = prePos;
                prePos--;
                d = 0;
            }
            else{
                refPos = prePos+1;
                d = 1;
            }
            
            
                while(d + preChToCheck && prePos  != 0xFFFFFFFE)
                {
                    
                    while(delFromRight[prevDelVar].nextRefPos > prePos + 2 && d + preChToCheck)
                    {
                        //delRegVariatns = true;
                        
                        
                        
                        newRefPos = vt[delFromRight[prevDelVar].vtNo].prePos < refPos ? vt[delFromRight[prevDelVar].vtNo].prePos : refPos;
                        
                        if(newRefPos < refPos)
                        {
                            if(vtType[delFromRight[prevDelVar].vtNo] == DEL)
                            {
                                insLen = 0;
                            }
                            else if(vtType[delFromRight[prevDelVar].vtNo] == SV)
                            {
                                LoadUInt(ins_sv_arr+vt[delFromRight[prevDelVar].vtNo].alt+sizeof(uint32_t), insLen);
                                //A_memcpy(&insLen, ins_sv_arr+vt[delFromRight[prevDelVar].vtNo].alt+sizeof(uint32_t), sizeof(uint32_t));
                            }
                            else
                            {
                                //delLen = 0;
                                insLen = 0;
                            }
                            
                            
                            ///
                            (vtToCheck)[pos_vtToCheck++] = delFromRight[prevDelVar].vtNo;
                            if(pos_vtToCheck == max_vtToCheck)
                            {
                                max_vtToCheck = max_vtToCheck * 2;
                                
                                temp_vtToCheck = ( uint32_t *) realloc (vtToCheck, max_vtToCheck * sizeof( uint32_t));
                                if (NULL == temp_vtToCheck)
                                {
                                    printf("SA: Error in realloc, the message is - %s\n max = %d\n", strerror(errno), max_vtToCheck);
                                    fflush(stdout);
                                    //free(result);
                                    exit(8);
                                }
                                else
                                {
                                    vtToCheck = temp_vtToCheck;
                                }
                            }
                            ///
                            
                            last = vtList[noVt/NO_BITS_UINT];
                            vtList[noVt/NO_BITS_UINT] = vtList[noVt/NO_BITS_UINT] | (1U << (NO_BITS_UINT -  (noVt&NO_BITS_UINT_MINUS_ONE) -  1U));
                            ++noVt;
                            if(noVt == vtListSize * NO_BITS_UINT)
                            {
                                vtListSize++;
                                uint32_t * tempList = (uint32_t *) realloc(vtList, vtListSize * sizeof(uint32_t));
                                if (NULL ==tempList) {
                                    printf("Error in realloc, the message is - %s\n", strerror(errno));
                                    exit(8);
                                }
                                else
                                {
                                    vtList = tempList;
                                }
                                vtList[vtListSize-1]=0;
                            }
                            
                            
                            
                            if(insLen == 0)
                            {
                                pos = checkVTconfigurationExact(newRefPos,  noVt, delFromRight[prevDelVar].vtNo, currVt, lastVt,firstVarRead1,  newRefPos, chCovered, read1, readLen, delta, preChToCheck,  &result, pos, max, prevDelCount, firstSV, lastSV, pos_vtToCheck, nextRefPos, gap, remainedToInsert, prevDelVar > 0 ? prevDelVar-1 : 0);
                            }
                            else
                            {
                                
                                ///
                                int32_t c=0;
                                goodPath=true;
                                for(ch = insLen - 1; ch >= 0  && c+1 <= preChToCheck+d; --ch)
                                {
                                    if(read1[preChToCheck-c-1+d] != ins_sv_arr[vt[delFromRight[prevDelVar].vtNo].alt+(sizeof(uint32_t)<<1)+ch])
                                    {
                                        goodPath=false;
                                        break;
                                    }
                                    ++c;
                                }
                                
                                if(goodPath)
                                {
                                    pos = checkVTconfigurationExact(newRefPos, noVt, delFromRight[prevDelVar].vtNo, currVt,  lastVt, firstVarRead1, newRefPos, chCovered+c, read1, readLen, delta, preChToCheck-c+d, &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,   nextRefPos, ch+1, remainedToInsert,  prevDelVar > 0 ? prevDelVar - 1 : 0);
                                    
                                }
                                
                                
                                
                                
                            }
                            
                            vtList[(noVt-1)/NO_BITS_UINT] = last;
                            
                        }
                        
                        
                        
                       prevDelVar = prevDelVar > 0 ? prevDelVar - 1 : 0; // --prevDelVar;
                        
                        
                    }
                    
                    if(delFromRight[prevDelVar].nextRefPos == prePos + 2)
                    {
                        if(vtType[delFromRight[prevDelVar].vtNo] == DEL)
                        {
                            insLen = 0;
                        }
                        else if(vtType[delFromRight[prevDelVar].vtNo] == SV)
                        {
                            LoadUInt(ins_sv_arr+vt[delFromRight[prevDelVar].vtNo].alt+sizeof(uint32_t), insLen);
                           
                        }
                        else
                        {
                            insLen = 0;
                        }
                        
                        
                        (vtToCheck)[(pos_vtToCheck)++] = delFromRight[prevDelVar].vtNo;
                        if(pos_vtToCheck == max_vtToCheck)
                        {
                            max_vtToCheck = max_vtToCheck * 2;
                            
                            temp_vtToCheck = ( uint32_t *) realloc (vtToCheck, max_vtToCheck * sizeof( uint32_t));
                            if (NULL == temp_vtToCheck)
                            {
                                printf("SA: Error in realloc, the message is - %s\n max = %d\n", strerror(errno), max_vtToCheck);
                                fflush(stdout);
                                //free(result);
                                exit(8);
                            }
                            else
                            {
                                vtToCheck = temp_vtToCheck;
                            }
                        }
                        last = vtList[noVt/NO_BITS_UINT];
                        vtList[noVt/NO_BITS_UINT] = vtList[noVt/NO_BITS_UINT] | (1U << (NO_BITS_UINT -  (noVt&NO_BITS_UINT_MINUS_ONE) -  1U));
                        ++noVt;
                        if(noVt == vtListSize * NO_BITS_UINT)
                        {
                            vtListSize++;
                            uint32_t * tempList = (uint32_t *) realloc(vtList, vtListSize * sizeof(uint32_t));
                            if (NULL ==tempList) {
                                printf("Error in realloc, the message is - %s\n", strerror(errno));
                                exit(8);
                            }
                            else
                            {
                                vtList = tempList;
                            }
                            vtList[vtListSize-1]=0;
                        }
                        
                        
                        if(vtType[delFromRight[prevDelVar].vtNo] == DEL)
                        {
                            
                            
                            pos = checkVTconfigurationExact(vt[delFromRight[prevDelVar].vtNo].prePos, noVt, delFromRight[prevDelVar].vtNo, currVt, lastVt, firstVarRead1, vt[delFromRight[prevDelVar].vtNo].prePos, chCovered, read1, readLen,delta, preChToCheck, &result, pos, max, prevDelCount, firstSV, lastSV, pos_vtToCheck,   nextRefPos, gap, remainedToInsert, prevDelVar > 0 ? prevDelVar - 1 : 0);
                            
                            
                        }
                        else //SV
                        {
                            if(insLen  == 0)
                            {
                                
                                
                                pos = checkVTconfigurationExact(vt[delFromRight[prevDelVar].vtNo].prePos, noVt, delFromRight[prevDelVar].vtNo, currVt, lastVt,firstVarRead1,  vt[delFromRight[prevDelVar].vtNo].prePos, chCovered, read1, readLen,delta, preChToCheck, &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  nextRefPos, gap, remainedToInsert, prevDelVar > 0 ? prevDelVar - 1 : 0);
                                
                            }
                            else    // deletion+insertion
                            {
                                int32_t c = 0;
                                goodPath=true;
                                for(ch = insLen -1 ; ch >=0  && c+1 <= preChToCheck+d; --ch)
                                {
                                    if(read1[preChToCheck-c-1+d] != ins_sv_arr[vt[delFromRight[prevDelVar].vtNo].alt+(sizeof(uint32_t)<<1)+ch])
                                    {
                                        goodPath=false;
                                        break;
                                    }
                                    ++c;
                                }
                                if(goodPath)
                                {
                                    
                                    pos = checkVTconfigurationExact(vt[delFromRight[prevDelVar].vtNo].prePos, noVt, delFromRight[prevDelVar].vtNo, currVt,  currVt, firstVarRead1, vt[delFromRight[prevDelVar].vtNo].prePos, chCovered+c, read1, readLen, delta, preChToCheck-c+d, &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,   nextRefPos, ch+1, remainedToInsert,  prevDelVar > 0 ? prevDelVar - 1: 0);
                                  
                                }
                            }
                            
                        }
                        
                        vtList[(noVt-1)/NO_BITS_UINT] = last;
                        prevDelVar = prevDelVar > 0 ? prevDelVar - 1 : 0;//prevDelVar--;
                    }
                    
                    
                    
                    //prePos--;
                    else if (vt[currVt].prePos == prePos)
                    {
                        (vtToCheck)[(pos_vtToCheck)++] = currVt;
                        if(pos_vtToCheck == max_vtToCheck)
                        {
                            max_vtToCheck = max_vtToCheck * 2;
                            
                            temp_vtToCheck = ( uint32_t *) realloc (vtToCheck, max_vtToCheck * sizeof( uint32_t));
                            if (NULL == temp_vtToCheck)
                            {
                                printf("SA: Error in realloc, the message is - %s\n max = %d\n", strerror(errno), max_vtToCheck);
                                fflush(stdout);
                                //free(result);
                                exit(8);
                            }
                            else
                            {
                                vtToCheck = temp_vtToCheck;
                            }
                        }
                        last = vtList[noVt/NO_BITS_UINT];
                        vtList[noVt/NO_BITS_UINT] = vtList[noVt/NO_BITS_UINT] | (1U << (NO_BITS_UINT -  (noVt&NO_BITS_UINT_MINUS_ONE) -  1U));
                        ++noVt;
                        if(noVt == vtListSize * NO_BITS_UINT)
                        {
                            vtListSize++;
                            uint32_t * tempList = (uint32_t *) realloc(vtList, vtListSize * sizeof(uint32_t));
                            if (NULL ==tempList) {
                                printf("Error in realloc, the message is - %s\n", strerror(errno));
                                exit(8);
                            }
                            else
                            {
                                vtList = tempList;
                            }
                            vtList[vtListSize-1]=0;
                        }
                        
                        if(vtType[currVt] == SNP)
                        {
                            if (vt[currVt].alt == read1[preChToCheck+d-1])
                                pos = checkVTconfigurationExact(prePos, noVt, currVt, currVt,  lastVt, firstVarRead1, refPos, chCovered, read1, readLen, delta, preChToCheck-1, &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,   nextRefPos, gap, remainedToInsert, prevDelVar);
                            
                            
                        }
                        else //(vtType[currVt] == INS)
                        {
                            
                            int32_t c = 0;
                            ch=0;
                            goodPath=true;
                            if (d == 0 && getRefFinal(refPos) != read1[preChToCheck-1])  //ref
                                goodPath=false;
                            else
                            {
                                LoadUInt(ins_sv_arr+vt[currVt].alt, insLen);
                                
                                for(ch = insLen-1; ch >= 0  && c+2 <= preChToCheck+d; --ch)
                                {
                                    if(read1[preChToCheck+d-2-c] != ins_sv_arr[vt[currVt].alt+sizeof(uint32_t)+ch])
                                    {
                                        goodPath=false;
                                        break;
                                    }
                                    ++c;
                                }
                            }
                            if(goodPath)
                                pos = checkVTconfigurationExact(prePos, noVt, currVt, currVt, lastVt, firstVarRead1,  refPos, chCovered, read1, readLen, delta, preChToCheck-c-1+d, &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  nextRefPos, ch+1, remainedToInsert, prevDelVar);
                            
                            
                        }
                        
                        firstVar = currVt;
                        vtList[(noVt-1)/NO_BITS_UINT] = last;
                        
                        if(d == 0)
                        {
                            
                            
                            if(getRefFinal(refPos) == read1[preChToCheck-1])
                            {
                                
                                
                                currVt = currVt > 0 ? currVt - 1 : noVar;
                                while (vtType[currVt] == SV || vtType[currVt] == DEL)
                                    currVt--;
                                
                                if(vt[currVt].prePos != prePos)
                                {
                                    refPos = prePos;
                                    prePos--;
                                    d =0;
                                    --preChToCheck;
                                }
                                else
                                {
                                    d = 1;
                                    --preChToCheck;
                                }
                                
                            }
                            else
                            {
                                currPathValid = false;
                                break;
                            }
                        }
                        else
                        {
                            vtList[(noVt-1)/NO_BITS_UINT] = last;
                            
                            currVt = currVt > 0 ? currVt - 1 : noVar;
                            while (vtType[currVt] == SV || vtType[currVt] == DEL)
                                currVt--;
                            
                            if(vt[currVt].prePos != prePos )
                            {
                                refPos = prePos;
                                prePos--;
                                d=0;
                            }
                            else
                            {
                                d =1;
                            }
                            
                        }
                        
                        
                        
                    }
                    else
                    {
                        if (getRefFinal(refPos) == read1[preChToCheck-1])
                        {
                            refPos = prePos;
                            prePos--;
                            d =0;
                            --preChToCheck;
                        }
                        else
                        {
                            currPathValid = false;
                            break;
                        }
                        
                    }
                }
                ++prePos;
        }
        if(preChToCheck > 0)
        {
            currPathValid = false;
        }
    }
    
    
    
    
    
    //deletions at the beginnig?
    if(currPathValid)
    {
        while(delFromRight[prevDelVar].nextRefPos >= prePos + 1 && prePos != 0xFFFFFFFF)  // can't be >, unless del in del region
        {
            (vtToCheck)[(pos_vtToCheck)++] = delFromRight[prevDelVar].vtNo;
            if(pos_vtToCheck == max_vtToCheck)
            {
                max_vtToCheck = max_vtToCheck * 2;
                
                temp_vtToCheck = ( uint32_t *) realloc (vtToCheck, max_vtToCheck * sizeof( uint32_t));
                if (NULL == temp_vtToCheck)
                {
                    printf("SA: Error in realloc, the message is - %s\n max = %d\n", strerror(errno), max_vtToCheck);
                    fflush(stdout);
                    //free(result);
                    exit(8);
                }
                else
                {
                    vtToCheck = temp_vtToCheck;
                }
            }
            last = vtList[noVt/NO_BITS_UINT];
            vtList[noVt/NO_BITS_UINT] = vtList[noVt/NO_BITS_UINT] | (1U << (NO_BITS_UINT -  (noVt&NO_BITS_UINT_MINUS_ONE) -  1U));
            ++noVt;
            if(noVt == vtListSize * NO_BITS_UINT)
            {
                vtListSize++;
                uint32_t * tempList = (uint32_t *) realloc(vtList, vtListSize * sizeof(uint32_t));
                if (NULL ==tempList) {
                    printf("Error in realloc, the message is - %s\n", strerror(errno));
                    exit(8);
                }
                else
                {
                    vtList = tempList;
                }
                vtList[vtListSize-1]=0;
            }
            
            
            if(vtType[delFromRight[prevDelVar].vtNo] == DEL)
            {
                
                
                pos = checkVTconfigurationExact(vt[delFromRight[prevDelVar].vtNo].prePos, noVt, delFromRight[prevDelVar].vtNo, currVt, lastVt, firstVarRead1, vt[delFromRight[prevDelVar].vtNo].prePos, chCovered, read1, readLen,delta, preChToCheck, &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  nextRefPos, gap, remainedToInsert,  prevDelVar > 0 ? prevDelVar - 1 : 0);
                
                
            }
            else //SV
            {
                if(insLen  == 0)
                {
                    
                    
                    pos = checkVTconfigurationExact(vt[delFromRight[prevDelVar].vtNo].prePos, noVt, delFromRight[prevDelVar].vtNo, currVt, lastVt,firstVarRead1,  vt[delFromRight[prevDelVar].vtNo].prePos, chCovered, read1, readLen,delta, preChToCheck, &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  nextRefPos, gap, remainedToInsert, prevDelVar > 0 ? prevDelVar - 1 : 0);
                    
                }
                else    // deletion+insertion
                {
                    
                    int32_t c = 0;
                    ch = 0;
                    goodPath=true;
                    for(ch = insLen -1 ; ch >=0  && c <= preChToCheck-1; --ch)
                    {
                        if(read1[preChToCheck-1-c] != ins_sv_arr[vt[delFromRight[prevDelVar].vtNo].alt+(sizeof(uint32_t)<<1)+ch])
                        {
                            
                            goodPath=false;
                        }
                        ++c;
                    }
                    if(goodPath)
                    {
                        pos = checkVTconfigurationExact(vt[delFromRight[prevDelVar].vtNo].prePos, noVt, delFromRight[prevDelVar].vtNo, currVt,  lastVt, firstVarRead1, vt[delFromRight[prevDelVar].vtNo].prePos, chCovered+c, read1, readLen, delta, preChToCheck, &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  nextRefPos, ch+1, remainedToInsert, prevDelVar > 0 ? prevDelVar - 1 : 0);
                        
                    }
                }
                
            }
            
            vtList[(noVt-1)/NO_BITS_UINT] = last;
            prevDelVar = prevDelVar > 0 ? prevDelVar - 1 : 0; //prevDelVar--;
        }
    }
    
    
    
    
    
    
    
    
    if(currPathValid)
    {
        
        
        bool suffix, vtFound;
        
        
        
        uint64_t  *tempResult;
        
        int32_t back, dels;
        
        
        if(pos_vtToCheck == 0 && prevDelCount == 0 && firstSV>lastSV)
        {
            //tempPos = inBitV;
            //flag
            result[pos++] = 0xFFFFFFFFFFFFFFFE;
            
            if (pos  == *max)
            {
                *max = (*max) * 2;
                
                
                tempResult = ( uint64_t  *) realloc (result, (*max) * sizeof( uint64_t ));
                if (NULL == tempResult)
                {
                    printf("SA: Error in realloc, the message is - %s\n max = %llu\n", strerror(errno), *max);
                    fflush(stdout);
                    //free(result);
                    exit(8);
                }
                else
                {
                    result = tempResult;
                }
            }
            
            // pos
            result[pos++] = prePos;
            if (pos == *max)
            {
                *max = (*max) * 2 ;
                
                
                tempResult = ( uint64_t  *) realloc (result, (*max) * sizeof( uint64_t ));
                if (NULL == tempResult)
                {
                    printf("SA: Error in realloc, the message is - %s\n max = %llu\n", strerror(errno), *max);
                    fflush(stdout);
                    //free(result);
                    exit(8);
                }
                else
                {
                    result = tempResult;
                }
            }
            
    
        }
        else
        {
            
            uint32_t aux;
            
            A_memset(resV, 0xFF, resSize<<3);
            
            
            dels = prevDelCount;
            back = 1;
            while(dels)
            {
                if(vtType[firstVarRead1 - back] == DEL)
                {
                    dels--;
                    if(vt[firstVarRead1 - back].prePos + vt[firstVarRead1 - back].alt > prePos)
                    {
                        aux = all_var[(firstVarRead1 - back)]*parts;
                        for(int32_t r=0; r < resSize; r++)
                            resV[r] = resV[r] & ~uniq_parts[uniq_var[aux+r/blocking]*blocking+(r)%blocking];
                        
                    }
                }
                back++;
            }
            
            if(firstSV < noSV)
                for(int32_t sv = lastSV; sv >= (int32_t)firstSV; sv--)
                {
                    LoadUInt(ins_sv_arr+vt[listSV[sv]].alt, delLen);
                    
                    if(vt[listSV[sv]].prePos +  delLen > prePos)
                    {
                        aux = all_var[listSV[sv]]*parts;
                        for(int32_t r=0; r < resSize; r++)
                            resV[r] = resV[r] & ~uniq_parts[uniq_var[aux+r/blocking]*blocking+(r)%blocking];
                    }
                }
            
            for(int32_t variant=0; variant < pos_vtToCheck; variant++)
            {
                aux = all_var[(vtToCheck)[variant]]*parts;
                suffix = vtList[variant/NO_BITS_UINT] & (1 << (NO_BITS_UINT - variant%NO_BITS_UINT - 1));
                if(suffix)
                    for(int32_t r=0; r < resSize; r++)
                        resV[r] = resV[r] & uniq_parts[uniq_var[aux+r/blocking]*blocking+(r)%blocking];
                else
                    for(int32_t r=0; r < resSize; r++)
                        resV[r] = resV[r] & ~uniq_parts[uniq_var[aux+r/blocking]*blocking+(r)%blocking];
                
            }
            
            vtFound = false;
            for ( int32_t r = 0; r < resSize; r++)
            {
                if(resV[r] > 0)
                {
                    vtFound = true;
                    break;
                }
            }
            
            if(vtFound)
            {
                
                result[pos++] = 0xFFFFFFFFFFFFFFFF;
                if (pos  == *max)
                {
                    *max = (*max) * 2;
                    
                    
                    tempResult = ( uint64_t  *) realloc (result, (*max) * sizeof( uint64_t ));
                    if (NULL == tempResult)
                    {
                        printf("SA: Error in realloc, the message is - %s\n max = %llu\n", strerror(errno), *max);
                        fflush(stdout);
                        exit(8);
                    }
                    else
                    {
                        result = tempResult;
                    }
                }
                
                // pos
                result[pos++] = prePos;
                if (pos  == *max)
                {
                    *max = (*max) * 2;
                    
                    
                    tempResult = ( uint64_t  *) realloc (result, (*max) * sizeof( uint64_t ));
                    if (NULL == tempResult)
                    {
                        printf("SA: Error in realloc, the message is - %s\n max = %llu\n", strerror(errno), *max);
                        fflush(stdout);
                        //free(result);
                        exit(8);
                    }
                    else
                    {
                        result = tempResult;
                    }
                }
                
                
                if (pos + resSize >= *max)
                {
                    *max = (*max) * 2 > (*max) + resSize ? (*max) * 2 : (*max) + 2*resSize ;
                    
                    
                    tempResult = ( uint64_t  *) realloc (result, (*max) * sizeof( uint64_t ));
                    if (NULL == tempResult)
                    {
                        printf("SA: Error in realloc, the message is - %s\n max = %llu\n", strerror(errno), *max);
                        fflush(stdout);
                        //free(result);
                        exit(8);
                    }
                    else
                    {
                        result = tempResult;
                    }
                }
                
                A_memcpy(result+pos, resV, resSize<<3); pos = pos + resSize;
                
            }
            
        }

        
        
    }
    
    //  free(vtList);
    *passed_result=result;
    return pos;
}


uint64_t checkVTconfiguration1(uint32_t prePos,  uint32_t noVt,  uint32_t firstVar, uint32_t currVt,  uint32_t lastVt,  uint32_t firstVarRead1, uint32_t refPos, uint32_t chCovered, unsigned char * read1, int32_t readLen, uint32_t delta,  int32_t preChToCheck,  uint64_t ** passed_result, uint64_t pos, uint64_t * max, uint32_t prevDelCount, uint32_t firstSV, uint32_t lastSV,  uint32_t  pos_vtToCheck,  unsigned char * wholeRead, uint32_t wholeReadLen,   int32_t leftPos, uint32_t rightPos, uint32_t error, uint32_t maxError,  uint32_t nextRefPos, uint32_t gap, uint32_t remainedToInsert, uint32_t startPosFoundExact, uint32_t prePosExact, uint32_t prevDelVar)

{
  
    
    
    uint64_t * result = *passed_result;
    
    
    int32_t tempError;
    int32_t ch;
    
    
    uint32_t tempLen;
    
    uint32_t * temp_vtToCheck;
    bool goodPath, currPathValid = true;;
    uint32_t last = 0, len = 0, delLen = 0, insLen = 0;
    uint32_t newRefPos;
    
    uint32_t patternLen = readLen-delta;
    
    
    //
    if(remainedToInsert && chCovered < patternLen)
    {
        if(vtType[lastVt] == INS)
        {
            LoadUInt(ins_sv_arr+vt[lastVt].alt, insLen);
            tempLen = vt[lastVt].alt+sizeof(uint32_t);
        }
        else if(vtType[lastVt] == SV)//SV
        {
            LoadUInt(ins_sv_arr+vt[lastVt].alt+sizeof(uint32_t), insLen);
            tempLen = vt[lastVt].alt+(sizeof(uint32_t)<<1);
        }
        else
        {
            cout << "Should be SV or INS! " <<  lastVt << " " << prePos <<" " << remainedToInsert << "\n";
            exit(9);
        }
        
        
        
        while (chCovered < patternLen && remainedToInsert > 0)
        {
            if ((ins_sv_arr[tempLen+insLen-remainedToInsert]) == wholeRead[leftPos+delta+chCovered+1])
            {
                chCovered++; remainedToInsert--;
            }
            else
            {
                currPathValid = false;
                break;
            }
        }
        if(!remainedToInsert)
        {
            lastVt++;
        }
    }
    
    while (currPathValid && chCovered < patternLen)
    {
        
        while(vt[currVt].prePos+1 < refPos && chCovered < patternLen) //currVt < noVar &&
        {
            //delRegVariatns = true;
            
            if(vtType[currVt] == DEL || vtType[currVt] == SV)
            {
                if(vtType[currVt] == DEL)
                {
                    delLen = vt[currVt].alt;
                    insLen = 0;
                }
                else
                {
                    LoadUInt(ins_sv_arr+vt[currVt].alt, delLen);
                    LoadUInt(ins_sv_arr+vt[currVt].alt+sizeof(uint32_t), insLen);
                }
                
                
                newRefPos = vt[currVt].prePos+1+delLen > refPos ? vt[currVt].prePos+1+delLen : refPos;
                
                if(newRefPos > refPos)
                {
                    ///
                    (vtToCheck)[pos_vtToCheck++] = currVt;
                    if(pos_vtToCheck == max_vtToCheck)
                    {
                        max_vtToCheck = max_vtToCheck * 2;
                        
                        temp_vtToCheck = ( uint32_t *) realloc (vtToCheck, max_vtToCheck * sizeof( uint32_t));
                        if (NULL == temp_vtToCheck)
                        {
                            printf("SA: Error in realloc, the message is - %s\n max = %d\n", strerror(errno), max_vtToCheck);
                            fflush(stdout);
                            //free(result);
                            exit(8);
                        }
                        else
                        {
                            vtToCheck = temp_vtToCheck;
                        }
                    }
                    ///
                    
                    
                    last = vtList[noVt/NO_BITS_UINT];
                    vtList[noVt/NO_BITS_UINT] = vtList[noVt/NO_BITS_UINT] | (1U << (NO_BITS_UINT -  (noVt&NO_BITS_UINT_MINUS_ONE) -  1U));
                    ++noVt;
                    if(noVt == vtListSize * NO_BITS_UINT)
                    {
                        vtListSize++;
                        uint32_t * tempList = (uint32_t *) realloc(vtList, vtListSize * sizeof(uint32_t));
                        if (NULL ==tempList) {
                            printf("Error in realloc, the message is - %s\n", strerror(errno));
                            exit(8);
                        }
                        else
                        {
                            vtList = tempList;
                        }
                        vtList[vtListSize-1]=0;
                    }
                    
                    
                    
                    if(insLen == 0)
                    {
                        pos = checkVTconfiguration1(prePos, noVt, firstVar, currVt+1, currVt+1, firstVarRead1, newRefPos, chCovered, read1, readLen, delta, preChToCheck,  &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  wholeRead, wholeReadLen, leftPos, rightPos, error, maxError, newRefPos, gap, remainedToInsert, startPosFoundExact, prePosExact, prevDelVar);
                    }
                    else
                    {
                        
                        ///
                        //int32_t ch;
                        goodPath=true;
                        for(ch = 0; ch < insLen && chCovered+ch < patternLen; ch++)
                        {
                            if(read1[chCovered+ch+delta] != ins_sv_arr[vt[currVt].alt+(sizeof(uint32_t)<<1)+ch])
                            {
                                goodPath=false;
                                break;
                            }
                        }
                        if(goodPath)
                        {
                            if(insLen-ch == 0)
                                pos = checkVTconfiguration1(prePos, noVt, firstVar, currVt+1,  currVt+1, firstVarRead1, newRefPos, chCovered+ch, read1, readLen, delta, preChToCheck, &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  wholeRead, wholeReadLen, leftPos, rightPos, error, maxError, newRefPos, gap, insLen-ch, startPosFoundExact, prePosExact, prevDelVar);
                            else
                                pos = checkVTconfiguration1(prePos, noVt, firstVar, currVt+1,  currVt, firstVarRead1, newRefPos, chCovered+ch, read1, readLen, delta, preChToCheck, &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  wholeRead, wholeReadLen, leftPos, rightPos, error, maxError, newRefPos, gap, insLen-ch, startPosFoundExact, prePosExact, prevDelVar);
                        }
                        
                        
                        ///
                    }
                    
                    vtList[(noVt-1)/NO_BITS_UINT] = last;
                    
                }
                
                
                
            }
            
            
            currVt++;
            lastVt++;
            
        }
        
        if( chCovered < patternLen)
        {
            if( vt[currVt].prePos+1 == refPos && refPos < refSize) //currVt < noVar &&
            {
                (vtToCheck)[(pos_vtToCheck)++] = currVt;
                if(pos_vtToCheck == max_vtToCheck)
                {
                    max_vtToCheck = max_vtToCheck * 2;
                    
                    temp_vtToCheck = ( uint32_t *) realloc (vtToCheck, max_vtToCheck * sizeof( uint32_t));
                    if (NULL == temp_vtToCheck)
                    {
                        printf("SA: Error in realloc, the message is - %s\n max = %d\n", strerror(errno), max_vtToCheck);
                        fflush(stdout);
                        //free(result);
                        exit(8);
                    }
                    else
                    {
                        vtToCheck = temp_vtToCheck;
                    }
                }
                ///
                
                
                
                last = vtList[noVt/NO_BITS_UINT];
                vtList[noVt/NO_BITS_UINT] = vtList[noVt/NO_BITS_UINT] | (1U << (NO_BITS_UINT -  (noVt&NO_BITS_UINT_MINUS_ONE) -  1U));
                ++noVt;
                
                if(noVt == vtListSize * NO_BITS_UINT)
                {
                    vtListSize++;
                    uint32_t * tempList = (uint32_t *) realloc(vtList, vtListSize * sizeof(uint32_t));
                    if (NULL ==tempList) {
                        printf("Error in realloc, the message is - %s\n", strerror(errno));
                        exit(8);
                    }
                    else
                    {
                        vtList = tempList;
                    }
                    vtList[vtListSize-1]=0;
                }
                
                
                
                switch (vtType[currVt])
                {
                    case SNP:
                        //introduce variant
                        
                        if(read1[chCovered+delta] == vt[currVt].alt)
                            pos = checkVTconfiguration1(prePos, noVt, firstVar, currVt+1,  currVt+1, firstVarRead1, refPos+1, chCovered+1, read1, readLen,delta, preChToCheck, &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  wholeRead, wholeReadLen, leftPos, rightPos, error, maxError, refPos+1, gap, remainedToInsert, startPosFoundExact, prePosExact, prevDelVar);
                        
                        break;
                    case DEL:
                        //introduce variant
                        
                        
                        
                        pos = checkVTconfiguration1(prePos, noVt, firstVar, currVt+1,  currVt+1, firstVarRead1, refPos+vt[currVt].alt, chCovered, read1, readLen,delta, preChToCheck, &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  wholeRead, wholeReadLen, leftPos, rightPos, error, maxError, refPos+vt[currVt].alt, gap, remainedToInsert, startPosFoundExact, prePosExact, prevDelVar);
                        
                        
                        break;
                    case INS:
                        //introduce variant
                        
                        //int32_t ch;
                        goodPath=true;
                        LoadUInt(ins_sv_arr+vt[currVt].alt, insLen);
                        //A_memcpy(&insLen, ins_sv_arr+vt[currVt].alt, sizeof(uint32_t));
                        for(ch = 0; ch <  insLen && chCovered+ch < patternLen; ch++)
                        {
                            if(read1[chCovered+ch+delta] != ins_sv_arr[vt[currVt].alt+sizeof(uint32_t)+ch])
                            {
                                goodPath=false;
                                break;
                            }
                        }
                        if(goodPath)
                        {  if(insLen-ch == 0)
                            pos = checkVTconfiguration1(prePos, noVt, firstVar, currVt+1, currVt+1,  firstVarRead1, refPos, chCovered+ch, read1, readLen, delta, preChToCheck, &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  wholeRead, wholeReadLen, leftPos, rightPos, error, maxError, refPos, gap, insLen-ch, startPosFoundExact, prePosExact, prevDelVar);
                        else
                            pos = checkVTconfiguration1(prePos, noVt, firstVar, currVt+1, currVt, firstVarRead1,  refPos, chCovered+ch, read1, readLen, delta, preChToCheck, &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  wholeRead, wholeReadLen, leftPos, rightPos, error, maxError, refPos, gap, insLen-ch, startPosFoundExact, prePosExact, prevDelVar);
                        }
                        
                        break;
                    case SV:
                        //introduce variant
                        //no insertion -> treat like simple deletion
                        LoadUInt(ins_sv_arr+vt[currVt].alt+sizeof(uint32_t), insLen);
                        //A_memcpy(&insLen, ins_sv_arr+vt[currVt].alt+sizeof(uint32_t), sizeof(uint32_t));
                        
                        if(insLen  == 0)
                        {
                            //  if(chCovered > 0)
                            {
                                LoadUInt(ins_sv_arr+vt[currVt].alt, delLen);
                                //A_memcpy(&delLen , ins_sv_arr+vt[currVt].alt, sizeof(uint32_t));
                                
                                
                                pos = checkVTconfiguration1(prePos, noVt, firstVar, currVt+1,  currVt+1, firstVarRead1, refPos+delLen , chCovered, read1, readLen, delta,preChToCheck, &result, pos,max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  wholeRead, wholeReadLen, leftPos, rightPos, error, maxError, refPos+delLen, gap, remainedToInsert, startPosFoundExact, prePosExact, prevDelVar);
                                
                            }
                        }
                        else    // deletion+insertion
                        {
                            
                            //int32_t ch;
                            goodPath=true;
                            for(ch = 0; ch < insLen && chCovered+ch < patternLen; ch++)
                            {
                                if(read1[chCovered+ch+delta] != ins_sv_arr[vt[currVt].alt+(sizeof(uint32_t)<<1)+ch])
                                {
                                    goodPath=false;
                                    break;
                                }
                            }
                            if(goodPath)
                            {
                                LoadUInt(ins_sv_arr+vt[currVt].alt, delLen);
                                //A_memcpy(&delLen, ins_sv_arr+vt[currVt].alt, sizeof(uint32_t));
                                
                                if(insLen-ch==0)
                                    pos = checkVTconfiguration1(prePos, noVt, firstVar, currVt+1,  currVt+1, firstVarRead1,  refPos+delLen, chCovered+ch,read1, readLen,delta, preChToCheck, &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  wholeRead, wholeReadLen, leftPos, rightPos, error, maxError, refPos+delLen, gap, insLen-ch, startPosFoundExact, prePosExact, prevDelVar);
                                else
                                    pos = checkVTconfiguration1(prePos, noVt, firstVar, currVt+1,  currVt, firstVarRead1,  refPos+delLen, chCovered+ch,read1, readLen,delta, preChToCheck, &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  wholeRead, wholeReadLen, leftPos, rightPos, error, maxError, refPos+delLen, gap, insLen-ch, startPosFoundExact, prePosExact, prevDelVar);
                                
                            }
                        }
                        
                        break;
                }
                
                currVt++;
                lastVt++;
                vtList[(noVt-1)/NO_BITS_UINT] = last;
                
            }
            else
            {
                if(refPos >= refSize)
                {
                    currPathValid = false;
                    break;
                }
                len = vt[currVt].prePos - refPos + 1 < patternLen - chCovered ? vt[currVt].prePos - refPos + 1 : patternLen - chCovered ;
                {
                    //if (refPos != data[0].size)
                    {
                        for(int32_t i = 0; i < len; ++i)
                            if(read1[chCovered+delta+i] != getRefFinal(refPos+i))
                            {
                                currPathValid = false;
                                break;
                            }
                        
                        if(!currPathValid)
                            break;
                        
                        chCovered += len;
                        refPos += len;
                        nextRefPos += len;
                        
                        
                        
                    }
                    
                }
                
            }
        }
    }
    int32_t d = 0;
    
    
    
    //delta region...
    if(preChToCheck && currPathValid)
    {
        
       // if(prePos > 0)
        {
                //Check beginning!!
            currVt = firstVar > 0 ? firstVar - 1 : noVar;
            while (vtType[currVt] == SV || vtType[currVt] == DEL)
                currVt--;
            
            
            
            
            
            if(vt[currVt].prePos != prePos)
            {
                refPos = prePos;
                prePos--;
                d = 0;
            }
            else{
                refPos = prePos+1;
                d = 1;
            }
            
            //if(prePos > 0)
            {
                while(d + preChToCheck && prePos != 0xFFFFFFFE)
                {
                    
                    while(delFromRight[prevDelVar].nextRefPos > prePos + 2 && d + preChToCheck)
                    {
                        
                        
                        
                        newRefPos = vt[delFromRight[prevDelVar].vtNo].prePos < refPos ? vt[delFromRight[prevDelVar].vtNo].prePos : refPos;
                        
                        if(newRefPos < refPos)
                        {
                            if(vtType[delFromRight[prevDelVar].vtNo] == DEL)
                            {
                                insLen = 0;
                            }
                            else if(vtType[delFromRight[prevDelVar].vtNo] == SV)
                            {
                                LoadUInt(ins_sv_arr+vt[delFromRight[prevDelVar].vtNo].alt+sizeof(uint32_t), insLen);
                            }
                            else
                            {
                                insLen = 0;
                            }
                            
                            
                            ///
                            (vtToCheck)[pos_vtToCheck++] = delFromRight[prevDelVar].vtNo;
                            if(pos_vtToCheck == max_vtToCheck)
                            {
                                max_vtToCheck = max_vtToCheck * 2;
                                
                                temp_vtToCheck = ( uint32_t *) realloc (vtToCheck, max_vtToCheck * sizeof( uint32_t));
                                if (NULL == temp_vtToCheck)
                                {
                                    printf("SA: Error in realloc, the message is - %s\n max = %d\n", strerror(errno), max_vtToCheck);
                                    fflush(stdout);
                                    //free(result);
                                    exit(8);
                                }
                                else
                                {
                                    vtToCheck = temp_vtToCheck;
                                }
                            }
                            ///
                            
                            last = vtList[noVt/NO_BITS_UINT];
                            vtList[noVt/NO_BITS_UINT] = vtList[noVt/NO_BITS_UINT] | (1U << (NO_BITS_UINT -  (noVt&NO_BITS_UINT_MINUS_ONE) -  1U));
                            ++noVt;
                            if(noVt == vtListSize * NO_BITS_UINT)
                            {
                                vtListSize++;
                                uint32_t * tempList = (uint32_t *) realloc(vtList, vtListSize * sizeof(uint32_t));
                                if (NULL ==tempList) {
                                    printf("Error in realloc, the message is - %s\n", strerror(errno));
                                    exit(8);
                                }
                                else
                                {
                                    vtList = tempList;
                                }
                                vtList[vtListSize-1]=0;
                            }
                            
                            
                            
                            if(insLen == 0)
                            {
                                pos = checkVTconfiguration1(newRefPos, noVt, delFromRight[prevDelVar].vtNo, currVt, lastVt,firstVarRead1,  newRefPos, chCovered, read1, readLen, delta, preChToCheck,  &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  wholeRead, wholeReadLen, leftPos, rightPos, error, maxError, nextRefPos, gap, remainedToInsert, startPosFoundExact, prePosExact, prevDelVar > 0 ? prevDelVar - 1 : 0);
                            }
                            else
                            {
                                
                                ///
                                int32_t c=0;
                                goodPath=true;
                                for(ch = insLen - 1; ch >= 0  && c+1 <= preChToCheck+d; --ch)
                                {
                                    if(read1[preChToCheck-c-1+d] != ins_sv_arr[vt[delFromRight[prevDelVar].vtNo].alt+(sizeof(uint32_t)<<1)+ch])
                                    {
                                        goodPath=false;
                                        break;
                                    }
                                    ++c;
                                }
                                
                                if(goodPath)
                                {
                                    pos = checkVTconfiguration1(newRefPos, noVt, delFromRight[prevDelVar].vtNo, currVt,  lastVt, firstVarRead1, newRefPos, chCovered+c, read1, readLen, delta, preChToCheck-c+d, &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  wholeRead, wholeReadLen, leftPos, rightPos, error, maxError, nextRefPos, ch+1, remainedToInsert, startPosFoundExact, prePosExact, prevDelVar > 0 ? prevDelVar - 1 : 0);
                              
                                }
                                
                                
                                ///
                                
                            }
                            
                            vtList[(noVt-1)/NO_BITS_UINT] = last;
                            
                        }
                        
                        
                        
                        prevDelVar = prevDelVar > 0 ? prevDelVar - 1 : 0; // --prevDelVar;
                        
                        
                    }
                    
                    if(delFromRight[prevDelVar].nextRefPos == prePos + 2)
                    {
                        if(vtType[delFromRight[prevDelVar].vtNo] == DEL)
                        {
                            insLen = 0;
                        }
                        else if(vtType[delFromRight[prevDelVar].vtNo] == SV)
                        {
                            LoadUInt(ins_sv_arr+vt[delFromRight[prevDelVar].vtNo].alt+sizeof(uint32_t), insLen);
                        }
                        else
                        {
                            insLen = 0;
                        }
                        
                        
                        (vtToCheck)[(pos_vtToCheck)++] = delFromRight[prevDelVar].vtNo;
                        if(pos_vtToCheck == max_vtToCheck)
                        {
                            max_vtToCheck = max_vtToCheck * 2;
                            
                            temp_vtToCheck = ( uint32_t *) realloc (vtToCheck, max_vtToCheck * sizeof( uint32_t));
                            if (NULL == temp_vtToCheck)
                            {
                                printf("SA: Error in realloc, the message is - %s\n max = %d\n", strerror(errno), max_vtToCheck);
                                fflush(stdout);
                                //free(result);
                                exit(8);
                            }
                            else
                            {
                                vtToCheck = temp_vtToCheck;
                            }
                        }
                        last = vtList[noVt/NO_BITS_UINT];
                        vtList[noVt/NO_BITS_UINT] = vtList[noVt/NO_BITS_UINT] | (1U << (NO_BITS_UINT -  (noVt&NO_BITS_UINT_MINUS_ONE) -  1U));
                        ++noVt;
                        if(noVt == vtListSize * NO_BITS_UINT)
                        {
                            vtListSize++;
                            uint32_t * tempList = (uint32_t *) realloc(vtList, vtListSize * sizeof(uint32_t));
                            if (NULL ==tempList) {
                                printf("Error in realloc, the message is - %s\n", strerror(errno));
                                exit(8);
                            }
                            else
                            {
                                vtList = tempList;
                            }
                            vtList[vtListSize-1]=0;
                        }
                        
                        
                        if(vtType[delFromRight[prevDelVar].vtNo] == DEL)
                        {
                            
                            
                            pos = checkVTconfiguration1(vt[delFromRight[prevDelVar].vtNo].prePos, noVt, delFromRight[prevDelVar].vtNo, currVt, lastVt, firstVarRead1, vt[delFromRight[prevDelVar].vtNo].prePos, chCovered, read1, readLen,delta, preChToCheck, &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  wholeRead, wholeReadLen, leftPos, rightPos, error, maxError, nextRefPos, gap, remainedToInsert, startPosFoundExact, prePos, prevDelVar > 0 ? prevDelVar - 1 : 0);
                            
                            
                        }
                        else //SV
                        {
                            if(insLen  == 0)
                            {
                                
                                
                                pos = checkVTconfiguration1(vt[delFromRight[prevDelVar].vtNo].prePos, noVt, delFromRight[prevDelVar].vtNo, currVt, lastVt,firstVarRead1,  vt[delFromRight[prevDelVar].vtNo].prePos, chCovered, read1, readLen,delta, preChToCheck, &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  wholeRead, wholeReadLen, leftPos, rightPos, error, maxError, nextRefPos, gap, remainedToInsert, startPosFoundExact, prePos, prevDelVar > 0 ? prevDelVar - 1 : 0);
                                
                            }
                            else    // deletion+insertion
                            {
                                int32_t c = 0;
                                goodPath=true;
                                for(ch = insLen -1 ; ch >=0  && c+1 <= preChToCheck+d; --ch)
                                {
                                    if(read1[preChToCheck-c-1+d] != ins_sv_arr[vt[delFromRight[prevDelVar].vtNo].alt+(sizeof(uint32_t)<<1)+ch])
                                    {
                                        goodPath=false;
                                        break;
                                    }
                                    ++c;
                                }
                                if(goodPath)
                                {
                                    
                                    pos = checkVTconfiguration1(vt[delFromRight[prevDelVar].vtNo].prePos, noVt, delFromRight[prevDelVar].vtNo, currVt,  currVt, firstVarRead1, vt[delFromRight[prevDelVar].vtNo].prePos, chCovered+c, read1, readLen, delta, preChToCheck-c+d, &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  wholeRead, wholeReadLen, leftPos, rightPos, error, maxError, nextRefPos, ch+1, remainedToInsert, startPosFoundExact, prePosExact, prevDelVar > 0 ? prevDelVar - 1 : 0);
                                
                                }
                            }
                            
                        }
                        
                        vtList[(noVt-1)/NO_BITS_UINT] = last;
                       prevDelVar = prevDelVar > 0 ? prevDelVar - 1 : 0;//prevDelVar--;
                    }
                    
                    
                    
                    //prePos--;
                    else if (vt[currVt].prePos == prePos)
                    {
                        (vtToCheck)[(pos_vtToCheck)++] = currVt;
                        if(pos_vtToCheck == max_vtToCheck)
                        {
                            max_vtToCheck = max_vtToCheck * 2;
                            
                            temp_vtToCheck = ( uint32_t *) realloc (vtToCheck, max_vtToCheck * sizeof( uint32_t));
                            if (NULL == temp_vtToCheck)
                            {
                                printf("SA: Error in realloc, the message is - %s\n max = %d\n", strerror(errno), max_vtToCheck);
                                fflush(stdout);
                                //free(result);
                                exit(8);
                            }
                            else
                            {
                                vtToCheck = temp_vtToCheck;
                            }
                        }
                        last = vtList[noVt/NO_BITS_UINT];
                        vtList[noVt/NO_BITS_UINT] = vtList[noVt/NO_BITS_UINT] | (1U << (NO_BITS_UINT -  (noVt&NO_BITS_UINT_MINUS_ONE) -  1U));
                        ++noVt;
                        if(noVt == vtListSize * NO_BITS_UINT)
                        {
                            vtListSize++;
                            uint32_t * tempList = (uint32_t *) realloc(vtList, vtListSize * sizeof(uint32_t));
                            if (NULL ==tempList) {
                                printf("Error in realloc, the message is - %s\n", strerror(errno));
                                exit(8);
                            }
                            else
                            {
                                vtList = tempList;
                            }
                            vtList[vtListSize-1]=0;
                        }
                        
                        if(vtType[currVt] == SNP)
                        {
                            if (vt[currVt].alt == read1[preChToCheck+d-1])
                                pos = checkVTconfiguration1(prePos, noVt, currVt, currVt,  lastVt, firstVarRead1, refPos, chCovered, read1, readLen, delta, preChToCheck-1, &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  wholeRead, wholeReadLen, leftPos, rightPos, error, maxError, nextRefPos, gap, remainedToInsert, startPosFoundExact, prePos, prevDelVar);
                            
                            
                        }
                        else //(vtType[currVt] == INS)
                        {
                            
                            int32_t c = 0;
                            ch=0;
                            goodPath=true;
                            if (d == 0 && getRefFinal(refPos) != read1[preChToCheck-1])  //ref
                                goodPath=false;
                            else
                              
                            {
                                LoadUInt(ins_sv_arr+vt[currVt].alt, insLen);
                               
                                for(ch = insLen-1; ch >= 0  && c+2 <= preChToCheck+d; --ch)
                                {
                                    if(read1[preChToCheck+d-2-c] != ins_sv_arr[vt[currVt].alt+sizeof(uint32_t)+ch])
                                    {
                                        goodPath=false;
                                        break;
                                    }
                                    ++c;
                                }
                            }
                            if(goodPath)
                                pos = checkVTconfiguration1(prePos, noVt, currVt, currVt, lastVt, firstVarRead1,  refPos, chCovered, read1, readLen, delta, preChToCheck-c-1+d, &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  wholeRead, wholeReadLen, leftPos, rightPos, error, maxError, nextRefPos, ch+1, remainedToInsert, startPosFoundExact, prePos, prevDelVar);
                            
                            
                        }
                        
                        firstVar = currVt;
                        vtList[(noVt-1)/NO_BITS_UINT] = last;
                        
                        if(d == 0)
                        {
                            
                            
                            if(getRefFinal(refPos) == read1[preChToCheck-1])
                            {
                                
                                
                                currVt = currVt > 0 ? currVt - 1 : noVar;
                                while (vtType[currVt] == SV || vtType[currVt] == DEL)
                                    currVt--;
                                
                                if(vt[currVt].prePos != prePos)
                                {
                                    refPos = prePos;
                                    prePos--;
                                    d =0;
                                    --preChToCheck;
                                }
                                else
                                {
                                    d = 1;
                                    --preChToCheck;
                                }
                                
                            }
                            else
                            {
                                currPathValid = false;
                                break;
                            }
                        }
                        else
                        {
                            vtList[(noVt-1)/NO_BITS_UINT] = last;
                            
                            currVt = currVt > 0 ? currVt - 1 : noVar;
                            while (vtType[currVt] == SV || vtType[currVt] == DEL)
                                currVt--;
                            
                            if(vt[currVt].prePos != prePos )
                            {
                                refPos = prePos;
                                prePos--;
                                d=0;
                            }
                            else
                            {
                                d =1;
                            }
                            
                        }
                        
                        
                        
                    }
                    else
                    {
                        if (getRefFinal(refPos) == read1[preChToCheck-1])
                        {
                            refPos = prePos;
                            prePos--;
                            d =0;
                            --preChToCheck;
                        }
                        else
                        {
                            currPathValid = false;
                            break;
                        }
                        
                    }
                }
                ++prePos;
                prePosExact = prePos;
            }
        }
        if(preChToCheck > 0)
        {
            currPathValid = false;
        }
    }
    //////////////////
    
    
    
    
    //with errors from here...
    
    if(currPathValid && rightPos < wholeReadLen) //check right side
    {
        currVt = lastVt;
        {
            
            if(remainedToInsert )
            {
                //currVt--;
                if(vtType[currVt] == INS)
                {
                    LoadUInt(ins_sv_arr+vt[currVt].alt, insLen);
                    tempLen = vt[currVt].alt+sizeof(uint32_t);
                }
                else if(vtType[currVt] == SV)//SV
                {
                    LoadUInt(ins_sv_arr+vt[currVt].alt+sizeof(uint32_t), insLen);
                    tempLen = vt[currVt].alt+(sizeof(uint32_t)<<1);
                }
                else
                {
                    cout << "  Should be SV or INS! " <<  lastVt << " " << prePos <<" " << remainedToInsert << "\n";
                    exit(9);
                }
                
                
                
                while (rightPos < wholeReadLen && remainedToInsert > 0)
                {
                    if ((ins_sv_arr[tempLen+insLen-remainedToInsert]) == wholeRead[rightPos])
                    {
                        ++rightPos; --remainedToInsert;
                    }
                    else
                    {
                        if(error == maxError)
                        {
                            currPathValid = false;
                            break;
                        }
                        else
                        {
                            ++error;
                            ++rightPos; --remainedToInsert;
                        }
                    }
                }
                lastVt++;
                currVt++;
                
            }
            
        }
        if(currPathValid)
            while(rightPos != wholeReadLen)
            {
                if(nextRefPos >= refSize)
                {
                    currPathValid = false;
                    break;
                }
                while(vt[currVt].prePos+1 < nextRefPos && rightPos != wholeReadLen)
                {
                 
                    
                    if(vtType[currVt] == DEL || vtType[currVt] == SV)
                    {
                        if(vtType[currVt] == DEL)
                        {
                            delLen = vt[currVt].alt;
                            insLen = 0;
                        }
                        else
                        {
                            LoadUInt(ins_sv_arr+vt[currVt].alt, delLen);
                            LoadUInt(ins_sv_arr+vt[currVt].alt+sizeof(uint32_t), insLen);
                           
                        }
                        
                        
                        newRefPos = vt[currVt].prePos+1+delLen > nextRefPos ? vt[currVt].prePos+1+delLen : nextRefPos;
                        
                        if(newRefPos > nextRefPos)
                        {
                            ///
                            (vtToCheck)[pos_vtToCheck++] = currVt;
                            if(pos_vtToCheck == max_vtToCheck)
                            {
                                max_vtToCheck = max_vtToCheck * 2;
                                
                                temp_vtToCheck = ( uint32_t *) realloc (vtToCheck, max_vtToCheck * sizeof( uint32_t));
                                if (NULL == temp_vtToCheck)
                                {
                                    printf("SA: Error in realloc, the message is - %s\n max = %d\n", strerror(errno), max_vtToCheck);
                                    fflush(stdout);
                                    //free(result);
                                    exit(8);
                                }
                                else
                                {
                                    vtToCheck = temp_vtToCheck;
                                }
                            }
                            ///
                            
                            
                            
                            
                            last = vtList[noVt/NO_BITS_UINT];
                            vtList[noVt/NO_BITS_UINT] = vtList[noVt/NO_BITS_UINT] | (1U << (NO_BITS_UINT -  (noVt&NO_BITS_UINT_MINUS_ONE) -  1U));
                            ++noVt;
                            if(noVt == vtListSize * NO_BITS_UINT)
                            {
                                vtListSize++;
                                uint32_t * tempList = (uint32_t *) realloc(vtList, vtListSize * sizeof(uint32_t));
                                if (NULL ==tempList) {
                                    printf("Error in realloc, the message is - %s\n", strerror(errno));
                                    exit(8);
                                }
                                else
                                {
                                    vtList = tempList;
                                }
                                vtList[vtListSize-1]=0;
                            }
                            
                            
                            
                            if(insLen == 0)
                            {
                                pos = checkVTconfiguration1(prePos, noVt, firstVar, currVt+1, currVt+1,firstVarRead1,  refPos, chCovered, read1, readLen, delta, preChToCheck,  &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  wholeRead, wholeReadLen, leftPos, rightPos, error, maxError, newRefPos, gap, remainedToInsert, startPosFoundExact, prePosExact, prevDelVar);
                            }
                            else
                            {
                                
                                ///
                                tempError = error;
                                goodPath=true;
                                for(ch = 0; ch < insLen && rightPos+ch != wholeReadLen; ch++)
                                {
                                    if(wholeRead[rightPos+ch] != ins_sv_arr[vt[currVt].alt+(sizeof(uint32_t)<<1)+ch])
                                    {
                                        if(tempError == maxError)
                                        {
                                            goodPath = false;
                                            break;
                                        }
                                        else
                                        {
                                            ++tempError;
                                        }
                                    }
                                }
                                if(goodPath)
                                    pos = checkVTconfiguration1(prePos, noVt, firstVar, currVt+1, currVt+1, firstVarRead1, refPos, chCovered, read1, readLen, delta, preChToCheck, &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  wholeRead, wholeReadLen, leftPos, rightPos+ch, tempError, maxError, newRefPos, gap, insLen-ch, startPosFoundExact, prePosExact, prevDelVar);
                                
                                
                                ///
                            }
                            
                            vtList[(noVt-1)/NO_BITS_UINT] = last;
                            
                        }
                        
                        
                        
                    }
                    
                    
                    currVt++;
                    lastVt++;
                    
                }
                
                if( rightPos != wholeReadLen)
                {
                    if( vt[currVt].prePos+1 == nextRefPos) //currVt < noVar &&
                    {
                        (vtToCheck)[(pos_vtToCheck)++] = currVt;
                        if(pos_vtToCheck == max_vtToCheck)
                        {
                            max_vtToCheck = max_vtToCheck * 2;
                            
                            temp_vtToCheck = ( uint32_t *) realloc (vtToCheck, max_vtToCheck * sizeof( uint32_t));
                            if (NULL == temp_vtToCheck)
                            {
                                printf("SA: Error in realloc, the message is - %s\n max = %d\n", strerror(errno), max_vtToCheck);
                                fflush(stdout);
                                //free(result);
                                exit(8);
                            }
                            else
                            {
                                vtToCheck = temp_vtToCheck;
                            }
                        }
                        ///
                        
                        
                        
                        last = vtList[noVt/NO_BITS_UINT];
                        vtList[noVt/NO_BITS_UINT] = vtList[noVt/NO_BITS_UINT] | (1U << (NO_BITS_UINT -  (noVt&NO_BITS_UINT_MINUS_ONE) -  1U));
                        ++noVt;
                        
                        if(noVt == vtListSize * NO_BITS_UINT)
                        {
                            vtListSize++;
                            uint32_t * tempList = (uint32_t *) realloc(vtList, vtListSize * sizeof(uint32_t));
                            if (NULL ==tempList) {
                                printf("Error in realloc, the message is - %s\n", strerror(errno));
                                exit(8);
                            }
                            else
                            {
                                vtList = tempList;
                            }
                            vtList[vtListSize-1]=0;
                        }
                        
                        
                        
                        switch (vtType[currVt])
                        {
                            case SNP:
                                //introduce variant
                                if(wholeRead[rightPos] == vt[currVt].alt)
                                {
                                    pos = checkVTconfiguration1(prePos, noVt, firstVar, currVt+1, currVt+1, firstVarRead1, refPos, chCovered, read1, readLen,delta, preChToCheck, &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  wholeRead, wholeReadLen, leftPos, rightPos+1, error, maxError, nextRefPos+1, gap, remainedToInsert, startPosFoundExact, prePosExact, prevDelVar);
                                }
                                else
                                {
                                    if(error < maxError)
                                    {
                                        pos = pos = checkVTconfiguration1(prePos, noVt, firstVar, currVt+1, currVt+1,firstVarRead1,  refPos, chCovered, read1, readLen,delta, preChToCheck, &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  wholeRead, wholeReadLen, leftPos, rightPos+1, error+1, maxError, nextRefPos+1, gap, remainedToInsert, startPosFoundExact, prePosExact, prevDelVar);
                                    }
                                    
                                }
                                
                                break;
                            case DEL:
                                //introduce variant
                               
                            
                                pos = checkVTconfiguration1(prePos, noVt, firstVar, currVt+1, currVt+1, firstVarRead1, refPos, chCovered, read1, readLen,delta, preChToCheck, &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  wholeRead, wholeReadLen, leftPos, rightPos, error, maxError, nextRefPos+vt[currVt].alt, gap, remainedToInsert, startPosFoundExact, prePosExact, prevDelVar);
                                
                            
                                break;
                            case INS:
                                //introduce variant
                              
                                tempError = error;
                                goodPath=true;
                                LoadUInt(ins_sv_arr+vt[currVt].alt, insLen);
                                for(ch = 0; ch <  insLen && rightPos+ch < wholeReadLen; ch++)
                                {
                                    if(wholeRead[rightPos+ch] != ins_sv_arr[vt[currVt].alt+sizeof(uint32_t)+ch])
                                    {
                                        if(tempError == maxError)
                                        {
                                            goodPath=false;
                                            break;
                                        }
                                        else
                                            ++tempError;
                                    }
                                }
                                if(goodPath)
                                    pos = checkVTconfiguration1(prePos, noVt, firstVar, currVt+1, currVt+1, firstVarRead1, refPos, chCovered, read1, readLen, delta, preChToCheck, &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  wholeRead, wholeReadLen, leftPos, rightPos+ch, tempError, maxError, nextRefPos, gap, insLen-ch, startPosFoundExact, prePosExact, prevDelVar);
                                
                                
                                
                                break;
                            case SV:
                                LoadUInt(ins_sv_arr+vt[currVt].alt+sizeof(uint32_t), insLen);
                                
                                if(insLen  == 0)
                                {
                                
                                    LoadUInt(ins_sv_arr+vt[currVt].alt, delLen);
                                   
                                    
                                    pos = checkVTconfiguration1(prePos, noVt, firstVar, currVt+1, currVt+1, firstVarRead1, refPos , chCovered, read1, readLen, delta,preChToCheck, &result, pos,max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  wholeRead, wholeReadLen, leftPos, rightPos, error, maxError, nextRefPos+delLen, gap, remainedToInsert, startPosFoundExact, prePosExact, prevDelVar);
                                    
                                    
                                }
                                else    // deletion+insertion
                                {
                                    
                                    tempError = error;
                                    goodPath=true;
                                    for(ch = 0; ch < insLen && rightPos+ch < wholeReadLen; ch++)
                                    {
                                        if(wholeRead[rightPos+ch] != ins_sv_arr[vt[currVt].alt+(sizeof(uint32_t)<<1)+ch])
                                        {
                                            if(tempError == maxError)
                                            {
                                                goodPath=false;
                                                break;
                                            }
                                            else
                                                ++tempError;
                                            
                                        }
                                    }
                                    if(goodPath)
                                    {
                                        LoadUInt(ins_sv_arr+vt[currVt].alt, delLen);
                                        
                                        pos = checkVTconfiguration1(prePos, noVt, firstVar, currVt+1, currVt+1, firstVarRead1, refPos, chCovered,read1, readLen,delta, preChToCheck, &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  wholeRead, wholeReadLen, leftPos, rightPos+ch, tempError, maxError, nextRefPos+delLen, gap, insLen-ch, startPosFoundExact, prePosExact, prevDelVar);
                                    }
                                }
                                
                                break;
                        }
                        
                      
                        //don't introduce variant
                        currVt++;
                        lastVt++;
                        vtList[(noVt-1)/NO_BITS_UINT] = last;
                        
                    }
                    else
                    {
                        len = vt[currVt].prePos - nextRefPos + 1 < wholeReadLen-rightPos ? vt[currVt].prePos - nextRefPos + 1 : wholeReadLen-rightPos ;
                        {
                            {
                                for(int32_t i = 0; i < len; ++i)
                                {
                                    if(wholeRead[rightPos+i] != getRefFinal(nextRefPos+i))
                                    {
                                        if(error == maxError)
                                        {
                                            currPathValid = false;
                                            break;
                                        }
                                        else
                                        {
                                            ++error;
                                        }
                                    }
                                }
                                
                                if(!currPathValid)
                                    break;
                                
                                rightPos += len;
                                nextRefPos += len;
                            }
                        }
                    }
                }
            }
    }
    
    if(currPathValid && leftPos >= 0) //check left side
    {
        
        ///check characters from unfinished left insertion!
        
        
        if(gap > 0)
        {
            if(vtType[firstVar] == INS)
                tempLen = vt[firstVar].alt+sizeof(uint32_t);
            else //SV
                tempLen = vt[firstVar].alt+(sizeof(uint32_t)<<1);
            
            
            
            
            while (leftPos >= 0 && gap > 0)
            {
                if ((ins_sv_arr[tempLen+gap-1]) == wholeRead[leftPos])
                {
                    
                    --leftPos;
                    --gap;
                }
                else
                {
                    if(error == maxError)
                    {
                        currPathValid = false;
                        break;
                    }
                    else
                    {
                        ++error;
                        --leftPos;
                        --gap;
                    }
                }
            }
            
        }
        
        
        
        if(currPathValid)
        {
            if(prePos != 0xFFFFFFFF)
            {
                while(leftPos >= 0 && currPathValid && prePos != 0xFFFFFFFF)
                {
                    
                    //Check beginning!!
                    currVt = firstVar > 0 ? firstVar - 1 : noVar;
                    while (vtType[currVt] == SV || vtType[currVt] == DEL)
                        currVt--;
                    
                    
                    
                    
                    if(vt[currVt].prePos != prePos )
                    {
                        refPos = prePos;
                        prePos--;
                        d = 0;
                    }
                    else{
                        refPos = prePos+1;
                        d = 1;
                    }
                    
                        while(d + (leftPos) >= 0 && prePos != 0xFFFFFFFE)
                        {
                            
                            while(delFromRight[prevDelVar].nextRefPos > prePos + 2 && d+leftPos>=0)
                            {
                                
                                newRefPos = vt[delFromRight[prevDelVar].vtNo].prePos < refPos ? vt[delFromRight[prevDelVar].vtNo].prePos : refPos;
                                
                                if(newRefPos < refPos)
                                {
                                    if(vtType[delFromRight[prevDelVar].vtNo] == DEL)
                                    {
                                        insLen = 0;
                                    }
                                    else if(vtType[delFromRight[prevDelVar].vtNo] == SV)
                                    {
                                        
                                        LoadUInt(ins_sv_arr+vt[delFromRight[prevDelVar].vtNo].alt+sizeof(uint32_t), insLen);
                                         }
                                    else
                                    {
                                        insLen = 0;
                                    }
                                    
                                    ///
                                    (vtToCheck)[pos_vtToCheck++] = delFromRight[prevDelVar].vtNo;
                                    if(pos_vtToCheck == max_vtToCheck)
                                    {
                                        max_vtToCheck = max_vtToCheck * 2;
                                        
                                        temp_vtToCheck = ( uint32_t *) realloc (vtToCheck, max_vtToCheck * sizeof( uint32_t));
                                        if (NULL == temp_vtToCheck)
                                        {
                                            printf("SA: Error in realloc, the message is - %s\n max = %d\n", strerror(errno), max_vtToCheck);
                                            fflush(stdout);
                                            //free(result);
                                            exit(8);
                                        }
                                        else
                                        {
                                            vtToCheck = temp_vtToCheck;
                                        }
                                    }
                                    ///
                                    
                                    last = vtList[noVt/NO_BITS_UINT];
                                    vtList[noVt/NO_BITS_UINT] = vtList[noVt/NO_BITS_UINT] | (1U << (NO_BITS_UINT -  (noVt&NO_BITS_UINT_MINUS_ONE) -  1U));
                                    ++noVt;
                                    if(noVt == vtListSize * NO_BITS_UINT)
                                    {
                                        vtListSize++;
                                        uint32_t * tempList = (uint32_t *) realloc(vtList, vtListSize * sizeof(uint32_t));
                                        if (NULL ==tempList) {
                                            printf("Error in realloc, the message is - %s\n", strerror(errno));
                                            exit(8);
                                        }
                                        else
                                        {
                                            vtList = tempList;
                                        }
                                        vtList[vtListSize-1]=0;
                                    }
                                    
                                    
                                    
                                    if(insLen == 0)
                                    {
                                        pos = checkVTconfiguration1(newRefPos, noVt, delFromRight[prevDelVar].vtNo, currVt, lastVt,firstVarRead1,  newRefPos, chCovered, read1, readLen, delta, preChToCheck,  &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  wholeRead, wholeReadLen, leftPos, rightPos, error, maxError, nextRefPos, gap, remainedToInsert, startPosFoundExact, prePosExact, prevDelVar > 0 ? prevDelVar - 1 : 0);
                                    }
                                    else
                                    {
                                        
                                        ///
                                        tempError = error;
                                        int32_t c=0;
                                        goodPath=true;
                                        for(ch = insLen - 1; ch >= 0  && c <= leftPos+d; --ch)
                                        {
                                            if(wholeRead[leftPos-c+d] != ins_sv_arr[vt[delFromRight[prevDelVar].vtNo].alt+(sizeof(uint32_t)<<1)+ch])
                                            {
                                                if(tempError == maxError)
                                                {
                                                    goodPath=false;
                                                    break;
                                                }
                                                else
                                                    ++tempError;
                                            }
                                            ++c;
                                        }
                                        
                                        if(goodPath)
                                        {
                                            pos = checkVTconfiguration1(newRefPos, noVt, delFromRight[prevDelVar].vtNo, currVt,  lastVt, firstVarRead1, newRefPos, chCovered+c, read1, readLen, delta, preChToCheck, &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  wholeRead, wholeReadLen, leftPos-c+d, rightPos, tempError, maxError, nextRefPos, ch+1, remainedToInsert, startPosFoundExact, prePosExact, prevDelVar > 0 ? prevDelVar - 1 : 0);
                                     
                                        }
                                        
                                        
                                        ///
                                        
                                    }
                                    
                                    vtList[(noVt-1)/NO_BITS_UINT] = last;
                                    
                                }
                                
                                
                                
                                prevDelVar = prevDelVar > 0 ? prevDelVar - 1 : 0;  // --prevDelVar;
                                
                            }
                            
                            if(d+leftPos>=0 && prePos != 0xFFFFFFFE)
                            {
                                if(delFromRight[prevDelVar].nextRefPos == prePos + 2)
                                {
                                    
                                    if(vtType[delFromRight[prevDelVar].vtNo] == DEL)
                                    {
                                        insLen = 0;
                                    }
                                    else if(vtType[delFromRight[prevDelVar].vtNo] == SV)
                                    {
                                        LoadUInt(ins_sv_arr+vt[delFromRight[prevDelVar].vtNo].alt+sizeof(uint32_t), insLen);
                                   
                                    }
                                    else
                                    {
                                        insLen = 0;
                                    }
                                    
                                    (vtToCheck)[(pos_vtToCheck)++] = delFromRight[prevDelVar].vtNo;
                                    if(pos_vtToCheck == max_vtToCheck)
                                    {
                                        max_vtToCheck = max_vtToCheck * 2;
                                        
                                        temp_vtToCheck = ( uint32_t *) realloc (vtToCheck, max_vtToCheck * sizeof( uint32_t));
                                        if (NULL == temp_vtToCheck)
                                        {
                                            printf("SA: Error in realloc, the message is - %s\n max = %d\n", strerror(errno), max_vtToCheck);
                                            fflush(stdout);
                                            //free(result);
                                            exit(8);
                                        }
                                        else
                                        {
                                            vtToCheck = temp_vtToCheck;
                                        }
                                    }
                                    last = vtList[noVt/NO_BITS_UINT];
                                    vtList[noVt/NO_BITS_UINT] = vtList[noVt/NO_BITS_UINT] | (1U << (NO_BITS_UINT -  (noVt&NO_BITS_UINT_MINUS_ONE) -  1U));
                                    ++noVt;
                                    if(noVt == vtListSize * NO_BITS_UINT)
                                    {
                                        vtListSize++;
                                        uint32_t * tempList = (uint32_t *) realloc(vtList, vtListSize * sizeof(uint32_t));
                                        if (NULL ==tempList) {
                                            printf("Error in realloc, the message is - %s\n", strerror(errno));
                                            exit(8);
                                        }
                                        else
                                        {
                                            vtList = tempList;
                                        }
                                        vtList[vtListSize-1]=0;
                                    }
                                    
                                    
                                    if(vtType[delFromRight[prevDelVar].vtNo] == DEL)
                                    {
                                        
                                        
                                        pos = checkVTconfiguration1(vt[delFromRight[prevDelVar].vtNo].prePos, noVt, delFromRight[prevDelVar].vtNo, currVt, lastVt, firstVarRead1, vt[delFromRight[prevDelVar].vtNo].prePos, chCovered, read1, readLen,delta, preChToCheck, &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  wholeRead, wholeReadLen, leftPos, rightPos, error, maxError, nextRefPos, gap, remainedToInsert, startPosFoundExact, prePosExact, prevDelVar > 0 ? prevDelVar - 1 : 0);
                                        
                                        
                                    }
                                    else //SV
                                    {
                                        if(insLen  == 0)
                                        {
                                            
                                            
                                            pos = checkVTconfiguration1(vt[delFromRight[prevDelVar].vtNo].prePos, noVt, delFromRight[prevDelVar].vtNo, currVt, lastVt,firstVarRead1,  vt[delFromRight[prevDelVar].vtNo].prePos, chCovered, read1, readLen,delta, preChToCheck, &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  wholeRead, wholeReadLen, leftPos, rightPos, error, maxError, nextRefPos, gap, remainedToInsert, startPosFoundExact, prePosExact, prevDelVar > 0 ? prevDelVar - 1 : 0);
                                            
                                        }
                                        else    // deletion+insertion
                                        {
                                            
                                            tempError = error;
                                            int32_t c = 0;
                                            goodPath=true;
                                            for(ch = insLen -1 ; ch >=0  && c <= leftPos+d; --ch)
                                            {
                                                if(wholeRead[leftPos-c+d] != ins_sv_arr[vt[delFromRight[prevDelVar].vtNo].alt+(sizeof(uint32_t)<<1)+ch])
                                                {
                                                    if(tempError == maxError)
                                                        goodPath=false;
                                                    else
                                                        ++tempError;
                                                }
                                                ++c;
                                            }
                                            if(goodPath)
                                            {
                                                pos = checkVTconfiguration1(vt[delFromRight[prevDelVar].vtNo].prePos, noVt, delFromRight[prevDelVar].vtNo, currVt,  lastVt, firstVarRead1, vt[delFromRight[prevDelVar].vtNo].prePos, chCovered+c, read1, readLen, delta, preChToCheck, &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  wholeRead, wholeReadLen, leftPos-c+d, rightPos, tempError, maxError, nextRefPos, ch+1, remainedToInsert, startPosFoundExact, prePosExact, prevDelVar > 0 ? prevDelVar - 1 : 0);
                                        
                                            }
                                        }
                                        
                                    }
                                    
                                    vtList[(noVt-1)/NO_BITS_UINT] = last;
                                    prevDelVar = prevDelVar > 0 ? prevDelVar - 1 : 0; // prevDelVar--;
                                }
                                else if (vt[currVt].prePos == prePos)
                                {
                                    (vtToCheck)[(pos_vtToCheck)++] = currVt;
                                    if(pos_vtToCheck == max_vtToCheck)
                                    {
                                        max_vtToCheck = max_vtToCheck * 2;
                                        
                                        temp_vtToCheck = ( uint32_t *) realloc (vtToCheck, max_vtToCheck * sizeof( uint32_t));
                                        if (NULL == temp_vtToCheck)
                                        {
                                            printf("SA: Error in realloc, the message is - %s\n max = %d\n", strerror(errno), max_vtToCheck);
                                            fflush(stdout);
                                            //free(result);
                                            exit(8);
                                        }
                                        else
                                        {
                                            vtToCheck = temp_vtToCheck;
                                        }
                                    }
                                    last = vtList[noVt/NO_BITS_UINT];
                                    vtList[noVt/NO_BITS_UINT] = vtList[noVt/NO_BITS_UINT] | (1U << (NO_BITS_UINT -  (noVt&NO_BITS_UINT_MINUS_ONE) -  1U));
                                    ++noVt;
                                    if(noVt == vtListSize * NO_BITS_UINT)
                                    {
                                        vtListSize++;
                                        uint32_t * tempList = (uint32_t *) realloc(vtList, vtListSize * sizeof(uint32_t));
                                        if (NULL ==tempList) {
                                            printf("Error in realloc, the message is - %s\n", strerror(errno));
                                            exit(8);
                                        }
                                        else
                                        {
                                            vtList = tempList;
                                        }
                                        vtList[vtListSize-1]=0;
                                    }
                                    
                                    if(vtType[currVt] == SNP)
                                    {
                                        if (vt[currVt].alt == wholeRead[leftPos])
                                        {
                                            pos = checkVTconfiguration1(prePos, noVt, currVt, currVt,  lastVt, firstVarRead1, refPos, chCovered, read1, readLen, delta, preChToCheck, &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  wholeRead, wholeReadLen, leftPos-1, rightPos, error, maxError, nextRefPos, gap, remainedToInsert, startPosFoundExact, prePosExact, prevDelVar);
                                        }
                                        else
                                        {
                                            if(error < maxError)
                                            {
                                                pos = checkVTconfiguration1(prePos, noVt, currVt, currVt,  lastVt, firstVarRead1, refPos, chCovered, read1, readLen, delta, preChToCheck, &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  wholeRead, wholeReadLen, leftPos-1, rightPos, error+1, maxError, nextRefPos, gap, remainedToInsert, startPosFoundExact, prePosExact, prevDelVar);
                                            }
                                        }
                                        
                                        
                                    }
                                    else //(vtType[currVt] == INS)
                                    {
                                        tempError = error;
                                        
                                        int32_t c = 0;
                                        ch = 0;
                                        goodPath=true;
                                        if (d == 0 && getRefFinal(refPos) != wholeRead[leftPos])
                                        {   if(tempError == maxError)
                                            goodPath=false;
                                        else
                                            ++tempError;
                                        }
                                        else
                                        {
                                            LoadUInt(ins_sv_arr+vt[currVt].alt, insLen);
                                            
                                            for(ch = insLen-1; ch >= 0  && c+1 <= leftPos+d; --ch)
                                            {
                                                if(wholeRead[leftPos+d-1-c] != ins_sv_arr[vt[currVt].alt+sizeof(uint32_t)+ch])
                                                {
                                                    if(tempError == maxError)
                                                    {
                                                        goodPath=false;
                                                        break;
                                                    }
                                                    else
                                                        ++tempError;
                                                    
                                                }
                                                ++c;
                                            }
                                        }
                                        if(goodPath)
                                            pos = checkVTconfiguration1(prePos, noVt, currVt, currVt, lastVt, firstVarRead1,  refPos, chCovered, read1, readLen, delta, preChToCheck, &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  wholeRead, wholeReadLen, leftPos+d-c-1, rightPos, tempError, maxError, nextRefPos, ch+1, remainedToInsert, startPosFoundExact, prePosExact, prevDelVar);
                                        
                                        
                                    }
                                    
                                    firstVar = currVt;
                                    vtList[(noVt-1)/NO_BITS_UINT] = last;
                                    
                                    if(d == 0)
                                    {
                                        
                                        
                                        if(getRefFinal(refPos) == wholeRead[leftPos])
                                        {
                                            
                                            
                                            currVt = currVt > 0 ? currVt - 1 : noVar;
                                            while (vtType[currVt] == SV || vtType[currVt] == DEL)
                                                currVt--;
                                            
                                            if(vt[currVt].prePos != prePos )//&&
                                            {
                                                refPos = prePos;
                                                prePos--;
                                                d =0;
                                                --leftPos;
                                            }
                                            else
                                            {
                                                d = 1;
                                                --leftPos;
                                            }
                                            
                                        }
                                        else
                                        {
                                            if(error == maxError)
                                            {
                                                currPathValid = false;
                                                break;
                                            }
                                            else
                                            {
                                                ++error;
                                                currVt = currVt > 0 ? currVt - 1 : noVar;
                                                while (vtType[currVt] == SV || vtType[currVt] == DEL)
                                                    currVt--;
                                                
                                                if(vt[currVt].prePos != prePos )
                                                {
                                                    refPos = prePos;
                                                    prePos--;
                                                    d =0;
                                                    leftPos--;
                                                }
                                                else
                                                {
                                                    d = 1;
                                                    leftPos--;
                                                }
                                                
                                            }
                                            
                                        }
                                    }
                                    else
                                    {
                                        vtList[(noVt-1)/NO_BITS_UINT] = last;
                                        
                                        currVt = currVt > 0 ? currVt - 1 : noVar;
                                        while (vtType[currVt] == SV || vtType[currVt] == DEL)
                                            currVt--;
                                        
                                        if(vt[currVt].prePos != prePos )
                                        {
                                            refPos = prePos;
                                            prePos--;
                                            d=0;
                                        }
                                        else
                                        {
                                            d =1;
                                        }
                                        
                                    }
                                    
                                    
                                    
                                }
                                else
                                {
                                    if (getRefFinal(refPos) == wholeRead[leftPos])
                                    {
                                        refPos = prePos;
                                        prePos--;
                                        d =0;
                                        --leftPos;
                                    }
                                    else
                                    {
                                        if(error == maxError)
                                        {
                                            currPathValid = false;
                                            break;
                                        }
                                        else
                                        {
                                            ++error;
                                            refPos = prePos;
                                            prePos--;
                                            d =0;
                                            leftPos--;
                                            
                                        }
                                    }
                                    
                                }
                            }
                        }
                        ++prePos;
                    
                }
                if(leftPos >= 0)
                {
                    currPathValid = false;
                }
            }
            else
            {
                currPathValid = false;
            }
        }
    }
    
    
    
    //deletions at the beginnig?
    if(currPathValid)
    {
        while(delFromRight[prevDelVar].nextRefPos >= prePos + 1 && prePos != 0xFFFFFFFF)  // can't be >, unless del in del region
        {
            (vtToCheck)[(pos_vtToCheck)++] = delFromRight[prevDelVar].vtNo;
            if(pos_vtToCheck == max_vtToCheck)
            {
                max_vtToCheck = max_vtToCheck * 2;
                
                temp_vtToCheck = ( uint32_t *) realloc (vtToCheck, max_vtToCheck * sizeof( uint32_t));
                if (NULL == temp_vtToCheck)
                {
                    printf("SA: Error in realloc, the message is - %s\n max = %d\n", strerror(errno), max_vtToCheck);
                    fflush(stdout);
                    //free(result);
                    exit(8);
                }
                else
                {
                    vtToCheck = temp_vtToCheck;
                }
            }
            last = vtList[noVt/NO_BITS_UINT];
            vtList[noVt/NO_BITS_UINT] = vtList[noVt/NO_BITS_UINT] | (1U << (NO_BITS_UINT -  (noVt&NO_BITS_UINT_MINUS_ONE) -  1U));
            ++noVt;
            if(noVt == vtListSize * NO_BITS_UINT)
            {
                vtListSize++;
                uint32_t * tempList = (uint32_t *) realloc(vtList, vtListSize * sizeof(uint32_t));
                if (NULL ==tempList) {
                    printf("Error in realloc, the message is - %s\n", strerror(errno));
                    exit(8);
                }
                else
                {
                    vtList = tempList;
                }
                vtList[vtListSize-1]=0;
            }
            
            
            if(vtType[delFromRight[prevDelVar].vtNo] == DEL)
            {
                
                
                pos = checkVTconfiguration1(vt[delFromRight[prevDelVar].vtNo].prePos, noVt, delFromRight[prevDelVar].vtNo, currVt, lastVt, firstVarRead1, vt[delFromRight[prevDelVar].vtNo].prePos, chCovered, read1, readLen,delta, preChToCheck, &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  wholeRead, wholeReadLen, leftPos, rightPos, error, maxError, nextRefPos, gap, remainedToInsert, startPosFoundExact, prePosExact, prevDelVar > 0 ? prevDelVar - 1 : 0);
                
                
            }
            else //SV
            {
                if(insLen  == 0)
                {
                    
                    
                    pos = checkVTconfiguration1(vt[delFromRight[prevDelVar].vtNo].prePos, noVt, delFromRight[prevDelVar].vtNo, currVt, lastVt,firstVarRead1,  vt[delFromRight[prevDelVar].vtNo].prePos, chCovered, read1, readLen,delta, preChToCheck, &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  wholeRead, wholeReadLen, leftPos, rightPos, error, maxError, nextRefPos, gap, remainedToInsert, startPosFoundExact, prePosExact, prevDelVar > 0 ? prevDelVar - 1 : 0);
                    
                }
                else    // deletion+insertion
                {
                    
                    tempError = error;
                    int32_t c = 0;
                    //ch = 0;
                    goodPath=true;
                    for(ch = insLen -1 ; ch >=0  && c <= leftPos; --ch)
                    {
                        if(wholeRead[leftPos-c] != ins_sv_arr[vt[delFromRight[prevDelVar].vtNo].alt+(sizeof(uint32_t)<<1)+ch])
                        {
                            if(tempError == maxError)
                                goodPath=false;
                            else
                                ++tempError;
                        }
                        ++c;
                    }
                    if(goodPath)
                    {
                        pos = checkVTconfiguration1(vt[delFromRight[prevDelVar].vtNo].prePos, noVt, delFromRight[prevDelVar].vtNo, currVt,  lastVt, firstVarRead1, vt[delFromRight[prevDelVar].vtNo].prePos, chCovered+c, read1, readLen, delta, preChToCheck, &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  wholeRead, wholeReadLen, leftPos-c+d, rightPos, tempError, maxError, nextRefPos, ch+1, remainedToInsert, startPosFoundExact, prePosExact, prevDelVar > 0 ? prevDelVar - 1 : 0);
                        
                    }
                }
                
            }
            
            vtList[(noVt-1)/NO_BITS_UINT] = last;
            prevDelVar = prevDelVar > 0 ? prevDelVar - 1 : 0;
          
        }
    }
    
    
    
    
    
    
    
    
    if(currPathValid)
    {
        
        
        bool suffix, vtFound;
      
        uint64_t  *tempResult;
        int32_t back, dels;

        if(pos_vtToCheck == 0 && prevDelCount == 0 && firstSV>lastSV)
        {
            result[pos++] = 0xFFFFFFFFFFFFFFFE;
            
            if (pos  == *max)
            {
                *max = (*max) * 2;
               
                tempResult = ( uint64_t  *) realloc (result, (*max) * sizeof( uint64_t ));
                if (NULL == tempResult)
                {
                    printf("SA: Error in realloc, the message is - %s\n max = %llu\n", strerror(errno), *max);
                    fflush(stdout);
                    //free(result);
                    exit(8);
                }
                else
                {
                    result = tempResult;
                }
            }
            
            // pos
            result[pos++] = prePos;
            if (pos == *max)
            {
                *max = (*max) * 2 ;
                
                tempResult = ( uint64_t  *) realloc (result, (*max) * sizeof( uint64_t ));
                if (NULL == tempResult)
                {
                    printf("SA: Error in realloc, the message is - %s\n max = %llu\n", strerror(errno), *max);
                    fflush(stdout);
                    //free(result);
                    exit(8);
                }
                else
                {
                    result = tempResult;
                }
            }
            
            
          
        }
        else
        {
            
            uint32_t aux;
            
            A_memset(resV, 0xFF, resSize<<3);
            
            
            dels = prevDelCount;
            back = 1;
            while(dels)
            {
                if(vtType[firstVarRead1 - back] == DEL)
                {
                    dels--;
                    if(vt[firstVarRead1 - back].prePos + vt[firstVarRead1 - back].alt > prePosExact)
                    {
                       
                        
                        aux = all_var[(firstVarRead1 - back)]*parts;
                        for(int32_t r=0; r < resSize; r++)
                            resV[r] = resV[r] & ~uniq_parts[uniq_var[aux+r/blocking]*blocking+(r)%blocking];
                        
                    }
                }
                back++;
            }
            
            if(firstSV < noSV)
                for(int32_t sv = lastSV; sv >= (int32_t)firstSV; sv--)
                {
                    LoadUInt(ins_sv_arr+vt[listSV[sv]].alt, delLen);
                    
                    if(vt[listSV[sv]].prePos +  delLen > prePosExact)
                    {

                        
                        aux = all_var[listSV[sv]]*parts;
                        for(int32_t r=0; r < resSize; r++)
                            resV[r] = resV[r] & ~uniq_parts[uniq_var[aux+r/blocking]*blocking+(r)%blocking];
                    }
                }
            
            for(int32_t variant=0; variant < pos_vtToCheck; variant++)
            {
                aux = all_var[(vtToCheck)[variant]]*parts;
                suffix = vtList[variant/NO_BITS_UINT] & (1 << (NO_BITS_UINT - variant%NO_BITS_UINT - 1));
                if(suffix)
                    for(int32_t r=0; r < resSize; r++)
                        resV[r] = resV[r] & uniq_parts[uniq_var[aux+r/blocking]*blocking+(r)%blocking];
                else
                    for(int32_t r=0; r < resSize; r++)
                        resV[r] = resV[r] & ~uniq_parts[uniq_var[aux+r/blocking]*blocking+(r)%blocking];
              
            }
            
            vtFound = false;
            for ( int32_t r = 0; r < resSize; r++)
            {
                if(resV[r] > 0)
                {
                    vtFound = true;
                    break;
                }
            }
            
            if(vtFound)
            {
                
                result[pos++] = 0xFFFFFFFFFFFFFFFF;
                if (pos  == *max)
                {
                    *max = (*max) * 2;
                    tempResult = ( uint64_t  *) realloc (result, (*max) * sizeof( uint64_t ));
                    if (NULL == tempResult)
                    {
                        printf("SA: Error in realloc, the message is - %s\n max = %llu\n", strerror(errno), *max);
                        fflush(stdout);
                        //free(result);
                        exit(8);
                    }
                    else
                    {
                        result = tempResult;
                    }
                }
                
                // pos
                result[pos++] = prePos;
                if (pos  == *max)
                {
                    *max = (*max) * 2;
                    
                    
                    tempResult = ( uint64_t  *) realloc (result, (*max) * sizeof( uint64_t ));
                    if (NULL == tempResult)
                    {
                        printf("SA: Error in realloc, the message is - %s\n max = %llu\n", strerror(errno), *max);
                        fflush(stdout);
                        //free(result);
                        exit(8);
                    }
                    else
                    {
                        result = tempResult;
                    }
                }
                
                
                if (pos + resSize >= *max)
                {
                    *max = (*max) * 2 > (*max) + resSize ? (*max) * 2 : (*max) + 2*resSize ;
        
                    
                    tempResult = ( uint64_t  *) realloc (result, (*max) * sizeof( uint64_t ));
                    if (NULL == tempResult)
                    {
                        printf("SA: Error in realloc, the message is - %s\n max = %llu\n", strerror(errno), *max);
                        fflush(stdout);
                        //free(result);
                        exit(8);
                    }
                    else
                    {
                        result = tempResult;
                    }
                }
                
                A_memcpy(result+pos, resV, resSize<<3); pos = pos + resSize;
                
            }
            
        }

        
        
    }
    
    //  free(vtList);
    *passed_result=result;
    return pos;
}


uint32_t getCoveredVariants(uint32_t  prevPos,  uint32_t * vtListLoc, uint32_t * firstVar, uint32_t gap, uint32_t * chCov, uint32_t * nextRefPos, uint32_t * currVT,  uint32_t * pos_vtToCheck, uint32_t  lenToCheck , uint32_t * remainedToInsert)
{
    uint32_t *temp_vtToCheck;
    *pos_vtToCheck = 0;
    uint32_t chCovered = 0, newCh;
    uint32_t refPos = prevPos+1;
    // bool firstVT = true;
    ///BINARY SEARCH FOR FIRST VT
    uint32_t vtNo = 0,  m, l, r; //l = 0, r = noVar-1, m;
    
    uint32_t temp = (prevPos != 0xFFFFFFFF)? prevPos/VT_LUT_STEP: 0;
    if(refSize >= VT_LUT_STEP)
    {
        l = VT_LUT[temp++];
        r = VT_LUT[temp];
    }
    else
    {
        l = 0;
        r = noVar;
    }

    

    
    while (l < r)
    {
        m = ((l+r) / 2);
        if ((vt[m].prePos < prevPos || vt[m].prePos == 0xFFFFFFFF) && prevPos != 0xFFFFFFFF)
        {
            l = m + 1;
        }
        else
        {
            r = m;
        }
    }
    vtNo = l;
    *firstVar = vtNo;
    
    lenToCheck =  lenToCheck < kl ?  lenToCheck : kl;
    
    uint32_t vtCovered = 0;
    uint32_t newRefPos, startPos;
    
    
    uint32_t insLen, delLen,  tempLen;
    while (chCovered < lenToCheck)
    {
        
        while( (startPos = vt[vtNo].prePos+1) < refPos && chCovered < lenToCheck) //vtNo < noVar &&
        {
            
            if(vtType[vtNo] == SV || vtType[vtNo] == DEL)
            {
                if(vtType[vtNo] == DEL)
                {
                    delLen = vt[vtNo].alt;
                    insLen = 0;
                }
                else
                {
                    
                    LoadUInt(ins_sv_arr+vt[vtNo].alt, delLen);
                    LoadUInt(ins_sv_arr+vt[vtNo].alt+sizeof(uint32_t), insLen);
                }
                
                
                newRefPos = startPos+delLen > refPos ? startPos+delLen : refPos;
                
                if(newRefPos > refPos)
                {
                    
                    if(vtListLoc[vtCovered/NO_BITS_UINT] & (1 << (NO_BITS_UINT - vtCovered%NO_BITS_UINT - 1))) //variant found/introduced
                    {
                        
                        if(insLen == 0)
                        {
                            refPos = newRefPos;
                            
                        }
                        else
                        {
                            refPos = newRefPos;
                            
                            
                            int32_t back = 1;
                            int32_t prevIns = 0;
                            if(vtType[vtNo-back] == SV)
                            {
                                LoadUInt(ins_sv_arr+vt[vtNo-back].alt+sizeof(uint32_t), tempLen);
                            }
                            while(vt[vtNo].prePos == vt[vtNo-back].prePos && (vtType[vtNo-back] == INS || (vtType[vtNo-back] == SV && tempLen != 0)))
                            {
                                prevIns =  tempLen > prevIns?  tempLen : prevIns;
                                back++;
                            }
                            if(prevIns == 0 || insLen > prevIns)
                            {
                                
                                if(chCovered + insLen - prevIns <= lenToCheck)
                                {
                                    chCovered = chCovered + insLen - prevIns;
                                    
                                }
                                else
                                {
                                    *remainedToInsert = insLen - (lenToCheck-chCovered);
                                    chCovered = lenToCheck;//chCovered + insLen-gap;
                                    
                                    
                                    vtNo--;
                                }
                            
                            }
                            
                        }
                        
                    }
                    vtToCheck[(*pos_vtToCheck)++] = vtNo;
                    if(*pos_vtToCheck == max_vtToCheck)
                    {
                        max_vtToCheck = (max_vtToCheck) * 2;
                        
                        temp_vtToCheck = ( uint32_t *) realloc (vtToCheck, (max_vtToCheck) * sizeof( uint32_t));
                        if (NULL == temp_vtToCheck)
                        {
                            printf("SA: Error in realloc, the message is - %s\n max = %d\n", strerror(errno), max_vtToCheck);
                            fflush(stdout);
                            //free(result);
                            exit(8);
                        }
                        else
                        {
                            vtToCheck = temp_vtToCheck;
                        }
                        
                    }
                    vtCovered++;
                }
            }
            
            vtNo++;
            //cout << "shouldn't end up here with suffix1 (no var in del)\n";
        }
        
        if(chCovered < lenToCheck)
        {
            if(startPos > refPos)  //no variant at currPos //vtNo >= noVar ||
            {
                
                newCh = vt[vtNo].prePos+1-refPos;
                newCh = chCovered + newCh < lenToCheck ?  newCh : lenToCheck - chCovered;
                chCovered = chCovered + newCh ;
                refPos = refPos + newCh ;
                gap = 0;
                
            }
            else
            {
                vtToCheck[(*pos_vtToCheck)++] = vtNo;
                if(*pos_vtToCheck == max_vtToCheck)
                {
                    max_vtToCheck = (max_vtToCheck) * 2;
                    
                    temp_vtToCheck = ( uint32_t *) realloc (vtToCheck, (max_vtToCheck) * sizeof( uint32_t));
                    if (NULL == temp_vtToCheck)
                    {
                        printf("SA: Error in realloc, the message is - %s\n max = %d\n", strerror(errno), max_vtToCheck);
                        fflush(stdout);
                        //free(result);
                        exit(8);
                    }
                    else
                    {
                        vtToCheck = temp_vtToCheck;
                    }
                    
                }
                
                
                if(vtListLoc[vtCovered/NO_BITS_UINT] & (1 << (NO_BITS_UINT - vtCovered%NO_BITS_UINT - 1))) //variant found/introduced
                {
                    switch (vtType[vtNo])
                    {
                        case SNP:
                            //introduce variant
                            
                            chCovered++;
                            refPos++;
                            vtCovered++;
                            vtNo++;
                            gap = 0;
                            
                            break;
                        case DEL:
                            //introduce variant
                            refPos = refPos + vt[vtNo].alt;
                            vtCovered++;
                            vtNo++;
                            gap = 0;
                            
                            break;
                        case INS:
                            //introduce variant
                            LoadUInt(ins_sv_arr+vt[vtNo].alt, insLen);
                            if(chCovered + insLen-gap <= lenToCheck)
                            {
                                chCovered = chCovered + insLen-gap;
                                gap = 0;
                                vtCovered++;
                                vtNo++;
                            }
                            else
                            {
                                *remainedToInsert = insLen - (lenToCheck-chCovered)-gap;
                                chCovered = lenToCheck;//chCovered + insLen-gap;
                                gap = 0;
                                vtCovered++;
                                //vtNo++;
                            }
                            
                            break;
                        case SV:
                            
                            
                            
                            LoadUInt(ins_sv_arr+vt[vtNo].alt, delLen);
                            LoadUInt(ins_sv_arr+vt[vtNo].alt+sizeof(uint32_t), insLen);
                            //introduce variant
                            if(insLen == 0)
                            {
                                refPos = refPos + delLen;
                                vtCovered++;
                                vtNo++;
                                gap = 0;
                            }
                            else
                            {
                                if(chCovered + insLen-gap <= lenToCheck)
                                {
                                    chCovered = chCovered + insLen-gap;
                                    gap = 0;
                                    vtCovered++;
                                    vtNo++;
                                    refPos = refPos + delLen;
                                }
                                else
                                {
                                    *remainedToInsert = insLen - (lenToCheck-chCovered)-gap;
                                    chCovered = lenToCheck;
                                    gap = 0;
                                    vtCovered++;
                                    refPos = refPos + delLen;
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
                gap = 0;
            }
        }
    }
    
    *nextRefPos = refPos;
    
    *chCov = chCovered;
    *currVT = vtNo;
    
    return vtCovered;
}


uint32_t getCoveredVariants(uint32_t  prevPos,  uint32_t vtListLoc, uint32_t * firstVar, uint32_t gap, uint32_t * chCov, uint32_t * nextRefPos, uint32_t * pos_vtToCheck, uint32_t  lenToCheck, uint32_t * remainedToInsert)
{
    uint32_t insLen, delLen;
    
    uint32_t *temp_vtToCheck;
    *pos_vtToCheck = 0;
    
    uint32_t chCovered = 0, newCh;
    uint32_t refPos = prevPos+1;
    
    ///BINARY SEARCH FOR FIRST VT
    uint32_t vtNo = 0,  m, l, r;

    
    uint32_t temp = (prevPos != 0xFFFFFFFF)? prevPos/VT_LUT_STEP: 0;
    if(refSize >= VT_LUT_STEP)
    {
        l = VT_LUT[temp++];
        r = VT_LUT[temp];
    }
    else
    {
        l = 0;
        r = noVar;
    }

    
    while (l < r)
    {
        m = ((l+r) / 2);
        if ((vt[m].prePos < prevPos || vt[m].prePos == 0xFFFFFFFF) && prevPos != 0xFFFFFFFF)
        {
            l = m + 1;
        }
        else
        {
            r = m;
        }
    }
    vtNo = l;
    *firstVar = vtNo;

    
    uint32_t vtCovered = 0;
    
    lenToCheck =  lenToCheck < kl ?  lenToCheck : kl;
    
    while (chCovered <  lenToCheck)
    {
        
  
        if(chCovered < lenToCheck)
        {
            if( vt[vtNo].prePos+1 > refPos)  //no variant at currPos //vtNo >= noVar ||
            {
                
                newCh = vt[vtNo].prePos+1-refPos;
                newCh = chCovered + newCh < lenToCheck ?  newCh : lenToCheck - chCovered;
                chCovered = chCovered + newCh ;
                refPos = refPos + newCh ;
  
                
            }
            else
            {
                (vtToCheck)[(*pos_vtToCheck)++] = vtNo;
                if(*pos_vtToCheck == max_vtToCheck)
                {
                    max_vtToCheck = (max_vtToCheck) * 2;
                    
                    temp_vtToCheck = ( uint32_t *) realloc (vtToCheck, (max_vtToCheck) * sizeof( uint32_t));
                    if (NULL == temp_vtToCheck)
                    {
                        printf("SA: Error in realloc, the message is - %s\n max = %d\n", strerror(errno), max_vtToCheck);
                        fflush(stdout);
                        //free(result);
                        exit(8);
                    }
                    else
                    {
                        vtToCheck = temp_vtToCheck;
                    }
                    
                }
                
                
                if(vtListLoc & (1 << (NO_BITS_UINT - vtCovered - 1))) //variant found/introduced
                {
                    switch (vtType[vtNo])
                    {
                        case SNP:
                            //introduce variant
                            
                            chCovered++;
                            refPos++;
                            vtCovered++;
                            vtNo++;
                            
                            break;
                        case DEL:
                            //introduce variant
                            refPos = refPos + vt[vtNo].alt;
                            vtCovered++;
                            vtNo++;
                            
                            break;
                        case INS:
                            //introduce variant
                            LoadUInt(ins_sv_arr+vt[vtNo].alt, insLen);
                             if(chCovered + insLen-gap <= lenToCheck)
                            {
                                chCovered = chCovered + insLen-gap;
                                gap = 0;
                                vtCovered++;
                                vtNo++;
                            }
                            else
                            {
                                *remainedToInsert = insLen - (lenToCheck-chCovered)-gap;
                                chCovered = lenToCheck;
                                gap = 0;
                                vtCovered++;
                                //vtNo++;
                            }
                            break;
                        case SV:
                            LoadUInt(ins_sv_arr+vt[vtNo].alt, delLen);
                            
                            LoadUInt(ins_sv_arr+vt[vtNo].alt + sizeof(uint32_t), insLen);
                            
                            //introduce variant
                            if(insLen == 0)
                            {
                                refPos = refPos + delLen;
                                vtCovered++;
                                vtNo++;
                            }
                            else
                            {
                                if(chCovered + insLen-gap <= lenToCheck)
                                {
                                    chCovered = chCovered + insLen-gap;
                                    gap = 0;
                                    vtCovered++;
                                    vtNo++;
                                    refPos = refPos + delLen;
                                }
                                else
                                {
                                    *remainedToInsert = insLen - (lenToCheck-chCovered)-gap;
                                    chCovered = lenToCheck;//chCovered + insLen-gap;
                                    gap = 0;
                                    vtCovered++;
                                    refPos = refPos + delLen;
                                    //vtNo++;
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
                gap = 0;
            }
        }
    }
    
    *nextRefPos = refPos;
    *chCov = chCovered;
    
    return vtCovered;
}








uint64_t checkDelRegion(uint32_t currfFirstVar, uint32_t prevFirstVar, int32_t leftPos, int32_t gapFirstVar, int32_t error, int32_t maxError, unsigned char * wholeRead,  uint32_t  pos_vtToCheck, uint32_t vtCovered, int32_t  prev_firstSV, int32_t prev_lastSV, uint32_t  prev_prevDelCount, uint64_t ** passed_result, uint64_t  pos, uint64_t *max, uint32_t lastVt, uint32_t nextRefPos, uint32_t remainedToInsert, uint32_t chCovered, unsigned char * read1, int32_t readLen, uint32_t delta,  int32_t preChToCheck,  uint32_t wholeReadLen, uint32_t rightPos)
{
    //remember current vtList
    uint32_t * lastVtList = (uint32_t *)malloc(sizeof(uint32_t)*vtListSize);
    for(int i = 0;i < vtListSize; i++)
        lastVtList[i] = vtList[i];
    
    uint64_t * result = *passed_result;
    
    uint32_t delLenFirst, delLen;
    uint64_t * tempResult;
    uint32_t * temp_vtToCheck;
    uint32_t tempLen;
    bool include = true;
    //introduce (with errors) characters from gap from current variant
    if(leftPos >=0 && gapFirstVar >=0)
    {
        tempLen = vt[prevFirstVar].alt+(sizeof(uint32_t)<<1);
        while (leftPos >= 0 && gapFirstVar >= 0)
        {
            if ((ins_sv_arr[tempLen+gapFirstVar]) == wholeRead[leftPos])
            {
                
                --leftPos;
                --gapFirstVar;
            }
            else
            {
                if(error == maxError)
                {
                    include = false;
                    break;
                }
                else
                {
                    ++error;
                    --leftPos;
                    --gapFirstVar;
                }
            }
        }
    }
    
    ////////////////////////////////////////////////////
    
    if(include)
    {
        // checkConfiguration with current variant as first - set
        (vtToCheck)[pos_vtToCheck++] = currfFirstVar;
        if(pos_vtToCheck == max_vtToCheck)
        {
            max_vtToCheck = max_vtToCheck * 2;
            
            temp_vtToCheck = ( uint32_t *) realloc (vtToCheck, max_vtToCheck * sizeof( uint32_t));
            if (NULL == temp_vtToCheck)
            {
                printf("SA: Error in realloc, the message is - %s\n max = %d\n", strerror(errno), max_vtToCheck);
                fflush(stdout);
                //free(result);
                exit(8);
            }
            else
            {
                vtToCheck = temp_vtToCheck;
            }
        }
        ///
        vtList[vtCovered/NO_BITS_UINT] = vtList[vtCovered/NO_BITS_UINT] | (1U << (NO_BITS_UINT -  (vtCovered&NO_BITS_UINT_MINUS_ONE) -  1U));
        vtCovered++;
        if(vtCovered == vtListSize * NO_BITS_UINT)
        {
            vtListSize++;
            uint32_t * tempList = (uint32_t *) realloc(vtList, vtListSize * sizeof(uint32_t));
            if (NULL ==tempList) {
                printf("Error in realloc, the message is - %s\n", strerror(errno));
                exit(8);
            }
            else
            {
                vtList = tempList;
            }
            vtList[vtListSize-1]=0;
            
            tempList = (uint32_t *) realloc(lastVtList, vtListSize * sizeof(uint32_t));
            if (NULL ==tempList) {
                printf("Error in realloc, the message is - %s\n", strerror(errno));
                exit(8);
            }
            else
            {
                lastVtList = tempList;
            }
            lastVtList[vtListSize-1]=0;
        }
        
        ////////////////////////////////////////////////////
        
        //add to variant list deletionts that should not be set (as they are removing prevFirstVar)
        int back = 1;
        LoadUInt(ins_sv_arr+vt[prevFirstVar].alt, delLenFirst);
        while(prev_prevDelCount)
        {
            if(vtType[prevFirstVar - back] == DEL)
            {
                prev_prevDelCount--;
                //if this deleton delete whole prevFirstVar (go beyond it)  and is after current deletion
                
                if(vt[prevFirstVar - back].prePos + vt[prevFirstVar - back].alt > vt[prevFirstVar].prePos + delLenFirst \
                   && vt[prevFirstVar - back].prePos  > vt[currfFirstVar].prePos)
                {
                    //add it to vt list (not set)
                    (vtToCheck)[pos_vtToCheck++] = prevFirstVar - back;
                    if(pos_vtToCheck == max_vtToCheck)
                    {
                        max_vtToCheck = max_vtToCheck * 2;
                        
                        temp_vtToCheck = ( uint32_t *) realloc (vtToCheck, max_vtToCheck * sizeof( uint32_t));
                        if (NULL == temp_vtToCheck)
                        {
                            printf("SA: Error in realloc, the message is - %s\n max = %d\n", strerror(errno), max_vtToCheck);
                            fflush(stdout);
                            //free(result);
                            exit(8);
                        }
                        else
                        {
                            vtToCheck = temp_vtToCheck;
                        }
                    }
                    /// not set (keep zero, while increasing pos_vtToCheck)
                    vtCovered++;
                    if(vtCovered == vtListSize * NO_BITS_UINT)
                    {
                        vtListSize++;
                        uint32_t * tempList = (uint32_t *) realloc(vtList, vtListSize * sizeof(uint32_t));
                        if (NULL ==tempList) {
                            printf("Error in realloc, the message is - %s\n", strerror(errno));
                            exit(8);
                        }
                        else
                        {
                            vtList = tempList;
                        }
                        vtList[vtListSize-1]=0;
                        
                        tempList = (uint32_t *) realloc(lastVtList, vtListSize * sizeof(uint32_t));
                        if (NULL ==tempList) {
                            printf("Error in realloc, the message is - %s\n", strerror(errno));
                            exit(8);
                        }
                        else
                        {
                            lastVtList = tempList;
                        }
                        lastVtList[vtListSize-1]=0;
                    }
                }
            }
            back++;
        }
        if(prev_firstSV < noSV)
            for(int32_t sv = prev_lastSV; sv >= (int32_t)prev_firstSV; sv--)
            {
                LoadUInt(ins_sv_arr+vt[listSV[sv]].alt, delLen);
                
                //if this deleton delete whole prevFirstVar (go beyond it) and is after current deletion
                if(vt[listSV[sv]].prePos +  delLen > vt[prevFirstVar].prePos + delLenFirst \
                   && vt[prevFirstVar - back].prePos  > vt[currfFirstVar].prePos)
                {
                    //add it to vt list (not set)
                    (vtToCheck)[pos_vtToCheck++] = prevFirstVar - back;
                    if(pos_vtToCheck == max_vtToCheck)
                    {
                        max_vtToCheck = max_vtToCheck * 2;
                        
                        temp_vtToCheck = ( uint32_t *) realloc (vtToCheck, max_vtToCheck * sizeof( uint32_t));
                        if (NULL == temp_vtToCheck)
                        {
                            printf("SA: Error in realloc, the message is - %s\n max = %d\n", strerror(errno), max_vtToCheck);
                            fflush(stdout);
                            //free(result);
                            exit(8);
                        }
                        else
                        {
                            vtToCheck = temp_vtToCheck;
                        }
                    }
                    /// not set (keep zero, while increasing pos_vtToCheck)
                    vtCovered++;
                    if(vtCovered == vtListSize * NO_BITS_UINT)
                    {
                        vtListSize++;
                        uint32_t * tempList = (uint32_t *) realloc(vtList, vtListSize * sizeof(uint32_t));
                        if (NULL ==tempList) {
                            printf("Error in realloc, the message is - %s\n", strerror(errno));
                            exit(8);
                        }
                        else
                        {
                            vtList = tempList;
                        }
                        vtList[vtListSize-1]=0;
                        
                        tempList = (uint32_t *) realloc(lastVtList, vtListSize * sizeof(uint32_t));
                        if (NULL ==tempList) {
                            printf("Error in realloc, the message is - %s\n", strerror(errno));
                            exit(8);
                        }
                        else
                        {
                            lastVtList = tempList;
                        }
                        lastVtList[vtListSize-1]=0;
                    }
                }
            }
        
        ////////////////////////////////////////////////////
        
        /** Caalculate new prevDelCount, firstSV, lastSV = to check in checkVTconfiguration1 **/
        /** Check if not in deleted region - find meaningful SVs and DELs **/
        int32_t lastSV = 0, firstSV = 0;
        uint32_t  prevDelCount = 0;
        bool firstVarNotZero = prevFirstVar > 0 ? true : false;
        
        
        int32_t cSV = prevSV[currfFirstVar];//noSV - 1;
        while(cSV >= 0 && vt[listSV[cSV]].prePos >= vt[currfFirstVar].prePos)
            cSV--;
        lastSV = cSV;
        if(cSV < 0)
        {
            firstSV = 0;
        }
        else
        {
            while (cSV > 0 && vt[currfFirstVar].prePos - vt[listSV[cSV]].prePos <= maxSVdel)
                cSV--;
            
            
            
            if(cSV < noSV)
            {
                LoadUInt(ins_sv_arr+vt[listSV[cSV]].alt, delLen);
                while (vt[listSV[cSV]].prePos +  delLen  <= vt[currfFirstVar].prePos)
                {
                    cSV++;
                    if(cSV < noSV)
                        LoadUInt(ins_sv_arr+vt[listSV[cSV]].alt, delLen);
                    else
                        break;
                }
            }
            
            
            firstSV = cSV;
        }
        prevDelCount = 0;
        //check how many meaningful DELs before
        back = 1;
        if(firstVarNotZero)
        {
            while (currfFirstVar >= back && vt[currfFirstVar].prePos - vt[currfFirstVar - back].prePos <= maxDel)
            {
                if(vtType[currfFirstVar - back] == DEL)
                {
                    prevDelCount++;
                }
                back++;
            }
            
            back--;
            while(back)
            {
                if(vtType[currfFirstVar - back] == DEL)
                {
                    if(vt[currfFirstVar - back].prePos + vt[currfFirstVar - back].alt <= vt[currfFirstVar].prePos)
                        prevDelCount--;
                    else
                        break;
                }
                back--;
            }
        }
        /************************************/
        
        
        int32_t preVDelNo = prevDeletion[currfFirstVar];
        while (delFromRight[preVDelNo].nextRefPos > vt[currfFirstVar].prePos + 1) {
            --preVDelNo;
        }
        
        if(maxError>0)
            pos = checkVTconfiguration1(vt[currfFirstVar].prePos, vtCovered, currfFirstVar, currfFirstVar+vtCovered, lastVt, currfFirstVar,  nextRefPos, chCovered, read1, readLen, delta, preChToCheck, &result, pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck, wholeRead, wholeReadLen, leftPos, rightPos, error, maxError, nextRefPos, gapFirstVar, remainedToInsert, leftPos+1, vt[currfFirstVar].prePos, preVDelNo);
        else
              pos = checkVTconfigurationExact(vt[currfFirstVar].prePos, vtCovered, currfFirstVar, currfFirstVar+vtCovered, lastVt,currfFirstVar,  nextRefPos, chCovered, read1, readLen, delta, preChToCheck, &result, pos, max, prevDelCount, firstSV, lastSV, pos_vtToCheck, nextRefPos, gapFirstVar, remainedToInsert,  preVDelNo);

        
        if (pos == *max)
        {
            *max = *max * 2;
            
            
            tempResult = ( uint64_t *) realloc (result, *max * sizeof( uint64_t));
            
            if (NULL == tempResult)
            {
                printf("SA: Error in realloc, the message is - %s\n max = %llu\n", strerror(errno), *max);
                fflush(stdout);
                //free(result);
                exit(8);
            }
            else
            {
                result = tempResult;
            }
        }
        
        
        /************************************/
        /* Check if current not in del region that will keep the initial insertion */
        
        
        
        
        uint32_t dels;
        LoadUInt(ins_sv_arr+vt[currfFirstVar].alt, delLenFirst);
        dels = prevDelCount;
        back = 1;
        while(dels)
        {
            if(vtType[currfFirstVar - back] == DEL)
            {
                dels--;
                //if it is indeed in del region of other DEL
                if(vt[currfFirstVar - back].prePos + vt[currfFirstVar - back].alt > vt[currfFirstVar].prePos \
                   && vt[currfFirstVar - back].prePos + vt[currfFirstVar - back].alt < vt[currfFirstVar].prePos + delLenFirst)
                {
                    
                    pos = checkDelRegion(currfFirstVar - back, currfFirstVar, leftPos, gapFirstVar, error, maxError, wholeRead, pos_vtToCheck, vtCovered, firstSV, lastSV, prevDelCount, &result, pos, max, lastVt, nextRefPos, remainedToInsert, chCovered, read1, readLen, delta, preChToCheck, wholeReadLen, rightPos);
                    
                    
                    
                    
                }
            }
            back++;
        }
        
        
        if(firstSV < noSV)
            for(int32_t sv = lastSV; sv >= (int32_t)firstSV; sv--)
            {
                LoadUInt(ins_sv_arr+vt[listSV[sv]].alt, delLen);
                
                //if it is indeed in del region of other DEL
                if(vt[listSV[sv]].prePos +  delLen > vt[currfFirstVar].prePos \
                   && vt[listSV[sv]].prePos +  delLen < vt[currfFirstVar].prePos + delLenFirst)
                {
                    pos = checkDelRegion(listSV[sv], currfFirstVar, leftPos, gapFirstVar, error, maxError, wholeRead, pos_vtToCheck, vtCovered, firstSV, lastSV, prevDelCount, &result, pos, max, lastVt, nextRefPos, remainedToInsert, chCovered, read1, readLen, delta, preChToCheck, wholeReadLen, rightPos);
                    
                }
            }
        
        
    }
    
    
    
    
    //restore vtList
    for(int i = 0;i < vtListSize; i++)
        vtList[i] = lastVtList[i];
    free(lastVtList);
    
    *passed_result=result;
    
    return pos;
}











