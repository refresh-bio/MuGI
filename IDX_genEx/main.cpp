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
#include <fstream>


#include "coding.h"



uint64_t *resV;
uint32_t * vtList;
uint32_t max_vtToCheck, vtListSize;
uint32_t * vtToCheck;
uint64_t resSize;

char * altVt;
using namespace std;
unsigned char * ins_sv_arr;

#define NO_BITS_UINT 32U // (uint32_t) 8*sizeof(uint32_t)

unsigned char* reference;
uint32_t refSize;

unsigned char *s_ch;
uint32_t kl; //k-mer length
uint32_t n_sparse;

uint32_t  *VT_LUT;


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


uint32_t maxDel = 0, maxSVdel = 0;
uint32_t * listSV, noSV = 0;
uint32_t * prevSV, * prevDeletion;


typedef enum {fasta, fastq, simlple_list} reads_mode;



int32_t main(int32_t argc, const char * argv[])
{
    
    if(!argv[1] || !argv[2]|| !argv[3])
    {
        cout << "Required parameters:\n[index-file] [out-file] [no_reads] {[mode]}" << endl;
        cout << "where:\n[index-file] path to file with the index," << endl;
        cout << "[out-file] output file with example reads," << endl;
        cout << "[no_reads] number of reads to create, "<< endl;
        cout << "[mode] optional,format of the output \"fastq\" or \"fasta\" or \"list\" (default: \"list\")." << endl;
        exit(1);
    }

    
    
    
    reads_mode mode;
    if(argv[4])
    {
        if(strcmp(argv[4], "fasta") == 0 || strcmp(argv[4], "FASTA") == 0)
            mode = fasta;
        else if(strcmp(argv[4], "fastq") == 0 || strcmp(argv[4], "FASTQ") == 0)
            mode = fastq;
        else if(strcmp(argv[4], "list") == 0 || strcmp(argv[4], "LIST") == 0)
            mode = simlple_list;
        else
        {
            cout << "[mode] should be: \"fastq\" or \"fasta\" or \"list\" (default: \"list\")." << endl;
            exit(1);
        }

    }
    else
        mode = simlple_list;
    
    cout << "Reading data from: " << argv[1] << endl;
    
    FILE * in = fopen(argv[1], "rb");
    if (in == NULL) {
        printf("Cannot open %s\n", argv[1]);
        printf("The message is - %s\n", strerror(errno));
        exit(8);
    }
    
    uint32_t lut_encoded_chars, vt_lut_step  ; //for LUT_ENCODED_CHARS
    fread(&kl, 1, sizeof(uint32_t), in);
    fread(&n_sparse, 1, sizeof(uint32_t), in);
    fread(&lut_encoded_chars, 1, sizeof(uint32_t), in);
    
    uint32_t compRefSize;
    fread(&compRefSize,  sizeof(uint32_t), 1, in);
    reference = (uchar *)malloc((uint32_t)(compRefSize+1) * sizeof(uchar));
    fread(reference, 1, compRefSize, in);
    refSize = compRefSize * 2;
    if((reference[compRefSize-1] & 0xF) == 0xF)
        refSize--; //last empty -ok
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
    
    
    
    //temp -> VT_LUT_STEP
    fread(&vt_lut_step, sizeof(uint32_t), 1, in);
    VT_LUT = (uint32_t * ) malloc ((uint32_t)(refSize/VT_LUT_STEP) * sizeof(uint32_t));
    fread(VT_LUT, sizeof(uint32_t), refSize/VT_LUT_STEP, in);
    
    
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
    fread(uniq_var, sizeof(uint32_t), cnt_uniq_var * parts, in);
    fread(uniq_parts, sizeof(unsigned long long), cnt_uniq_parts * blocking, in);
    
     /********************* GENERATE FILE WITH EXAMPLE READS *********************************************/
    //cout << "The value of RAND_MAX is " <<  RAND_MAX << endl;
    //srand((unsigned)time(0));
    uint32_t random_prePos, random_individual, random_length, random_gap;
    
    // uint32_t lowest_pos = 0, highest_pos = refSize-1, lowest_ind = 0, highest_ind = inBitV-1, lowest_len = 75, highest_len = 100, lowest_gap = 0, highest_gap = 1;
    uint32_t lowest_pos = 0, highest_pos = refSize-1, lowest_ind = 0, highest_ind = inBitV-1, lowest_len = 120, highest_len = 170, lowest_gap = 0, highest_gap = 1, errToIntroduce, tempPos, tempCharNo;
    
    uint32_t range_pos=(highest_pos-lowest_pos)+1, range_ind = (highest_ind-lowest_ind)+1, range_len = (highest_len-lowest_len)+1, range_gap = (highest_gap-lowest_gap)+1;
    
    
    uint64_t a;
    char *tempRead;
    
    FILE * exReads;
    
    exReads = fopen(argv[2], "w");
    printf("File with example reads: %s\n", argv[2]);
    if (exReads == NULL) {
        printf("Cannot open %s\n", argv[2]);
        printf("The message is - %s\n", strerror(errno));
        exit(8);
    }
    
    int32_t noReads = 1000000;
    if(argv[3])
         noReads = atoi(argv[3]);

   // ofstream myfile;
   // myfile.open ("readsDesc.txt");

    char * quality = (char *)malloc((highest_len+1)*sizeof(char));
    memset(quality, 'I', highest_len);
    quality[highest_len] = '\0';
    
    printf("Preparing %d reads of length: %d to %d\nProbability of introducing change: 0.05 * len\n", noReads, lowest_len, highest_len);
    for(int32_t index = 0; index < noReads; index++)
    {
        //        a = rand();
        a = (uint64_t )range_pos*(uint64_t)rand();
        a = a/(uint64_t)(RAND_MAX + 1.0);
        random_prePos = (uint32_t)lowest_pos+(uint32_t)(a);
        
        a = (uint64_t )range_ind*(uint64_t)rand();
        a = a/(uint64_t)(RAND_MAX + 1.0);
        random_individual = (uint32_t)lowest_ind+(uint32_t)(a);
        
        a = (uint64_t )range_len*(uint64_t)rand();
        a = a/(uint64_t)(RAND_MAX + 1.0);
        random_length = (uint32_t)lowest_len+(uint32_t)(a);
        
        
        a = (uint64_t )(random_length*0.05)*(uint64_t)rand();
        errToIntroduce = (uint32_t)(a/(uint64_t)(RAND_MAX + 1.0));
        
        
        
        a = (uint64_t )range_gap*(uint64_t)rand();
        a = a/(uint64_t)(RAND_MAX + 1.0);
        random_gap= (uint32_t)lowest_gap+(uint32_t)(a);
        
     //   myfile << index << " " << random_prePos << " " << random_individual << " " << random_length << " " << random_gap << " " << endl;
        tempRead = decode(random_prePos, random_length, random_individual, 0xFFFFFFFF, 0);
        while(tempRead == NULL || strchr (tempRead,'N') || strchr(tempRead,'n') || strchr (tempRead,'R') || strchr(tempRead,'Y') || strchr (tempRead,'S')|| strchr (tempRead,'W') || strchr(tempRead,'K') || strchr(tempRead,'D') || strchr (tempRead,'M') || strchr(tempRead,'U'))
        {
            if(tempRead!=NULL)
                free(tempRead);
            a = (uint64_t )range_pos*(uint64_t)rand();
            a = a/(uint64_t)(RAND_MAX + 1.0);
            random_prePos = (uint32_t)lowest_pos+(uint32_t)(a);
            tempRead = decode(random_prePos, random_length, random_individual, 0xFFFFFFFF, 0);
                   //     myfile << random_prePos << endl;
        }
        
   //     myfile << errToIntroduce << endl;
        
        for(int32_t i = 0; i < errToIntroduce; i++)
        {
            
            a = (uint64_t )(random_length)*(uint64_t)rand();
            tempPos = (uint32_t)(a/(uint64_t)(RAND_MAX + 1.0));
   //         myfile << tempPos << endl;
            
            a = (uint64_t )(4)*(uint64_t)rand();
            tempCharNo = (uint32_t)(a/(uint64_t)(RAND_MAX + 1.0));
            
            switch(tempCharNo)
            {
                case 0:
                    tempRead[tempPos] = 'A';
                    break;
                case 1:
                    tempRead[tempPos] = 'C';
                    break;
                case 2:
                    tempRead[tempPos] = 'G';
                    break;
                case 3:
                    tempRead[tempPos] = 'T';
                    break;
                default:
                    cout << "wrong random search!\n";
            }
        }
        
        if(mode == fasta)
            fprintf(exReads, ">%d\n", index);
        else if(mode == fastq)
        {
            fprintf(exReads, "@%d\n", index);
        }
        
        fprintf(exReads, "%s\n", tempRead);
        
        if(mode == fastq)
        {
            quality[random_length] = '\0';
            fprintf(exReads, "+\n%s\n", quality);
            quality[random_length] = 'I';
        }

        
        if(strlen(tempRead) > 170)
            cout << tempRead << endl;
        
        
  //      myfile << tempRead << endl;
        free(tempRead);
        
    }
    /*****************************************************************************************************/

   // myfile.close();
    
    free(reference);
    free(vt);
    free(VT_LUT);
    free(vtType);
    free(ins_sv_arr);
    free(prevSV);
    free(listSV);
    free(prevDeletion);
    free(delFromRight);
    
    delete all_var ;
    delete uniq_var;
    delete uniq_parts;
    
    return 0;
}
























