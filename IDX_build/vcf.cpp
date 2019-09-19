/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */

#include "vcf.h"

VCF::VCF()
{
    file  = NULL;
    file_size = 0;
    no_gen = 0;
    no_var = 0;
}

bool VCF::Open(string file_name)
{
    char  data[200], ch;
    size_t tmpPt;
    
	if((file = my_fopen(file_name.c_str(), "r")) == NULL)
    {
        return false;
    }
    
    setvbuf(file, NULL, _IOFBF, 1 << 22);
	my_fseek(file, 0, SEEK_END);
	file_size = my_ftell(file);
	my_fseek(file, 0, SEEK_SET);
    

    //go through the initial header
    do {
        fscanf(file, "%s", data);
    }
    while (strncmp(data, "#CHROM", 6) != 0);
    fscanf(file, "%s %s %s %s %s %s %s %s", data, data, data, data, data, data, data, data);
    
    //count no genomes
    tmpPt = ftell(file);
    no_gen = 0;
    int32 c = 1;
    while ((c = getc(file)) != '\n' && c != '\r')
    {
        fscanf(file, "%s", data);
        no_gen++;
    }
    fseek(file, tmpPt, SEEK_SET);
    
    genomesNames = new string[no_gen*Data::ploidity];
    
    if(Data::ploidity == 1)
    {
        for (int32 i=0; i < no_gen; i++)
        {
            fscanf(file, "%s", data);
            genomesNames[i] = data;
            
        }
    }
    else if(Data::ploidity == 2)
    {
        for (int32 i=0; i < no_gen*2; i = i + 2)
        {
            fscanf(file, "%s", data);
            
            genomesNames[i] = data;
            genomesNames[i].append("-1");
            genomesNames[i+1] = data;
            genomesNames[i+1].append("-2");
            
        }
    }
    else
    {
        cout << "Data::ploidity " << Data::ploidity << " is not supported" << endl;
        return false;
    }
    
    ch=getc(file);
    tmpPt = ftell(file);
  
  
    while (EOF != (fscanf(file, "%*[^\n]"), fscanf(file, "%*c")))
        ++no_var;
   
    
    fseek(file, tmpPt, SEEK_SET);
    
	return true;
}


// ********************************************************************************************
bool VCF::Close()
{
	if(!file)
		return false;
    
	fclose(file);
	file = NULL;
	
	return true;
}

//********************************************************************************************


uint32_t VCF::get_no_var()
{
    return no_var;
}
//********************************************************************************************
uint32_t VCF::get_no_gen()
{
    return no_gen;
}
//********************************************************************************************

string * VCF::get_gen_names()
{
    return genomesNames;
}

//********************************************************************************************
/*
string  VCF::get_vl_file_name()
{
    return vl_file_name;
}*/

//********************************************************************************************
bool VCF::ProcessVCFData(VariantList * vl, BitVectors * bv)
{
    if(Data::ploidity == 2)
        ProcessVCFDataPloidity2(vl, bv);
    else if(Data::ploidity == 1)
        ProcessVCFDataPloidity1(vl, bv);
    else
    {
        cout << "Data::ploidity " << Data::ploidity << " is not supported" << endl;
        return false;
    }
    
    cout << "Processed 100% of the input vcf" << endl;
    
    bv->freeBitPos();
    return true;
    
}
//********************************************************************************************

bool VCF::ProcessVCFDataPloidity2(VariantList * vl, BitVectors * bv)
{
    
   /* FILE * vl_tmp_file;
    if( (vl_tmp_file = fopen(vl_file_name.c_str(), "w")) == NULL)
    {
        cout << "Error when opening temporary file!\n";
        return false;
    }*/
    
    uint32_t vt_count = 0, logTreshold = no_var/5;
    int32 noThreads;
    char **line;
    int32 chl = 1;
    int32  l, i, j;
    int32 writeno;
    char  data[200];
    uint32_t n;
    char type[10],  temp[1000000];
    
    
#pragma omp parallel
    {
        noThreads = omp_get_num_threads();
    }
    
    
    int32_t *bitsCalculated;
    bitsCalculated = (int32_t *) malloc(noThreads * sizeof(int32_t));
    if(bitsCalculated == NULL)
    {
        fprintf(stderr, "Out of memory\n");
        exit(9);
    }

    
    unsigned char *** tempBits1, *** tempBits2;
    int line_max[NO_LINES_PER_THREAD*noThreads];
    
    int char_max[noThreads];
    
    line = (char** ) malloc(NO_LINES_PER_THREAD*noThreads*sizeof(char*));
    if(line == NULL)
    {
        fprintf(stderr, "Out of memory\n");
        exit(9);
    }
    for (i=0; i < NO_LINES_PER_THREAD*noThreads; i++)
    {
        line[i] = (char *) malloc (INITIAL_VT_LINE_WIDTH * sizeof(char));
        if(line[i] == NULL)
        {
            fprintf(stderr, "Out of memory\n");
            exit(9);
        }
        line_max[i] = INITIAL_VT_LINE_WIDTH;
    }
    
    tempBits1 = (unsigned char ***) malloc(noThreads * sizeof(unsigned char**));
    if(tempBits1 == NULL)
    {
        fprintf(stderr, "Out of memory\n");
        exit(9);
    }
    tempBits2 = (unsigned char ***) malloc(noThreads * sizeof(unsigned char**));
    if(tempBits2 == NULL)
    {
        fprintf(stderr, "Out of memory\n");
        exit(9);
    }
    for (i = 0; i < noThreads; i++)
    {
        tempBits1[i] = (unsigned char**)malloc(no_gen * sizeof( unsigned char*));
        if(tempBits1[i] == NULL)
        {
            fprintf(stderr, "Out of memory\n");
            exit(9);
        }
        tempBits2[i] = (unsigned char**)malloc(no_gen * sizeof(unsigned  char*));
        if(tempBits2[i] == NULL)
        {
            fprintf(stderr, "Out of memory\n");
            exit(9);
        }
        
        char_max[i] = 1+(NO_LINES_PER_THREAD-1)/8;
        
        for (j = 0; j < no_gen; j++)
        {
            tempBits1[i][j] = (unsigned char*)malloc(((NO_LINES_PER_THREAD-1)/8+1) * sizeof(unsigned char));
            if(tempBits1[i][j] == NULL)
            {
                fprintf(stderr, "Out of memory\n");
                exit(9);
            }
            tempBits2[i][j] = (unsigned char*)malloc(((NO_LINES_PER_THREAD-1)/8+1) * sizeof(unsigned char));
            if(tempBits2[i][j] == NULL)
            {
                fprintf(stderr, "Out of memory\n");
                exit(9);
            }
            
            
        }
        
    }

    
    int * VLbufPos = (int*) malloc (noThreads * sizeof (int));
    char ** VLbuffer;
    VLbuffer = (char **) malloc (noThreads *sizeof(char*));
    int * buf_max;
    buf_max = (int *) malloc (noThreads * sizeof(int));
    

    
    
    char **chrom, **name, **refC, **altList, **alt, **filter, **info, **format, **qual, **genotype, **endPos;
    
    
    chrom = (char **) malloc (noThreads *sizeof(char*));
    name = (char **) malloc (noThreads *sizeof(char*));
    refC = (char **) malloc (noThreads *sizeof(char*));
    altList = (char **) malloc (noThreads *sizeof(char*));
    alt = (char **) malloc (noThreads *sizeof(char*));
    filter = (char **) malloc (noThreads *sizeof(char*));
    info = (char **) malloc (noThreads *sizeof(char*));
    format = (char **) malloc (noThreads *sizeof(char*));
    qual = (char **) malloc (noThreads *sizeof(char*));
    genotype = (char **) malloc (noThreads *sizeof(char*));
    endPos = (char **) malloc (noThreads *sizeof(char*));
    
    for (i = 0; i < noThreads; i++)
    {
        VLbuffer[i] = (char *) malloc (NO_LINES_PER_THREAD*100*sizeof(char));
        buf_max[i] = NO_LINES_PER_THREAD*100;
        
        
        chrom[i] = (char *) malloc (15 * sizeof(char));
        name[i] = (char *) malloc (1000 *sizeof(char));
        refC[i] = (char *) malloc (1000000 *sizeof(char));
        altList[i] = (char *) malloc (1000000 *sizeof(char));
        alt[i] = (char *) malloc (1000000 *sizeof(char));
        filter[i] = (char *) malloc (20 *sizeof(char));
        info[i] = (char *) malloc (1000 *sizeof(char));
        format[i] = (char *) malloc (100 *sizeof(char));
        qual[i] = (char *) malloc (20 *sizeof(char));
        genotype[i] = (char *) malloc (200 *sizeof(char));
        endPos[i] = (char *) malloc (25 *sizeof(char));
    }

    
    //chl = getc(file); //gets '\n'
    
    int move, lastLine;
    char *ptr_s;
    
    while (chl != EOF)
    {
        
        if(vt_count > logTreshold)
        {
            cout << "Processed " << (((int)(((float)(logTreshold)/(float)(no_var))*100))+1)/20*20 << "% of the input vcf" << endl;
            logTreshold = logTreshold + no_var/5;
        }
        
        
        
        lastLine = NO_LINES_PER_THREAD * noThreads;
        for (i = 0; i < lastLine; i++)
        {
            
            
            move = 0;
            while(1)
            {
                
                if (fgets(line[i]+move, line_max[i]-move, file)  && !feof(file) ) {
                    ptr_s = strchr(line[i], '\n');
                    if (ptr_s) {
                        *ptr_s = '\0';
                        break;
                    }
                    else
                    {
                        line[i] = (char *) realloc(line[i], (line_max[i]*2)* sizeof(char));
                        move = line_max[i] - 1;
                        line_max[i] = line_max[i]*2;
                    }
                }
                else
                {
                    line[i][0] = '\n';
                    lastLine = i - 1;
                    break;
                }
            };
        }
        
#pragma omp parallel
        {
#pragma omp for
            for(int t = 0; t < noThreads; t++)
            {
                uint64_t pos,  currentPos = 0, currChar, offset;
                uint32_t  gt_pos, noAlt, a, gt1, gt2;
                
                char* pp, * rr;
                pp = new char[200];
                rr = new char[10];
                
                short currentLine;
                uint32_t mv, ptr;
                int j;
                long int shift;
                bool variantOK;
                char *strtokState1, *strtokState2, *strtokState3;
                
                int returnValue;
                bool isInsertion;
                char * chr_ptr;
                uint64_t temp_int, delLen;
                
                currentPos = 0;
                bitsCalculated[t] = 0;
                
                for (int g = 0; g < no_gen; g++)
                {
                    for (int b = 0; b < char_max[t]; b++)
                    {
                        tempBits1[t][g][b] = 0;
                        tempBits2[t][g][b] = 0;
                    }
                }
                VLbufPos[t] = 0;
                
                for (j = 0; j < NO_LINES_PER_THREAD; j++)
                {
                    mv = 0, ptr = 0;
                    variantOK = true;
                    currentLine = NO_LINES_PER_THREAD*t + j;
                    if(currentLine > lastLine)
                    {
                        break;
                    }
                    sscanf(line[currentLine], "%s %lld %s %s %s %s %s %s %s%n", chrom[t], &pos, name[t], refC[t], altList[t], qual[t], filter[t], info[t], format[t], &ptr);
                    if(ptr == 0)
                        break;
                    
                    
                    
                    //adjust positions according to input fasta file (reference sequence)
                    offset = 0;
                    for(int f = 0; f < Data::fasta->sequences.size(); f++)
                    {
                        if(strcmp(chrom[t], Data::fasta->sequences[f].chrom.c_str()) == 0) //variant from f chromosome
                        {
                            //cout << Data::fasta->sequences[f].start_pos << " " << pos << endl;
                            //pos = pos + Data::fasta->sequences[f].start_pos;
                            offset = Data::fasta->sequences[f].start_pos;
                            //cout << Data::fasta->sequences[f].start_pos << " " << pos << endl;
                            break;
                        }
                    }
                    
                    
                    
                    
                    
                    
                    
                    pp = strtok_r(format[t], ":", &strtokState1);
                    gt_pos = 0;
                    while (strcmp(pp,"GT") != 0) {
                        pp = strtok_r(NULL, ":", &strtokState1);
                        gt_pos++;
                    }
                    
                    currChar = 0;
                    noAlt = 1;
                    while(altList[t][currChar] != '\0')
                    {
                        if(altList[t][currChar] == ',')
                            noAlt++;
                        currChar++;
                    };
                    
                    
                    
                    currChar = 0;
                    while (altList[t][currChar] != '\0')
                    {
                        a = 0;
                        
                        while (altList[t][currChar] != ',' && altList[t][currChar] != '\0')
                        {
                            alt[t][a] = altList[t][currChar];
                            a++;
                            currChar++;
                        }
                        alt[t][a] = '\0';
                        if (altList[t][currChar] == ',')
                            currChar++;
                        
                        ///*** just to check if all variants are OK
                        if (strcmp("PASS", filter[t]) != 0 && strcmp(".", filter[t]) != 0)
                        {
                            printf("*** WARNING ***\n");
                            printf("\tPosition: %llu. Variant's FILTER value is neither 'PASS' nor '.': %s\n", pos+offset, filter[t]);
                            printf("\tVariant will be processed anyway. If it should not, it should be removed from the VCF file\n");
                            
                        }
                        
                        if (strstr(info[t], "VT=SNP") != NULL || strcmp(chrom[t], "Y") == 0)  //VT is a SNP or chrY(only SNPs)
                        {
                            if (strlen(refC[t]) != 1)
                            {
                                printf("*** UNSUPPORTED CASE, variant will be omitted. ***\n");
                                printf("\tPosition: %llu. It is a SNP and REF length is not 1.\n", pos+offset);
                                variantOK = false;
                                break;
                            }
                            if (strlen(alt[t]) != 1)
                            {
                                printf("*** UNSUPPORTED CASE, variant will be omitted. ***\n");
                                printf("\tPosition: %llu. It is a SNP and ALT length is not 1\n", pos+offset);
                                
                                variantOK = false;
                                break;
                            }
                            
                            
                            
                            while(1)
                            {
                                returnValue = snprintf(VLbuffer[t] + VLbufPos[t], buf_max[t] - VLbufPos[t], "%llu SNP %s\n", pos+offset, alt[t]);
                                if(returnValue > 0 && returnValue < buf_max[t] - VLbufPos[t])
                                {
                                    VLbufPos[t] +=  returnValue;
                                    break;
                                }
                                else
                                {
                                    VLbuffer[t] = (char *) realloc(VLbuffer[t], (buf_max[t]*2)* sizeof(char));
                                    buf_max[t] = buf_max[t]*2;
                                }
                            }
                        }
                        else if (strstr(info[t], "VT=INDEL") != NULL)  //VT is an INDEL
                        {
                            if ((255&alt[t][0]) != (255&refC[t][0]))
                            {
                                printf("*** UNSUPPORTED CASE, variant will be omitted. ***\n");
                                printf("\tPosition: %llu. It is an INDEL and first characters of ALT and REF do not match.\n", pos+offset);
                                variantOK = false;
                                break;
                            }
                            
                            if(strlen(refC[t])==1) //insertion
                            {
                                while(1)
                                {
                                    returnValue = snprintf(VLbuffer[t] + VLbufPos[t], buf_max[t] - VLbufPos[t], "%llu INS %s\n", pos+offset, alt[t]+1);
                                    if(returnValue > 0 && returnValue < buf_max[t] - VLbufPos[t])
                                    {
                                        VLbufPos[t] +=  returnValue;
                                        break;
                                    }
                                    else
                                    {
                                        VLbuffer[t] = (char *) realloc(VLbuffer[t], (buf_max[t]*2)* sizeof(char));
                                        buf_max[t] = buf_max[t]*2;
                                    }
                                }
                                
                            }
                            else //if (strlen(alt[t])==1) //deletion
                            {
                                
                                while(1)
                                {
                                    returnValue = snprintf(VLbuffer[t] + VLbufPos[t], buf_max[t] - VLbufPos[t],  "%lld DEL %ld\n", pos+offset, strlen(refC[t])-1);
                                    if(returnValue > 0 && returnValue < buf_max[t] - VLbufPos[t])
                                    {
                                        VLbufPos[t] +=  returnValue;
                                        break;
                                    }
                                    else
                                    {
                                        VLbuffer[t] = (char *) realloc(VLbuffer[t], (buf_max[t]*2)* sizeof(char));
                                        buf_max[t] = buf_max[t]*2;
                                    }
                                }
                                
                            }
                            
                        }
                        else if (strstr(info[t], "VT=SV") != NULL)  //VT is a SV
                        {
                            if (strstr(info[t], "SVTYPE=DEL") != NULL)
                            {
                                endPos[t][0] = '\0';
                                if (strstr(info[t], ";END=") != NULL)
                                {
                                    chr_ptr = strstr(info[t], ";END=");
                                    temp_int = 5;
                                    while (chr_ptr[temp_int] != ';' &&  chr_ptr[temp_int] != '\0')
                                    {
                                        endPos[t][temp_int-5] =  chr_ptr[temp_int];
                                        temp_int++;
                                    }
                                    endPos[t][temp_int-5]='\0';
                                }
                                else if (strncmp(info[t], "END=", 4) == 0)
                                {
                                    temp_int=4;
                                    while (info[t][temp_int] != ';' &&  info[t][temp_int] != '\0')
                                    {
                                        endPos[t][temp_int-4] =  info[t][temp_int];
                                        temp_int++;
                                    }
                                    endPos[j-4][0]='\0';
                                }
                                else
                                {
                                    printf("*** UNSUPPORTED CASE, variant will be omitted. ***\n");
                                    printf("\tPosition: %lld. SVTYPE=DEL and no END in info: %s\n", pos+offset, info[t]);
                                    
                                    variantOK = false;
                                    break;
                                }
                                delLen = (uint64_t)(atol(endPos[t])-pos);
                                
                                if (strncmp(alt[t], "<DEL>", 5) == 0) {
                                    isInsertion=false;
                                }
                                else if (strlen(alt[t]) == 1 && alt[t][0]==refC[t][0])
                                {
                                    isInsertion=false;
                                }
                                else if (strlen(alt[t]) > 1 && alt[t][0]==refC[t][0])
                                {
                                    isInsertion=true;
                                }
                                else {
                                    printf("*** UNSUPPORTED CASE, variant will be omitted. ***\n");
                                    printf("\tPosition: %lld. SVTYPE=DEL and for given REF, wrong ALT data: %s\n", pos+offset, alt[t]);
                                    
                                    variantOK = false;
                                    break;
                                }
                                
                                if (isInsertion)
                                {
                                    
                                    while(1)
                                    {
                                        returnValue = snprintf(VLbuffer[t] + VLbufPos[t], buf_max[t] - VLbufPos[t],  "%lld SV %llu %s\n", pos+offset, delLen, alt[t]+1);
                                        if(returnValue > 0 && returnValue < buf_max[t] - VLbufPos[t])
                                        {
                                            VLbufPos[t] +=  returnValue;
                                            break;
                                        }
                                        else
                                        {
                                            VLbuffer[t] = (char *) realloc(VLbuffer[t], (buf_max[t]*2)* sizeof(char));
                                            buf_max[t] = buf_max[t]*2;
                                        }
                                    }
                                    
                                }
                                else
                                {
                                    
                                    while(1)
                                    {
                                        returnValue = snprintf(VLbuffer[t] + VLbufPos[t], buf_max[t] - VLbufPos[t],  "%lld SV %llu -\n", pos+offset, delLen);
                                        if(returnValue > 0 && returnValue < buf_max[t] - VLbufPos[t])
                                        {
                                            VLbufPos[t] +=  returnValue;
                                            break;
                                        }
                                        else
                                        {
                                            VLbuffer[t] = (char *) realloc(VLbuffer[t], (buf_max[t]*2)* sizeof(char));
                                            buf_max[t] = buf_max[t]*2;
                                        }
                                    }
                                    
                                    
                                }
                                
                                
                            }
                            else
                            {
                                printf("*** UNSUPPORTED CASE, variant will be omitted. ***\n");
                                printf("\tUnsupported variant VT=SV at position: %lld  and not SVTYPE=DEL, INFO field:  %s\n",  pos+offset, info[t]);
                                
                                variantOK = false;
                                break;
                            }
                            
                        }
                        else
                        {
                            printf("*** UNSUPPORTED CASE, variant will be omitted. ***\n");
                            printf("\tUnsupported variant at position: %lld  (neither SNP, INDEL nor SV), INFO field: %s \n", pos+offset, info[t]);
                            
                            variantOK = false;
                            break;
                        }
                        
                        
                    }
                    
                    if(!variantOK)
                    {
                        no_var--;
                        continue;
                    }
                    // if there is a change in number of bitsCalculated that requires another unsigned char for storage, reallocate
                    if ((bitsCalculated[t]+noAlt) > char_max[t]*8)
                    {
                        for (int g = 0; g < no_gen; g++)
                        {
                            
                            tempBits1[t][g] = (unsigned char*) realloc(tempBits1[t][g], (char_max[t] + 1)* sizeof(unsigned char));
                            tempBits2[t][g] = (unsigned char*) realloc(tempBits2[t][g], (char_max[t] + 1)* sizeof(unsigned char));
                            tempBits1[t][g][char_max[t]] = 0;
                            tempBits2[t][g][char_max[t]] = 0;
                            
                        }
                        char_max[t] = char_max[t] + 1;
                        
                    }
                    
                    
                    
                    
                    for (int g = 0; g < no_gen; g++)
                    {
                        
                        sscanf(line[currentLine]+ptr, "%s%n", genotype[t], &mv);
                        ptr = ptr + mv;
                        
                        
                        pp = strtok_r(genotype[t], ":", &strtokState2);
                        for(int m=0; m < gt_pos; m++)
                            pp = strtok_r( NULL, ":", &strtokState2);
                        
                        
                        //get first/left chromosome genotype
                        rr = strtok_r(pp, "|" , &strtokState3);
                        if(rr == NULL)
                            rr = pp;
                        if (strcmp(rr,".") == 0)
                            rr[0] = '0';
                        gt1 = atoi(rr);
                        
                        if (gt1 > noAlt)
                        {
                            printf("*** WARNING***\n");
                            printf("\tPosition: %lld, GT data for sample no %d is wrong: %s, no variant will be placed in the consensus sequence here\n", pos+offset, i+1, data);
                            gt1 = 0;
                        }
                        //get second/right chromosome genotype
                        rr = strtok_r(NULL, "|" , &strtokState3);
                        if(rr == NULL)
                            rr = pp;
                        if (strcmp(rr,".") == 0)
                            rr[0] = '0';
                        gt2 = atoi(rr);
                        if (gt2 > noAlt)
                        {
                            printf("*** WARNING***\n");
                            printf("\tPosition: %lld, GT data for sample no %d is wrong: %s, no variant will be placed in the consensus sequence here\n", pos+offset, i+1, data);
                            gt2 = 0;
                        }
                        
                        
                        
                        for (a = 0; a < noAlt; a++)
                        {
                            if(gt1 == a+1)
                            {
                                shift = (8*sizeof(unsigned char)-1-(bitsCalculated[t]+a)%(8*sizeof(unsigned char)));
                                tempBits1[t][g][(bitsCalculated[t]+a)/(8*sizeof(unsigned char))] = tempBits1[t][g][(bitsCalculated[t]+a)/(8*sizeof(unsigned char))]|((unsigned char)1<<shift);
                                
                            }
                            
                            if(gt2 == a+1)
                            {
                                shift = (8*sizeof(unsigned char)-1-(bitsCalculated[t]+a)%(8*sizeof(unsigned char)));
                                tempBits2[t][g][(bitsCalculated[t]+a)/(8*sizeof(unsigned char))] = tempBits2[t][g][(bitsCalculated[t]+a)/(8*sizeof(unsigned char))]|((unsigned char)1<<shift);
                            }
                        }
                        
                    }
                    bitsCalculated[t] = bitsCalculated[t] + noAlt;
                }
            }
        }
        
        
        
        
        for (int t = 0; t < noThreads; t++)
        {
            for (int g = 0; g < no_gen; g++)
            {
                for (int b = 0; b <= (bitsCalculated[t]-1)/8; b++)
                {
                    
                    writeno = 8;
                    if(bitsCalculated[t]<8)
                        writeno = bitsCalculated[t];
                    if((b == (bitsCalculated[t]-1)/8 ) && bitsCalculated[t]%8 != 0)
                        writeno = bitsCalculated[t]%8;
                    
                    bv->writeBits(g*2, tempBits1[t][g][b], writeno);
                    //writeBits(genomesNames[g*2], bv.bit_buffer[g*2], bufPos1, tempBits1[t][g][b], writeno);
                    bv->writeBits(g*2+1, tempBits2[t][g][b], writeno);
                    //writeBits(genomesNames[g*2+1], bv.bit_buffer[g*2+1], bufPos2, tempBits2[t][g][b], writeno);
                }
            }
            
            l = 0;
            
            for (int v = 0; v < bitsCalculated[t]; v++)
            {
                /* fprintf(vl_tmp_file, "%d ", vt_count + v);
                 while(VLbuffer[t][l] != '\n')
                 {
                 putc(VLbuffer[t][l++], vl_tmp_file);
                 
                 }
                 putc(VLbuffer[t][l++], vl_tmp_file);
                 
                 */
                {
                    
                    sscanf(VLbuffer[t]+l, "%u %s %n",  &vl->vt[vt_count].prePos, type, &n);
                    l += n;
                    
                    vl->vt[vt_count].prePos = vl->vt[vt_count].prePos - 1; //go from 1-base to 0-base
                    
                    vl->prevSV[vt_count] = vl->no_sv - 1;
                    if (strstr(type, "SNP") != NULL)  //VT is a SNP
                    {
                        vl->vtType[vt_count] = SNP;
                        vl->vtAux[vt_count].alt = new char[1];
                        sscanf(VLbuffer[t]+l, "%s\n%n", vl->vtAux[vt_count].alt, &n);
                        l += n;
                        vl->vtAux[vt_count].delLen = 0;
                        vl->vt[vt_count].prePos = vl->vt[vt_count].prePos - 1;  //so it points to previous position (one before SNP) as for other variants
                    }
                    else if (strstr(type, "INS") != NULL)  //VT is an INS
                    {
                        vl->vtType[vt_count] = INS;
                        sscanf(VLbuffer[t]+l, "%s\n%n", temp, &n);
                        l += n;
                        vl->vtAux[vt_count].alt = new char[strlen(temp)+1];
                        vl->no_chars_inserted += strlen(temp);
                        memcpy(vl->vtAux[vt_count].alt , temp, strlen(temp));
                        vl->vtAux[vt_count].alt[strlen(temp)]='\0';
                        //cout << "ins: " << i << " " << strlen(temp) <<  " " <<  temp << " " << vtAux[vt_count].alt << endl;
                        vl->vtAux[vt_count].delLen = 0;
                        vl->no_ins++;
                        
                    }
                    else if (strstr(type, "DEL") != NULL)  //VT is an DEL
                    {
                        vl->vtType[vt_count] = DEL;
                        sscanf(VLbuffer[t]+l, "%u\n%n", &vl->vtAux[vt_count].delLen, &n);
                        l += n;
                        vl->vtAux[vt_count].alt = new char[1];
                        vl->vtAux[vt_count].alt[0] = '-';
                        
                        if(vl->vtAux[vt_count].delLen > vl->max_del)
                            vl->max_del = vl->vtAux[vt_count].delLen;
                        
                        vl->no_del++;
                    }
                    else if (strstr(type, "SV") != NULL)  //VT is an SV
                    {
                        vl->vtType[vt_count] = SV;
                        sscanf(VLbuffer[t]+l, "%u %s\n%n", &vl->vtAux[vt_count].delLen, temp, &n);
                        l += n;
                        vl->vtAux[vt_count].alt = new char[strlen(temp)+1];
                        memcpy(vl->vtAux[vt_count].alt , temp, strlen(temp));
                        vl->vtAux[vt_count].alt[strlen(temp)]='\0';
                        
                        
                        if(vl->vtAux[vt_count].alt[0] != '-')
                            vl->no_chars_inserted += strlen(temp);
                        
                        // cout << "ins: " << i << " " << strlen(temp) <<  " " <<  temp << " " << vtAux[vt_count].alt << endl;
                        vl->no_sv++;
                        
                        if(vl->vtAux[vt_count].delLen > vl->max_sv_del)
                            vl->max_sv_del = vl->vtAux[vt_count].delLen;
                        
                    }
                }
                
                vt_count++;
            }
            
            //vt_count += bitsCalculated[t];
        }
        
        
        if (chl != EOF)
            chl = getc(file);
        fseek( file, -1, SEEK_CUR);
    }
    
   
    for (i = 0; i < NO_LINES_PER_THREAD * noThreads; i++)
        free(line[i]);
    free(line);

    
    for (i = 0; i < noThreads; i++)
    {
        for (j = 0; j < no_gen; j++)
        {
            free(tempBits1[i][j]);
            free(tempBits2[i][j]);
        }
        
        free(tempBits1[i]);
        free(tempBits2[i]);
        free(VLbuffer[i]);
        free(chrom[i]);
        free(name[i]);
        free(refC[i]);
        free(altList[i]);
        free(alt[i]);
        free(filter[i]);
        free(info[i]);
        free(format[i]);
        free(qual[i]);
        free(genotype[i]);
        free(endPos[i]);
        
       
    }
    
    free(tempBits1);
    free(tempBits2);
    free(VLbuffer);
    free(VLbufPos);
    free(chrom);
    free(name);
    free(refC);
    free(altList);
    free(alt);
    free(filter);
    free(info);
    free(format);
    free(qual);
    free(genotype);
    free(endPos);
    free(bitsCalculated);
    free(buf_max);
    
    //fclose(vl_tmp_file);
    
    return true;
}

// *******************************************************************************************

// *******************************************************************************************

bool VCF::ProcessVCFDataPloidity1(VariantList * vl, BitVectors * bv)
{
    
    
    char * tmpChrPtr;
    unsigned char * tmpUChrPtr;
    //int32_t * tmpIntPtr;
    /*FILE * vl_tmp_file;
     if( (vl_tmp_file = fopen(vl_file_name.c_str(), "w")) == NULL)
     {
     cout << "Error when opening temporary file!\n";
     return false;
     }*/
    
    uint32_t vt_count = 0, logTreshold = no_var/5;
    int32 noThreads;
    char **line,  chl = 1;
    int32  l, i, j;
    int32 writeno;
    
    uint32_t n;
    char type[10],  temp[1000000];
    
#pragma omp parallel
    {
        noThreads = omp_get_num_threads();
    }
    
    uint32_t *bitsCalculated;
    bitsCalculated = (uint32_t *) malloc(noThreads * sizeof(uint32_t));
    if(bitsCalculated == NULL)
    {
        fprintf(stderr, "Out of memory\n");
        exit(9);
    }
    
    unsigned char *** tempBits;
    int line_max[NO_LINES_PER_THREAD*noThreads];
    
    int char_max[noThreads];
    
    line = (char** ) malloc(NO_LINES_PER_THREAD*noThreads*sizeof(char*));
    if(line == NULL)
    {
        fprintf(stderr, "Out of memory\n");
        exit(9);
    }
    for (i=0; i < NO_LINES_PER_THREAD*noThreads; i++)
    {
        line[i] = (char *) malloc (INITIAL_VT_LINE_WIDTH * sizeof(char));
        if(line[i] == NULL)
        {
            fprintf(stderr, "Out of memory\n");
            exit(9);
        }
        line_max[i] = INITIAL_VT_LINE_WIDTH;
    }
    
    tempBits = ( unsigned char ***) malloc(noThreads * sizeof(unsigned char**));
    if(tempBits == NULL)
    {
        fprintf(stderr, "Out of memory\n");
        exit(9);
    }
    
    for (i = 0; i < noThreads; i++)
    {
        tempBits[i] = (unsigned char**)malloc(no_gen * sizeof(unsigned char*));
        if(tempBits[i] == NULL)
        {
            fprintf(stderr, "Out of memory\n");
            exit(9);
        }
        
        
        char_max[i] = 1+(NO_LINES_PER_THREAD-1)/8;
        
        for (j = 0; j < no_gen; j++)
        {
            tempBits[i][j] = (unsigned char*)malloc(char_max[i]  * sizeof( unsigned char));
            if(tempBits[i][j] == NULL)
            {
                fprintf(stderr, "Out of memory\n");
                exit(9);
            }
        }
        
    }
    
    int * VLbufPos = (int*) malloc (noThreads * sizeof (int));
    char ** VLbuffer;
    VLbuffer = (char **) malloc (noThreads *sizeof(char*));
    int * buf_max;
    buf_max = (int *) malloc (noThreads * sizeof(int));
    
    
    char **chrom, **name, **refC, **altList, **alt, **filter, **info, **format, **qual, **genotype, **endPos;
    
    
    chrom = (char **) malloc (noThreads *sizeof(char*));
    name = (char **) malloc (noThreads *sizeof(char*));
    refC = (char **) malloc (noThreads *sizeof(char*));
    altList = (char **) malloc (noThreads *sizeof(char*));
    alt = (char **) malloc (noThreads *sizeof(char*));
    filter = (char **) malloc (noThreads *sizeof(char*));
    info = (char **) malloc (noThreads *sizeof(char*));
    format = (char **) malloc (noThreads *sizeof(char*));
    qual = (char **) malloc (noThreads *sizeof(char*));
    genotype = (char **) malloc (noThreads *sizeof(char*));
    endPos = (char **) malloc (noThreads *sizeof(char*));
    
    
    for (i = 0; i < noThreads; i++)
    {
        VLbuffer[i] = (char *) malloc (NO_LINES_PER_THREAD*100*sizeof(char));
        buf_max[i] = NO_LINES_PER_THREAD*100;
        
        
        chrom[i] = (char *) malloc (15 * sizeof(char));
        name[i] = (char *) malloc (1000 *sizeof(char));
        refC[i] = (char *) malloc (1000000 *sizeof(char));
        altList[i] = (char *) malloc (1000000 *sizeof(char));
        alt[i] = (char *) malloc (1000000  *sizeof(char));
        filter[i] = (char *) malloc (20 *sizeof(char));
        info[i] = (char *) malloc (1000 *sizeof(char));
        format[i] = (char *) malloc (100 *sizeof(char));
        qual[i] = (char *) malloc (20 *sizeof(char));
        genotype[i] = (char *) malloc (200 *sizeof(char));
        endPos[i] = (char *) malloc (25 *sizeof(char));
    }
    
    //chl = getc(file); //gets '\n'
    
    int move, lastLine;
    char *ptr_s;
    
    
    
    
    while (chl != EOF)
    {
        if(vt_count > logTreshold)
        {
            cout << "Processed " << (((int)(((float)(logTreshold)/(float)(no_var))*100))+1)/20*20 << "% of the input vcf" << endl;
            logTreshold = logTreshold + no_var/5;
        }
        
        
        
        lastLine = NO_LINES_PER_THREAD * noThreads;
        for (i = 0; i < lastLine; i++)
        {
            
            
            move = 0;
            while(1)
            {
                if (fgets(line[i]+move, line_max[i]-move, file)  && !feof(file) ) {
                    ptr_s = strchr(line[i], '\n');
                    if (ptr_s) {
                        *ptr_s = '\0';
                        break;
                    }
                    else
                    {
                        
                        tmpChrPtr = (char *) realloc(line[i], (line_max[i]*2)* sizeof(char));
                        if (NULL ==tmpChrPtr) {
                            printf("Error in realloc, the message is - %s\n", strerror(errno));
                            exit(8);
                        }
                        else
                        {
                            line[i] = tmpChrPtr;
                            move = line_max[i] - 1;
                            line_max[i] = line_max[i]*2;
                        }
                        
                        
                        
                    }
                }
                else
                {
                    line[i][0] = '\n';
                    lastLine = i - 1;
                    break;
                }
            };
            
        }
        
        //#pragma omp parallel
        {
            //#pragma omp for
            for(int t = 0; t < noThreads; t++)
            {
                uint64_t pos,  currentPos = 0, currChar, offset;
                uint32_t  gt_pos, noAlt, a, gt;
                
                char* pp, * rr;
                pp = new char[200];
                rr = new char[10];
                
                short currentLine;
                uint32 mv, ptr;
                int64_t j;
                int64_t shift;
                bool variantOK;
                char *strtokState1, *strtokState2;
                
                int64_t returnValue;
                bool isInsertion;
                char * chr_ptr;
                uint64_t temp_int, delLen;
                
                currentPos = 0;
                bitsCalculated[t] = 0;
                
                for (int g = 0; g < no_gen; g++)
                {
                    for (int b = 0; b < char_max[t]; b++)
                    {
                        tempBits[t][g][b] = 0x00000000;
                    }
                }
                VLbufPos[t] = 0;
                
                for (j = 0; j < NO_LINES_PER_THREAD; j++)
                {
                    mv = 0, ptr = 0;
                    variantOK = true;
                    currentLine = NO_LINES_PER_THREAD*t + j;
                    if(currentLine > lastLine)
                    {
                        break;
                    }
                    sscanf(line[currentLine], "%s %lld %s %s %s %s %s %s %s%n", chrom[t], &pos, name[t], refC[t], altList[t], qual[t], filter[t], info[t], format[t], &ptr);
                    
                    if(ptr == 0)
                        break;
                    
                    //adjust positions according to input fasta file (reference sequence)
                    offset = 0;
                    for(int f = 0; f < Data::fasta->sequences.size(); f++)
                    {
                        if(strcmp(chrom[t], Data::fasta->sequences[f].chrom.c_str()) == 0) //variant from f chromosome
                        {
                            offset =  Data::fasta->sequences[f].start_pos;
                            break;
                        }
                    }
                    
                    pp = strtok_r(format[t], ":", &strtokState1);
                    gt_pos = 0;
                    while (strcmp(pp,"GT") != 0) {
                        pp = strtok_r(NULL, ":", &strtokState1);
                        gt_pos++;
                    }
                    
                    currChar = 0;
                    noAlt = 1;
                    while(altList[t][currChar] != '\0')
                    {
                        if(altList[t][currChar] == ',')
                            noAlt++;
                        currChar++;
                    };
                    
                    
                    // newVTno = 0;
                    currChar = 0;
                    while (altList[t][currChar] != '\0')
                    {
                        a = 0;
                        
                        while (altList[t][currChar] != ',' && altList[t][currChar] != '\0')
                        {
                            alt[t][a] = altList[t][currChar];
                            a++;
                            currChar++;
                        }
                        alt[t][a] = '\0';
                        if (altList[t][currChar] == ',')
                            currChar++;
                        
                        ///*** just to check if all variants are OK
                        if (strcmp("PASS", filter[t]) != 0 && strcmp(".", filter[t]) != 0)
                        {
                            printf("*** WARNING ***\n");
                            printf("\tPosition: %lld. Variant's FILTER value is neither 'PASS' nor '.': %s\n", pos+offset, filter[t]);
                            printf("\tVariant will be processed anyway. If it should not, it should be removed from the VCF file\n");
                            
                        }
                        
                        if ((strstr(info[t], "VT=SNP") != NULL) || ( strstr(chrom[t], "Y") != NULL && strlen(refC[t]) == 1 && strlen(alt[t]) == 1))  //VT is a SNP or chrY(only SNPs)
                        {
                            if (strlen(refC[t]) != 1)
                            {
                                printf("*** UNSUPPORTED CASE, variant will be omitted. ***\n");
                                printf("\tPosition: %lld. It is a SNP and REF length is not 1.\n", pos+offset);
                                variantOK = false;
                                break;
                            }
                            if (strlen(alt[t]) != 1)
                            {
                                printf("*** UNSUPPORTED CASE, variant will be omitted. ***\n");
                                printf("\tPosition: %lld. It is a SNP and ALT length is not 1\n", pos+offset);
                                
                                
                                variantOK = false;
                                break;
                            }
                            
                            
                            while(1)
                            {
                                
                                
                                returnValue = snprintf(VLbuffer[t] + VLbufPos[t], buf_max[t] - VLbufPos[t], "%lld SNP %s\n", pos+offset, alt[t]);
                                if(returnValue > 0 && returnValue < buf_max[t] - VLbufPos[t])
                                {
                                    VLbufPos[t] +=  returnValue;
                                    
                                    break;
                                }
                                else
                                {
                                    tmpChrPtr = (char *) realloc(VLbuffer[t], (buf_max[t]*2)* sizeof(char));
                                    
                                    if (NULL ==tmpChrPtr) {
                                        printf("Error in realloc, the message is - %s\n", strerror(errno));
                                        exit(8);
                                    }
                                    else
                                    {
                                        VLbuffer[t] = tmpChrPtr;
                                        buf_max[t] = buf_max[t]*2;
                                    }
                                    
                                    
                                }
                            }
                        }
                        else if (strstr(info[t], "VT=INDEL") != NULL)  //VT is an INDEL
                        {
                            if ((255&alt[t][0]) != (255&refC[t][0]))
                            {
                                printf("*** UNSUPPORTED CASE, variant will be omitted. ***\n");
                                printf("\tPosition: %lld. It is an INDEL and first characters of ALT and REF do not match.\n", pos+offset);
                                variantOK = false;
                                break;
                            }
                            
                            if(strlen(refC[t])==1) //insertion
                            {
                                while(1)
                                {
                                    returnValue = snprintf(VLbuffer[t] + VLbufPos[t], buf_max[t] - VLbufPos[t], "%lld INS %s\n", pos+offset, alt[t]+1);
                                    if(returnValue > 0 && returnValue < buf_max[t] - VLbufPos[t])
                                    {
                                        VLbufPos[t] +=  returnValue;
                                        
                                        break;
                                    }
                                    else
                                    {
                                        tmpChrPtr = (char *) realloc(VLbuffer[t], (buf_max[t]*2)* sizeof(char));
                                        if (NULL ==tmpChrPtr) {
                                            printf("Error in realloc, the message is - %s\n", strerror(errno));
                                            exit(8);
                                        }
                                        else
                                        {
                                            VLbuffer[t] = tmpChrPtr;
                                            buf_max[t] = buf_max[t] * 2;
                                        }
                                        
                                        
                                    }
                                }
                                
                                
                            }
                            else //if (strlen(alt[t])==1) //deletion
                            {
                                while(1)
                                {
                                    returnValue = snprintf(VLbuffer[t] + VLbufPos[t], buf_max[t] - VLbufPos[t],  "%lld DEL %ld\n", pos+offset, strlen(refC[t])-1);
                                    if(returnValue > 0 && returnValue < buf_max[t] - VLbufPos[t])
                                    {
                                        VLbufPos[t] +=  returnValue;
                                        
                                        break;
                                    }
                                    else
                                    {
                                        tmpChrPtr = (char *) realloc(VLbuffer[t], (buf_max[t]*2)* sizeof(char));
                                        if (NULL ==tmpChrPtr) {
                                            printf("Error in realloc, the message is - %s\n", strerror(errno));
                                            exit(8);
                                        }
                                        else
                                        {
                                            VLbuffer[t] = tmpChrPtr;
                                            buf_max[t] = buf_max[t] * 2;
                                        }
                                        
                                    }
                                }
                                
                            }
                            
                        }
                        else if (strstr(info[t], "VT=SV") != NULL)  //VT is a SV
                        {
                            if (strstr(info[t], "SVTYPE=DEL") != NULL)
                            {
                                endPos[t][0] = '\0';
                                if (strstr(info[t], ";END=") != NULL)
                                {
                                    chr_ptr = strstr(info[t], ";END=");
                                    temp_int = 5;
                                    while (chr_ptr[temp_int] != ';' &&  chr_ptr[temp_int] != '\0')
                                    {
                                        endPos[t][temp_int-5] =  chr_ptr[temp_int];
                                        temp_int++;
                                    }
                                    endPos[t][temp_int-5]='\0';
                                }
                                else if (strncmp(info[t], "END=", 4) == 0)
                                {
                                    temp_int=4;
                                    while (info[t][temp_int] != ';' &&  info[t][temp_int] != '\0')
                                    {
                                        endPos[t][temp_int-4] =  info[t][temp_int];
                                        temp_int++;
                                    }
                                    endPos[j-4][0]='\0';
                                }
                                else
                                {
                                    printf("*** UNSUPPORTED CASE, variant will be omitted. ***\n");
                                    printf("\tPosition: %lld. SVTYPE=DEL and no END in info: %s\n", pos+offset, info[t]);
                                    
                                    variantOK = false;
                                    break;
                                }
                                delLen = (uint64_t)(atol(endPos[t]) - pos);
                                
                                if (strncmp(alt[t], "<DEL>", 5) == 0) {
                                    isInsertion=false;
                                }
                                else if (strlen(alt[t]) == 1 && alt[t][0]==refC[t][0])
                                {
                                    isInsertion=false;
                                }
                                else if (strlen(alt[t]) > 1 && alt[t][0]==refC[t][0])
                                {
                                    isInsertion=true;
                                }
                                else {
                                    printf("*** UNSUPPORTED CASE, variant will be omitted. ***\n");
                                    printf("\tPosition: %lld. SVTYPE=DEL and for given REF, wrong ALT data: %s\n", pos+offset, alt[t]);
                                    
                                    variantOK = false;
                                    break;
                                }
                                
                                if (isInsertion)
                                {
                                    
                                    while(1)
                                    {
                                        returnValue = snprintf(VLbuffer[t] + VLbufPos[t], buf_max[t] - VLbufPos[t],  "%lld SV %llu %s\n", pos+offset, delLen, alt[t]+1);
                                        if(returnValue > 0 && returnValue < buf_max[t] - VLbufPos[t])
                                        {
                                            VLbufPos[t] +=  returnValue;
                                            break;
                                        }
                                        else
                                        {
                                            tmpChrPtr = (char *) realloc(VLbuffer[t], (buf_max[t]*2)* sizeof(char));
                                            if (NULL ==tmpChrPtr) {
                                                printf("Error in realloc, the message is - %s\n", strerror(errno));
                                                exit(8);
                                            }
                                            else
                                            {
                                                VLbuffer[t] = tmpChrPtr;
                                                buf_max[t] = buf_max[t] * 2;
                                            }
                                        }
                                    }
                                    
                                }
                                else
                                {
                                    
                                    while(1)
                                    {
                                        returnValue = snprintf(VLbuffer[t] + VLbufPos[t], buf_max[t] - VLbufPos[t],  "%lld SV %llu -\n", pos+offset, delLen);
                                        if(returnValue > 0 && returnValue < buf_max[t] - VLbufPos[t])
                                        {
                                            VLbufPos[t] +=  returnValue;
                                            break;
                                        }
                                        else
                                        {
                                            tmpChrPtr = (char *) realloc(VLbuffer[t], (buf_max[t]*2)* sizeof(char));
                                            if (NULL ==tmpChrPtr) {
                                                printf("Error in realloc, the message is - %s\n", strerror(errno));
                                                exit(8);
                                            }
                                            else
                                            {
                                                VLbuffer[t] = tmpChrPtr;
                                                buf_max[t] = buf_max[t] * 2;
                                            }
                                        }
                                    }
                                    
                                    
                                }
                                
                                
                                
                            }
                            else
                            {
                                printf("*** UNSUPPORTED CASE, variant will be omitted. ***\n");
                                printf("\tUnsupported variant VT=SV at position: %lld  and not SVTYPE=DEL, INFO field:  %s\n",  pos+offset, info[t]);
                                
                                variantOK = false;
                                break;
                            }
                            
                        }
                        else
                        {
                            printf("*** UNSUPPORTED CASE, variant will be omitted. ***\n");
                            printf("\tUnsupported variant at position: %lld  (neither SNP, INDEL nor SV), INFO field: %s \n", pos+offset, info[t]);
                            
                            
                            variantOK = false;
                            break;
                        }
                        
                        
                    }
                    
                    if(!variantOK)
                        continue;
                    
                    
                    
                    // if there is a change in number of bitsCalculated that requires another unsigned char for storage, reallocate
                    if ((bitsCalculated[t]+noAlt) > char_max[t]*8)
                    {
                        for (int g = 0; g < no_gen; g++)
                        {
                            // cout << g << " " << t << " " << char_max[t] << endl;
                            tmpUChrPtr = (unsigned char*) realloc(tempBits[t][g], (char_max[t] + 1)* sizeof(unsigned char));
                            
                            if (NULL ==tmpUChrPtr) {
                                printf("Error in realloc, the message is - %s\n", strerror(errno));
                                exit(8);
                            }
                            else
                            {
                                tempBits[t][g] = tmpUChrPtr;
                                tempBits[t][g][char_max[t]] = 0x0;
                            }
                        }
                        char_max[t] = char_max[t] + 1;
                        
                    }
                    
                    for (int g = 0; g < no_gen; g++)
                    {
                        
                        sscanf(line[currentLine]+ptr, "%s%n", genotype[t], &mv);
                        ptr = ptr + mv;
                        
                        
                        pp = strtok_r(genotype[t], ":", &strtokState2);
                        for(int m=0; m < gt_pos; m++)
                            pp = strtok_r( NULL, ":", &strtokState2);
                        
                        
                        
                        gt = atoi(pp);
                        
                        if (gt > noAlt)
                        {
                            printf("*** WARNING***\n");
                            printf("\tPosition: %lld, GT data for sample no %d is wrong: %s, no variant will be placed in the consensus sequence here\n",  pos+offset, currentLine, info[t]);
                            cout << genotype[t] << " " << pp <<  " " << gt << " " << noAlt << " "<< g <<endl;
                            fflush(stdout);
                            gt = 0;
                        }
                        
                        
                        
                        for (a = 0; a < noAlt; a++)
                        {
                            if(gt == a+1)
                            {
                                shift = (8*sizeof(unsigned char)-1-(bitsCalculated[t]+a)%(8*sizeof(unsigned char)));
                                
                                tempBits[t][g][(bitsCalculated[t]+a)/(8*sizeof(unsigned char))] = tempBits[t][g][(bitsCalculated[t]+a)/(8*sizeof(unsigned char))]|((unsigned char)1<<shift);
                            }
                            
                        }
                        
                    }
                    bitsCalculated[t] = bitsCalculated[t] + noAlt;
                }
            }
        }
        
        
        
        
        for (int t = 0; t < noThreads; t++)
        {
            for (int g = 0; g < no_gen; g++)
            {
                for (int b = 0; b <= ((int32_t)bitsCalculated[t]-1.0)/8.0; b++)
                {
                    
                    writeno = 8;
                    if(bitsCalculated[t]<8)
                        writeno = bitsCalculated[t];
                    if((b == (bitsCalculated[t]-1)/8 ) && bitsCalculated[t]%8 != 0)
                        writeno = bitsCalculated[t]%8;
                    
                    bv->writeBits(g, tempBits[t][g][b], writeno);
                    
                }
            }
            
            l = 0;
            
            for (int v = 0; v < bitsCalculated[t]; v++)
            {
                /* fprintf(vl_tmp_file, "%d ", vt_count + v);
                 while(VLbuffer[t][l] != '\n')
                 {
                 putc(VLbuffer[t][l++], vl_tmp_file);
                 
                 }
                 putc(VLbuffer[t][l++], vl_tmp_file);
                 
                 */
                {
                    
                    sscanf(VLbuffer[t]+l, "%u %s %n",  &vl->vt[vt_count].prePos, type, &n);
                    l += n;
                    
                    vl->vt[vt_count].prePos = vl->vt[vt_count].prePos - 1; //go from 1-base to 0-base
                    
                    vl->prevSV[vt_count] = vl->no_sv - 1;
                    if (strstr(type, "SNP") != NULL)  //VT is a SNP
                    {
                        vl->vtType[vt_count] = SNP;
                        vl->vtAux[vt_count].alt = new char[1];
                        sscanf(VLbuffer[t]+l, "%s\n%n", vl->vtAux[vt_count].alt, &n);
                        l += n;
                        vl->vtAux[vt_count].delLen = 0;
                        vl->vt[vt_count].prePos = vl->vt[vt_count].prePos - 1;  //so it points to previous position (one before SNP) as for other variants
                    }
                    else if (strstr(type, "INS") != NULL)  //VT is an INS
                    {
                        vl->vtType[vt_count] = INS;
                        sscanf(VLbuffer[t]+l, "%s\n%n", temp, &n);
                        l += n;
                        vl->vtAux[vt_count].alt = new char[strlen(temp)+1];
                        vl->no_chars_inserted += strlen(temp);
                        memcpy(vl->vtAux[vt_count].alt , temp, strlen(temp));
                        vl->vtAux[vt_count].alt[strlen(temp)]='\0';
                        //cout << "ins: " << i << " " << strlen(temp) <<  " " <<  temp << " " << vtAux[vt_count].alt << endl;
                        vl->vtAux[vt_count].delLen = 0;
                        vl->no_ins++;
                        
                    }
                    else if (strstr(type, "DEL") != NULL)  //VT is an DEL
                    {
                        vl->vtType[vt_count] = DEL;
                        sscanf(VLbuffer[t]+l, "%u\n%n", &vl->vtAux[vt_count].delLen, &n);
                        l += n;
                        vl->vtAux[vt_count].alt = new char[1];
                        vl->vtAux[vt_count].alt[0] = '-';
                        
                        if(vl->vtAux[vt_count].delLen > vl->max_del)
                            vl->max_del = vl->vtAux[vt_count].delLen;
                        
                        vl->no_del++;
                    }
                    else if (strstr(type, "SV") != NULL)  //VT is an SV
                    {
                        vl->vtType[vt_count] = SV;
                        sscanf(VLbuffer[t]+l, "%u %s\n%n", &vl->vtAux[vt_count].delLen, temp, &n);
                        l += n;
                        vl->vtAux[vt_count].alt = new char[strlen(temp)+1];
                        memcpy(vl->vtAux[vt_count].alt , temp, strlen(temp));
                        vl->vtAux[vt_count].alt[strlen(temp)]='\0';
                        
                        
                        if(vl->vtAux[vt_count].alt[0] != '-')
                            vl->no_chars_inserted += strlen(temp);
                        
                        // cout << "ins: " << i << " " << strlen(temp) <<  " " <<  temp << " " << vtAux[vt_count].alt << endl;
                        vl->no_sv++;
                        
                        if(vl->vtAux[vt_count].delLen > vl->max_sv_del)
                            vl->max_sv_del = vl->vtAux[vt_count].delLen;
                        
                    }
                }
                
                vt_count++;
            }
            
            //vt_count += bitsCalculated[t];
            
        }
        
        
        if (chl != EOF)
            chl = getc(file);
        fseek( file, -1, SEEK_CUR);
        //  cout << "vt_count " << vt_count << endl;
    }
    
    
    for (i = 0; i < NO_LINES_PER_THREAD * noThreads; i++)
        free(line[i]);
    free(line);
    for (i = 0; i < noThreads; i++)
    {
        for (j = 0; j < no_gen; j++)
        {
            free(tempBits[i][j]);
        }
        
        free(tempBits[i]);
        free(VLbuffer[i]);
        free(chrom[i]);
        free(name[i]);
        free(refC[i]);
        free(altList[i]);
        free(alt[i]);
        free(filter[i]);
        free(info[i]);
        free(format[i]);
        free(qual[i]);
        free(genotype[i]);
        free(endPos[i]);
        
        
    }
    
    free(tempBits);
    free(VLbuffer);
    free(VLbufPos);
    free(chrom);
    free(name);
    free(refC);
    free(altList);
    free(alt);
    free(filter);
    free(info);
    free(format);
    free(qual);
    free(genotype);
    free(endPos);
    free(bitsCalculated);
    free(buf_max);
    
    // fclose(vl_tmp_file);
    
    return true;
    
    
}

