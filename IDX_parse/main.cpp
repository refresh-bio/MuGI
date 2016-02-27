/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */


#define NAME_MAX 1000

#include <iostream>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

using namespace std;


typedef enum   {hexa, full, old} out_mode;
int main(int argc, const char * argv[])
{

    out_mode mode = hexa;
    if(!argv[1] || !argv[2] || !argv[3])
    {
        printf("Wrong parameters! Usage:\n./parse [dict] [file_to_parse] [output_file] [mode]\nPossible modes:\n0 - hex\n1 - full\n");
        exit(1);
    }
    FILE * RES = fopen(argv[2], "rb");
    if (RES == NULL) {
        printf("Cannot open %s\n", argv[2]);
        printf("The message is - %s\n", strerror(errno));
        exit(8);
    }
    FILE * OUT = fopen(argv[3], "w");
    if (OUT == NULL) {
        printf("Cannot open %s\n", argv[3]);
        printf("The message is - %s\n", strerror(errno));
        exit(8);
    }

    ///**************** Read dictionary file  *******/
    unsigned char eol_type;
    FILE * DICT = fopen(argv[1], "rt");
    if (DICT == NULL) {
        printf("Cannot open %s\n", argv[1]);
        printf("The message is - %s\n", strerror(errno));
        exit(8);
    }
    
    int no_chr;
    int chr;
    fscanf(DICT, "%d\n", &no_chr);
    
    string * chr_name_long = new string[no_chr];
    string * chr_name = new string[no_chr];
    uint32_t * chr_start = new uint32_t[no_chr];
   
    
    uint32_t * chr_len = new uint32_t[no_chr];
    unsigned char c;
    for(chr = 0; chr < no_chr; chr++)
    {
        while(true)
        {
            c = getc(DICT);
            chr_name_long[chr].push_back(c);
            if(c == '\n' || c == '\r')
                break;
        }
        
        
        eol_type = c;
        
        
        for(int i = 0; ; ++i)
        {
            c = chr_name_long[chr].c_str()[i];
            if(c == ' ' || c == eol_type )
                break;
            
            chr_name[chr].push_back(c);
            
            
        }
        
        fscanf(DICT, "%u\t%u\n", &chr_start[chr], &chr_len[chr]);
        
    }
    int no_ind;
    fscanf(DICT, "%d\n", &no_ind);
    string * ind_name = new string[no_ind];
    for(int id = 0; id < no_ind; id++)
    {
        while(true)
        {
            c = getc(DICT);
            if(c == '\n' || c == '\r')
                break;
            ind_name[id].push_back(c);
        }
    }
    
    //END Read dictionary file ////////
    
    
    if(argv[4])
    {
        if(atoi(argv[4]) == 0)
            mode = hexa;
        else if(atoi(argv[4]) == 1)
            mode = full;
        /*else if(atoi(argv[4]) == 2)
            mode = old;*/
        else
        {
            cout << "mode " << argv[4] << " not supported" << endl;
            cout << "Possible modes:\n0 - hex (default)\n1 - full\n";
        }
    }
    
    
    
    uint32_t noID = no_ind;
    uint32_t nameMax = NAME_MAX;
    char * name = (char*) malloc (sizeof(char) * nameMax);
    char * tempName;
    uint64_t read64, help;
    uint32_t read32;
    size_t ret;
    
    uint32_t nameLen;
    
    ret = fread(&read32, sizeof(uint32_t), 1, RES);
    uint64_t i = 0;
    if(mode == hexa)
    {
        while(read32 == 0xFFFFFFFF && ret != 0)
        {
            
            
            ret = fread(&nameLen, sizeof(uint32_t), 1, RES);
            
            if(nameLen == 0xFFFFFFFF)
            {
                fprintf(OUT, "%llu\n", i);
            }
            else
            {
                if(nameLen > nameMax)
                {
                    nameMax = nameLen + 1;
                    tempName = (char*) realloc(name, sizeof(char) * nameMax);
                    if(!tempName)
                    {
                        printf("Memory error\n");
                        return 1;
                    }
                    else
                        name = tempName;
                }
                
                fread(name, sizeof(char), nameLen, RES);
                fprintf(OUT, "%s\n", name);
            }
            
            
            
            
            
            //pos
            ret = fread(&read32, sizeof(uint32_t), 1, RES);
            
            
            
            if(read32 == 0xFFFFFFFF) //not found
            {
                ret = fread(&read32, sizeof(uint32_t), 1, RES);
                i++;
                continue;
            }
            
            while(read32 != 0xFFFFFFFF && ret != 0)
            {
                //fprintf(OUT, "%u\n", read32);
                chr = 0;
                if(read32 != 0xFFFFFFFE)
                {
                    while(read32 >= (chr_start[chr] + chr_len[chr]))
                        chr++;
                    fprintf(OUT, "chr %s, pre_pos %u\n", chr_name[chr].c_str(), read32 - chr_start[chr]);
                }
                else
                    fprintf(OUT, "chr %s, pre_pos %i\n", chr_name[chr].c_str(), -1);
            
                
               
                

                
                
                int32_t r = 0;
                
                ret = fread(&read32, sizeof(uint32_t), 1, RES);
                if(read32 == 0xFFFFFFFF)
                    fprintf(OUT, "ALL\n");/*for(int32_t r=0; r < 35; r++)
                    {
                        fprintf(OUT, "0xFFFFFFFFFFFFFFFF");
                    }*/
                else
                {
                    for(r=0; r < noID/64; r++)
                    {
                        ret = fread(&read64, sizeof(uint64_t), 1, RES);
                        fprintf(OUT, "%016llX", read64);
                    }
                    // last only 8 individuals (2184-(64*34))
                    ret = fread(&read64, sizeof(uint64_t), 1, RES);
                    for(r = 0; r < (noID-(64*(noID/64)))/4; r++)
                    {
                        help = (read64 >> ((15-r)*4)) & 0xF;
                       // help = help & read64;
                        fprintf(OUT, "%llX", help);
                    }
                    switch ((noID-(64*(noID/64)))%4)
                    {
                        case 0:
                            break;
                        case 1:
                            help = (read64 >> ((15-r)*4)) & 0x8;
                            fprintf(OUT, "%llX", help);
                            break;
                        case 2:
                            help = (read64 >> ((15-r)*4)) & 0xC;
                            fprintf(OUT, "%llX", help);
                            break;
                        case 3:
                            help = (read64 >> ((15-r)*4)) & 0xE;
                            fprintf(OUT, "%llX", help);
                            break;
                    }
                    
                    fprintf(OUT, "\n");
                }
                ret = fread(&read32, sizeof(uint32_t), 1, RES);
            }
            
            
            i++;
            
        }
    }
    else if (mode == full)
    {
        while(read32 == 0xFFFFFFFF && ret != 0)
        {
            
            
            ret = fread(&nameLen, sizeof(uint32_t), 1, RES);
            
            if(nameLen == 0xFFFFFFFF)
            {
                fprintf(OUT, "#%llu\n", i);
            }
            else
            {
                if(nameLen > nameMax)
                {
                    nameMax = nameLen + 1;
                    tempName = (char*) realloc(name, sizeof(char) * nameMax);
                    if(!tempName)
                    {
                        printf("Memory error\n");
                        return 1;
                    }
                    else
                        name = tempName;
                }
                
                fread(name, sizeof(char), nameLen, RES);
                fprintf(OUT, "%s\n", name);
            }
            
            
            
            
            
            //pos
            ret = fread(&read32, sizeof(uint32_t), 1, RES);
            
            if(read32 == 0xFFFFFFFF) //not found
            {
                ret = fread(&read32, sizeof(uint32_t), 1, RES);
                i++;
                continue;
            }
            
            while(read32 != 0xFFFFFFFF && ret != 0)
            {
                chr = 0;
                if(read32 != 0xFFFFFFFE)
                {       while(read32 > chr_start[chr] + chr_len[chr])
                        chr++;
                    fprintf(OUT, "chr %s, pre_pos %u\n", chr_name[chr].c_str(), read32-chr_start[chr]);
                }
                else
                     fprintf(OUT, "chr %s, pre_pos %i\n", chr_name[chr].c_str(), -1);
            
                int32_t r =0;
                
                ret = fread(&read32, sizeof(uint32_t), 1, RES);
                if(read32 == 0xFFFFFFFF)
                    fprintf(OUT, "ALL\n");
                else
                {
                    for(r=0; r < noID/64; r++)
                    {
                        ret = fread(&read64, sizeof(uint64_t), 1, RES);
                    
                        for(int s = 0; s < 64; s++)
                            if((1UL << (63-s)) & read64)
                                fprintf(OUT, "%s ", ind_name[r*64+s].c_str());
                        
                        // fprintf(OUT, "%016llX", read64);
                    }
                    // last [only 8 individuals for 1000GP (2184-(64*34))]
                    ret = fread(&read64, sizeof(uint64_t), 1, RES);
                    for(r = 0; r < (noID-(64*(noID/64))); r++)
                    {
                        if((1UL << (63-r)) & read64)
                            fprintf(OUT, "%s ", ind_name[noID/64*64+r].c_str());
                    }
                    
                    fflush(stdout);
                    fprintf(OUT, "\n");
                }
                ret = fread(&read32, sizeof(uint32_t), 1, RES);
            }
            
            
            i++;
            
        }
    }
    
  
    return 0;
}

