/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */

#ifndef _DEFS_H
#define _DEFS_H

#define DEFAULT_MAX_ERROR 0


#ifdef WIN32
#define my_fopen	fopen
#define my_fseek	_fseeki64
#define my_ftell	_ftelli64
typedef unsigned short int uint16;
typedef int int32;
typedef short int int16;
typedef unsigned  int uint32;
typedef long long int64;
typedef unsigned  long long uint64;
typedef unsigned char uchar;
#else
#define my_fopen	fopen
#define my_fseek	fseek
#define my_ftell	ftell
#define _TCHAR	char
#define _tmain	main

typedef unsigned short int uint16;
typedef int int32;
typedef short int int16;
typedef unsigned  int uint32;
typedef long long int64;
typedef unsigned long  long uint64;
typedef unsigned char uchar;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#endif 

#endif
