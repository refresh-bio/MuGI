//#include "stdafx.h"
/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

  The source codes are based on codes written by Dennis and published:
    http://allmybrain.com/2008/06/10/timing-cc-code-on-linux/
 
  Version: 0.3
  Date   : 2012-11-29
*/

#ifdef WIN32
#include <windows.h>
#endif

#include <cstdio> // NULL
#include "timer.h"


#ifdef WIN32
double CStopWatch::LIToSecs( LARGE_INTEGER & L) {
	return ((double)L.QuadPart /(double)frequency.QuadPart);
}

CStopWatch::CStopWatch(){
	timer.start.QuadPart=0;
	timer.stop.QuadPart=0;	
	QueryPerformanceFrequency( &frequency );
}

void CStopWatch::startTimer( ) {
    QueryPerformanceCounter(&timer.start);
}

void CStopWatch::stopTimer( ) {
    QueryPerformanceCounter(&timer.stop);
}


double CStopWatch::getElapsedTime() {
	LARGE_INTEGER time;
	time.QuadPart = timer.stop.QuadPart - timer.start.QuadPart;
    return LIToSecs( time) ;
}
#else

void CStopWatch::startTimer( ) {
	gettimeofday(&(timer.start),NULL);
}

void CStopWatch::stopTimer( ) {
	gettimeofday(&(timer.stop),NULL);
}

double CStopWatch::getElapsedTime() {	
	timeval res;
	timersub(&(timer.stop),&(timer.start),&res);
   return (double) ((double)res.tv_sec + (double)res.tv_usec/1000000.0)*1000;
}

#endif