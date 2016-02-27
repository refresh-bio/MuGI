/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */

#ifndef _STAT_H
#define _STAT_H

#include "timer.h"

//**************************************************************
// Functions producing statistics of usage
#ifdef __WIN32__
#include <windows.h>
/*typedef struct {
 LARGE_INTEGER start;
 LARGE_INTEGER stop;
 } stopWatch;
 */

class TimeMeasure {
    
private:
    stopWatch timer;
    LARGE_INTEGER frequency;
public:
    TimeMeasure()
    {
        timer.start.QuadPart=0;
        timer.stop.QuadPart=0;
        QueryPerformanceFrequency( &frequency );
    }
    
    double LIToSecs( LARGE_INTEGER & L) {
        return ((double)L.QuadPart /(double)frequency.QuadPart);
    }
    
    void start_clock( )
    {
        QueryPerformanceCounter(&timer.start);
    }
    void stop_clock( )
    {
        QueryPerformanceCounter(&timer.stop);
    };
    //**************************************************************
    // Returns the measured time
    double get_time(void)
    {
        LARGE_INTEGER time;
        time.QuadPart = timer.stop.QuadPart - timer.start.QuadPart;
        return 1000.0 * LIToSecs( time) ;
    }
    
    //**************************************************************
    // Resets the measured time
    void reset_clock(void)
    {
    }
};
#endif

//#pragma once
//#include "config.h"


#include <iostream>
#include <stdlib.h>
#include <algorithm>

#include <vector>
#include <numeric>


typedef unsigned int  t_symbol;

typedef struct
{
    int32_t lcts;
    double time;
} t_stat;


template<typename T> class SDStat
{
	std::vector<T> stat;
public:
	SDStat()
	{
		stat.clear();
	};
	
	void add(T x)
	{
		stat.push_back(x);
	};
	
	void reset(void)
	{
		stat.clear();
	};
    
    void sort(void)
	{
		std::sort(stat.begin(), stat.end());
	};
	
	T get_median(void)
	{
		if(stat.empty())
			return 0;
        
		
        //	std::sort(stat.begin(), stat.end());
		return stat[stat.size() / 2];
	};
    
    
    T get_time_for_percent(float value)
	{
		if(stat.empty())
			return 0;
		if(value < 100)
            return stat[stat.size() * value / 100];
        else
            return stat[stat.size() - 1];
	};
    
    
       
	double get_average(void)
	{
		if(stat.empty())
			return 0.0;
		
		T r = std::accumulate(stat.begin(), stat.end(), (T) 0);
		return (double) r / stat.size();
	};
	double get_std_dev(double mean)
	{
		if(stat.empty())
			return 0.0;
		
		double r = 0.0;
		typename std::vector<T>::iterator p;
		for(p= stat.begin(); p != stat.end(); ++p)
			r += (*p - mean) * (*p - mean);
		return sqrt(r / stat.size());
	}
};



#ifdef __APPLE__
//**************************************************************
// Klasa s≥uøπca do pomiaru czasu
class TimeMeasure
{
	double val;
	clock_t t1, t2;
	
public:
	//**************************************************************
	TimeMeasure()
	{};
	
	//**************************************************************
	// Starts time measure
	void start_clock(void)
	{
		t1 = clock();
		val = 0.0;
	};
	
	//**************************************************************
	// Stops time measure
	double stop_clock(void)
	{
		t2 = clock();
		//std::cout << "t1=" << t1 << " t2=" << t2 << "\n";
		//std::cout << " t2-t1 = " << (double) ((double)t2 - (double)t1) << "\n";
		val += ((double)1000.0 * (double) ((double)t2 - (double)t1))/(double)CLOCKS_PER_SEC;
		return val;
	};
	
	//**************************************************************
	// Continues time measure
	void continue_clock(void)
	{
		t1 = clock();
	};
	
	//**************************************************************
	// Returns the measured time
	double get_time(void)
	{
		return val;
	}
	
	//**************************************************************
	// Resets the measured time
	void reset_clock(void)
	{
		val = 0.0;
	}
};

#endif



#ifdef __linux__
//**************************************************************
// Klasa s≥uøπca do pomiaru czasu
class TimeMeasure
{
	double val;
	struct timespec t1, t2;
	
public:
	//**************************************************************
	TimeMeasure()
	{};
	
	//**************************************************************
	// Starts time measure
	void start_clock(void)
	{
		clock_gettime( CLOCK_REALTIME, &t1);
		val = 0.0;
	};
	
	//**************************************************************
	// Stops time measure
	double stop_clock(void)
	{
		clock_gettime( CLOCK_REALTIME, &t2);
		//std::cout << "t1=" << t1 << " t2=" << t2 << "\n";
		//std::cout << " t2-t1 = " << (double) ((double)t2 - (double)t1) << "\n";
		val += ((t2.tv_sec - t1.tv_sec) + (double)( t2.tv_nsec - t1.tv_nsec )/ (double)1000000000L)*1000;
		return val;
	};
	
	//**************************************************************
	// Continues time measure
	void continue_clock(void)
	{
		clock_gettime( CLOCK_REALTIME, &t1);
	};
	
	//**************************************************************
	// Returns the measured time
	double get_time(void)
	{
		return val;
	}
	
	//**************************************************************
	// Resets the measured time
	void reset_clock(void)
	{
		val = 0.0;
	}
};

#endif

#endif
