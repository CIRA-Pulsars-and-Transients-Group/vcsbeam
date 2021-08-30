/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#ifndef __PERFORMANCE_H__
#define __PERFORMANCE_H__

#define PERFORMANCE_MAX_NUM_STOPWATCHES  16
#define PERFORMANCE_MAX_START_STOP       4096

#define PERFORMANCE_NO_STOPWATCH_FOUND  -1

#include <stdio.h>
#include <time.h>

typedef struct logger_stopwatch_t
{
    char    *name;
    double   values[PERFORMANCE_MAX_START_STOP];
    int      nstart_stops;
    bool     running;
    double   total, total_sq; // For calculation of stats
} logger_stopwatch;

typedef struct logger_t
{
    double            begintime;
    FILE             *fout;
    logger_stopwatch  stopwatches[PERFORMANCE_MAX_NUM_STOPWATCHES];
    int               nstopwatches;
    int               world_rank;
} logger;

/******************************************************
 * Functions for memory management and initialisation *
 ******************************************************/

logger *create_logger( FILE * fout, int world_rank );
void destroy_logger( logger *log );

/**************************************************
 * Functions for manipulating user-defined stopwatchs *
 **************************************************/

void logger_add_stopwatch( logger *log, const char *stopwatch_name );
void logger_start_stopwatch( logger *log, const char *stopwatch_name );
void logger_stop_stopwatch( logger *log, const char *stopwatch_name );

/********************************************************************
 * Functions for writing out different kinds of messages to the log *
 ********************************************************************/

void logger_timed_message( logger *log, const char *message );
void logger_message( logger *log, const char *message );
void logger_stopwatch_report_stats( logger *log, const char *stopwatch_name );
void logger_report_all_stats( logger *log );

#endif
