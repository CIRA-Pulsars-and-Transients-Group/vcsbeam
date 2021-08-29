/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "performance.h"


/************************************************
 * Helper functions not listed in performance.h *
 ************************************************/

double now()
/* Get the precise current time in seconds as a double
 */
{
    struct timespec t;
    clock_gettime( CLOCK_REALTIME, &t );
    return (double)t.tv_sec + (double)t.tv_nsec/1000000000L;
}

int get_stopwatch_idx( logger *log, const char *stopwatch_name )
/* Find the stopwatch with the given name
 */
{
    int idx;
    for (idx = 0; idx < log->nstopwatches; idx++)
    {
        if (strcmp( stopwatch_name, log->stopwatches[idx].name ) == 0)
            return idx;
    }

    return PERFORMANCE_NO_STOPWATCH_FOUND;
}

inline double calc_stopwatch_mean( logger_stopwatch *sw )
{
    return sw->total / (double)sw->nstart_stops;
}

inline double calc_stopwatch_std( logger_stopwatch *sw )
{
    double mean = calc_stopwatch_mean( sw );
    return sqrt(sw->total_sq / (double)sw->nstart_stops - mean*mean);
}

void write_stopwatch_stats_str( logger *log, logger_stopwatch *sw )
{
    char report[56 + strlen(sw->name)];
    sprintf( report, "Stopwatch \"%s\": "
            "Total = %.3lf s,  Mean = %.3lf +/- %.3lf",
            sw->name,
            sw->total,
            calc_stopwatch_mean( sw ),
            calc_stopwatch_std( sw )
           );
    logger_message( log, report );
}

/******************************************************
 * Functions for memory management and initialisation *
 ******************************************************/

logger *create_logger( FILE *fout )
{
    // Allocate memory
    logger *log = (logger *)malloc( sizeof(logger) );

    // Set values
    log->begintime    = now();
    log->fout         = fout;
    log->nstopwatches = 0;

    int i;
    for (i = 0; i < PERFORMANCE_MAX_NUM_STOPWATCHES; i++)
    {
        log->stopwatches[i].name         = NULL;
        log->stopwatches[i].nstart_stops = 0;
        log->stopwatches[i].running      = false;
        log->stopwatches[i].total        = 0.0;
        log->stopwatches[i].total_sq     = 0.0;
    }

    // Return pointer to new struct
    return log;
}


void destroy_logger( logger *log )
{
    int i;
    for (i = 0; i < log->nstopwatches; i++)
    {
        free( log->stopwatches[i].name );
    }

    free( log );
}


/*******************************************************
 * Functions for manipulating user-defined stopwatches *
 *******************************************************/

void logger_add_stopwatch( logger *log, const char *stopwatch_name )
{
    // Make sure the maximum allowed number of stopwatches has not been reached
    if (log->nstopwatches == PERFORMANCE_MAX_NUM_STOPWATCHES)
    {
        fprintf( stderr, "warning: logger_add_stopwatch: "
                "maximum number of stopwatches (%d) reached, \"%s\" not added\n",
                PERFORMANCE_MAX_NUM_STOPWATCHES, stopwatch_name );
        return;
    }

    // Make sure no stopwatch with the given name already exists
    if (get_stopwatch_idx( log, stopwatch_name ) != PERFORMANCE_NO_STOPWATCH_FOUND)
    {
        fprintf( stderr, "warning: logger_add_stopwatch: "
                "\"%s\" already exists, not added\n", stopwatch_name );
        return;
    }

    // Add the stopwatch!
    log->stopwatches[log->nstopwatches++].name = (char *)malloc( strlen(stopwatch_name) + 1 );
    strcpy( log->stopwatches[log->nstopwatches++].name, stopwatch_name );
}


void logger_start_stopwatch( logger *log, const char *stopwatch_name )
{
    int idx = get_stopwatch_idx( log, stopwatch_name );

    // Make sure a stopwatch with the given name exists
    if (idx == PERFORMANCE_NO_STOPWATCH_FOUND)
    {
        fprintf( stderr, "warning: logger_start_stopwatch: "
                "\"%s\" not found\n", stopwatch_name );
        return;
    }

    logger_stopwatch *sw = &(log->stopwatches[idx]);

    // Check to see if the stopwatch is already running
    if (log->stopwatches[idx].running)
    {
        fprintf( stderr, "warning: logger_start_stopwatch: "
                "\"%s\" is already running\n", sw->name );
        return;
    }

    // Check that we haven't exhausted the number of allowed runs
    if (log->stopwatches[idx].nstart_stops == PERFORMANCE_MAX_START_STOP)
    {
        fprintf( stderr, "warning: logger_start_stopwatch: "
                "\"%s\" has reached max allowed number of runs\n", stopwatch_name );
        return;
    }

    // Put the current time in the next "value"
    int n = sw->nstart_stops;

    sw->values[n] = now();
    sw->running   = true;
}


void logger_stop_stopwatch( logger *log, const char *stopwatch_name )
{
    int idx = get_stopwatch_idx( log, stopwatch_name );

    // Make sure a stopwatch with the given name exists
    if (idx == PERFORMANCE_NO_STOPWATCH_FOUND)
    {
        fprintf( stderr, "warning: logger_stop_stopwatch: "
                "\"%s\" not found\n", stopwatch_name );
        return;
    }

    logger_stopwatch *sw = &(log->stopwatches[idx]);

    // Check to see if the stopwatch is running
    if (sw->running == false)
    {
        fprintf( stderr, "warning: logger_stop_stopwatch: "
                "\"%s\" is not running\n", sw->name );
        return;
    }

    // Record the time and increment the count
    int    n   = sw->nstart_stops;
    double val = now() - sw->values[n];

    sw->values[n]  = val;
    sw->running    = false;
    sw->total     += val;
    sw->total_sq  += val*val;
    sw->nstart_stops++;
}


/********************************************************************
 * Functions for writing out different kinds of messages to the log *
 ********************************************************************/

void logger_timed_message( logger *log, const char *message )
{
    double current_time = now() - log->begintime;
    fprintf( log->fout, "[%f]  %s\n", current_time, message );
}

void logger_message( logger *log, const char *message )
{
    fprintf( log->fout, "%s\n", message );
}

void logger_stopwatch_report_stats( logger *log, const char *stopwatch_name )
{
    int idx = get_stopwatch_idx( log, stopwatch_name );

    // Make sure a stopwatch with the given name exists
    if (idx == PERFORMANCE_NO_STOPWATCH_FOUND)
    {
        fprintf( stderr, "warning: logger_stop_stopwatch: "
                "\"%s\" not found\n", stopwatch_name );
        return;
    }

    logger_stopwatch *sw = &(log->stopwatches[idx]);

    write_stopwatch_stats_str( log, sw );
}

void logger_report_all_stats( logger *log )
{
    int i;
    for (i = 0; i < log->nstopwatches; i++)
    {
        write_stopwatch_stats_str( log, &(log->stopwatches[i]) );
    }
}
