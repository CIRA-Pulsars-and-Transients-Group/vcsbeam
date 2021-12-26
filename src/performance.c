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

#include "vcsbeam.h"


/********************************************
 * Helper functions not listed in vcsbeam.h *
 *******************************************/

/**
 * Gets the precise current time.
 *
 * @return The current time in seconds
 */
double now()
{
    struct timespec t;
    clock_gettime( CLOCK_REALTIME, &t );
    return (double)t.tv_sec + (double)t.tv_nsec/1000000000L;
}

/**
 * Finds the stopwatch with the given name.
 *
 * @param log The logger object to be searched
 * @param[in] The string to be searched
 * @return The index into `log&rarr;stopwatches` where the stopwatch name
 *         matches `stopwatch_name`
 */
int get_stopwatch_idx( logger *log, const char *stopwatch_name )
{
    int idx;
    for (idx = 0; idx < log->nstopwatches; idx++)
    {
        if (strcmp( stopwatch_name, log->stopwatches[idx].name ) == 0)
            return idx;
    }

    return PERFORMANCE_NO_STOPWATCH_FOUND;
}

/**
 * Finds the maximum length of stopwatch names.
 *
 * @param log A logger object with stopwatches
 * @return The length of the longest name in the stopwatches in `log`
 */
int get_stopwatch_max_name_length( logger *log )
{
    int length, max_length = 0;
    int idx;
    for (idx = 0; idx < log->nstopwatches; idx++)
    {
        length = strlen( log->stopwatches[idx].name );
        if (length > max_length)
            max_length = length;
    }

    return max_length;
}

/**
 * Calculates the mean duration of a stopwatch's runs.
 *
 * @param sw A stopwatch object
 * @return The mean duration of the runs of `sw`
 */
double calc_stopwatch_mean( logger_stopwatch *sw )
{
    return sw->total / (double)sw->nstart_stops;
}

/**
 * Calculates the standard devation of a stopwatch's runs.
 *
 * @param sw A stopwatch object
 * @return The standard deviation of the runs of `sw`
 */
double calc_stopwatch_std( logger_stopwatch *sw )
{
    double mean = calc_stopwatch_mean( sw );
    return sqrt(sw->total_sq / (double)sw->nstart_stops - mean*mean);
}

/**
 * Prints a stopwatch's statisticss (mean and std).
 *
 * @param log       A logger object
 * @param sw        The stopwatch object whose statistics are to be printed
 * @param pad_size  The field width for the stopwatch name (for alignment
 *                  purposes)
 *
 * Setting `pad_size = 0` is the same as setting it to the number of
 * characters in the stopwatch's name.
 */
void write_stopwatch_stats_str( logger *log, logger_stopwatch *sw, int pad_size )
{
    if (pad_size == 0)
        pad_size = strlen( sw->name );

    char report[56 + strlen(sw->name)];
    sprintf( report, "Stopwatch %-*s: "
            "Total = %.3lf s,  Mean = %.3lf +/- %.3lf",
            pad_size,
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

/**
 * Creates a new logger object.
 *
 * @param fout        The file stream where to print logger messages
 * @param world_rank  The rank of the current MPI process
 * @return A pointer to a newly created logger object
 *
 * If `world_rank` is a non-negative integer, the value of `world_rank` will
 * be included in log messages.
 * To avoid printing the rank, set `world_rank` to `PERFORMANCE_NO_MPI`.
 */
logger *create_logger( FILE *fout, int world_rank )
{
    // Allocate memory
    logger *log = (logger *)malloc( sizeof(logger) );

    // Set values
    log->begintime    = now();
    log->fout         = fout;
    log->nstopwatches = 0;
    log->world_rank   = world_rank;

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

/**
 * Destroys a logger object.
 *
 * @param log The logger object to be destroyed and its memory freed.
 */
void destroy_logger( logger *log )
{
    int i;
    for (i = 0; i < log->nstopwatches; i++)
    {
        free( log->stopwatches[i].name );
        if (log->stopwatches[i].description != NULL)
            free( log->stopwatches[i].description );
    }

    free( log );
}


/*******************************************************
 * Functions for manipulating user-defined stopwatches *
 *******************************************************/

/**
 * Adds a stopwatch to the logger object.
 *
 * @param log             The logger object to be added to
 * @param stopwatch_name  A string specifying the new stopwatch's name
 * @param description     A string specifying the new stopwatch's purpose
 */
void logger_add_stopwatch( logger *log, const char *stopwatch_name, const char *description )
{
    // Stopwatch name cannot be NULL
    if (stopwatch_name == NULL)
    {
        fprintf( stderr, "warning: logger_add_stopwatch: "
                "stopwatch_name cannot be NULL\n" );
        return;
    }

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

    // Set the stopwatch name...
    log->stopwatches[log->nstopwatches].name = (char *)malloc( strlen(stopwatch_name) + 1 );
    strcpy( log->stopwatches[log->nstopwatches].name, stopwatch_name );

    // ...and description
    if (description != NULL)
    {
        log->stopwatches[log->nstopwatches].description = (char *)malloc( strlen(description) + 1 );
        strcpy( log->stopwatches[log->nstopwatches].description, description );
    }
    else
        log->stopwatches[log->nstopwatches].description = NULL;

    // Increment the count of stopwatches
    log->nstopwatches++;
}

/**
 * Starts a stopwatch.
 *
 * @param log                The logger object containing the stopwatch to be
 *                           started
 * @param stopwatch_name     The name of the stopwatch to be started
 * @param print_description  Print a timed message with the stopwatch's
 *                           description
 *
 * If the stopwatch is already running, this has no effect.
 */
void logger_start_stopwatch( logger *log, const char *stopwatch_name, bool print_description )
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
        return;

    // Check that we haven't exhausted the number of allowed runs
    if (log->stopwatches[idx].nstart_stops == PERFORMANCE_MAX_START_STOP)
    {
        fprintf( stderr, "warning: logger_start_stopwatch: "
                "\"%s\" has reached max allowed number of runs\n", stopwatch_name );
        return;
    }

    // Print stopwatch description, if requested
    if (print_description)
        logger_timed_message( log, log->stopwatches[idx].description );

    // Put the current time in the next "value"
    int n = sw->nstart_stops;

    sw->values[n] = now();
    sw->running   = true;
}

/**
 * Stops a stopwatch.
 *
 * @param log                The logger object containing the stopwatch to be
 *                           stopped
 * @param stopwatch_name     The name of the stopwatch to be stopped
 *
 * If the stopwatch is already stopped, this has no effect.
 */
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
        return;

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

/**
 * Prints a timestamped log message.
 *
 * @param log      A logger object
 * @param message  The message to be printed
 *
 * The message will be printed to the file stream `log->fout`.
 *
 * If `log->world_rank` is non-negative, the format of the message is
 * ```
 * [<timestamp>] [MPI:<world_rank>]  <message>
 * ```
 * Otherwise, the format is
 * ```
 * [timestamp]  <message>
 * ```
 */
void logger_timed_message( logger *log, const char *message )
{
    double current_time = now() - log->begintime;
    if (log->world_rank < 0)
        fprintf( log->fout, "[%f]  %s\n", current_time, message );
    else
        fprintf( log->fout, "[%f] [MPI:%2d]  %s\n", current_time, log->world_rank, message );
}

/**
 * Prints a log message.
 *
 * @param log      A logger object
 * @param message  The message to be printed
 *
 * The message is printed to `log->fout`, including a trailing newline
 * character.
 */
void logger_message( logger *log, const char *message )
{
    fprintf( log->fout, "%s\n", message );
}

/**
 * Prints a stopwatch's statistics.
 *
 * @param log            A logger object containing the stopwatch to be
 *                       reported
 * @param stopwatch_name The name of the stopwatch to be reported
 *
 * A wrapper for write_stopwatch_stats_str(), where `pad_size` is set to zero,
 * and the statistics are only printed if there has been at least one run for
 * the specified stopwatch.
 */
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

    // Only write something if this stopwatch was actually used
    if (sw->nstart_stops > 0)
        write_stopwatch_stats_str( log, sw, 0 );
}

/**
 * Prints all stopwatches' statistics.
 *
 * @param log  the logger object whose stopwatches are to be reported
 *
 * This function prints the statistics of all stopwatches which have been run
 * at least once.
 * The output is aligned, with each stopwatch name being padded by the same
 * amount.
 *
 * Example:
 * ```
 * TO DO
 * ```
 *
 * @todo Fill in an example of stats output for logger_report_all_stats().
 */
void logger_report_all_stats( logger *log )
{
    int max_length = get_stopwatch_max_name_length( log );
    int i;
    for (i = 0; i < log->nstopwatches; i++)
    {
        if (log->stopwatches[i].nstart_stops > 0)
            write_stopwatch_stats_str( log, &(log->stopwatches[i]), max_length );
    }
}
