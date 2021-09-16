/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#ifndef __VCSBEAM_H__
#define __VCSBEAM_H__

#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <cuComplex.h>

/*******************
 *                 *
 *   PERFORMANCE   *
 *                 *
 *******************/

#define PERFORMANCE_MAX_NUM_STOPWATCHES  16
#define PERFORMANCE_MAX_START_STOP       4096

#define PERFORMANCE_NO_STOPWATCH_FOUND  -1

#define PERFORMANCE_NO_MPI  -1

typedef struct logger_stopwatch_t
{
    char    *name;
    char    *description;
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

/*******************************************************
 * Functions for manipulating user-defined stopwatches *
 *******************************************************/

void logger_add_stopwatch( logger *log, const char *stopwatch_name, const char *description );
void logger_start_stopwatch( logger *log, const char *stopwatch_name, bool print_description );
void logger_stop_stopwatch( logger *log, const char *stopwatch_name );

/********************************************************************
 * Functions for writing out different kinds of messages to the log *
 ********************************************************************/

void logger_timed_message( logger *log, const char *message );
void logger_message( logger *log, const char *message );
void logger_stopwatch_report_stats( logger *log, const char *stopwatch_name );
void logger_report_all_stats( logger *log );

/*******************
 *                 *
 *     FILTER      *
 *                 *
 *******************/

typedef enum filter_type_t
{
    ANALYSIS_FILTER,
    SYNTHESIS_FILTER
} filter_type;

typedef struct pfb_filter_t
{
    double          *coeffs;
    int              ncoeffs;
    int              ntaps;
    int              nchans; // = size/ntaps
    cuDoubleComplex *twiddles; // twiddle factors
    filter_type      type;
} pfb_filter;



/* Function: ROOTS_OF_UNITY
 *
 * Creates a complex-valued array containing the N roots of unity.
 * The caller should free this memory (via free()).
 */
cuDoubleComplex *roots_of_unity( int N );



/* LOAD_FILTER_COEFFICIENTS
 *
 * Load a set of filter coefficients
 * Inputs:
 *   FILTERNAME - string specifying a filter. There should be a corresponding
 *                file in the RUNTIME_DIR called FILTERNAME.dat.
 *   TYPE       - Whether it's an ANALYSIS_FILTER or a SYNTHESIS_FILTER
 *   NCHANS     - the number of channels that this filter will be applied to.
 *                For both ANALYSIS and SYNTHESIS filters, this should be
 *                the number of ANALYSIS channels.
 * Outputs:
 *   [return value] - pointer to the newly allocated struct containing the
 *                    coefficients. This should be freed using
 *                    free_pfb_filter()
 */
pfb_filter *load_filter_coefficients( char *filtername, filter_type type, int nchans );



/* FREE_PFB_FILTER
 *
 * Free the memory allocated in load_filter_coefficients()
 */
void free_pfb_filter( pfb_filter *filter );


#endif
