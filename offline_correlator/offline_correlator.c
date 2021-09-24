/* --------------------------- header secton ----------------------------*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#include "fourbit.h"
#include "xgpu.h"
#include "corr_utils.h"

/* SM: These should not be defined here, but should be taken directly from
 * VCSTOOLS. However, VCSTOOLS does not currently have any facility to export
 * this kind of info externally.
 */
#define NFREQUENCY   128ul
#define NSTATION     128ul
#define NPOL           2ul
#define NBIT           4ul
#define NTIMESTEPS 10000ul
#define NDIM           2ul  /* SM: re+im = 2 parts to complex number */

#define NSAMPLES       (NTIMESTEPS*NFREQUENCY*NSTATION*NPOL*NDIM)
#define NCMPLX_SAMPLES (NTIMESTEPS*NFREQUENCY*NSTATION*NPOL)
#define NSAMPLES_PER_TIMESTEP  (NSAMPLES/NTIMESTEPS)

/* SM: This is apparently a fixed HDU header size for the output gpubox files
 */
#define HEADER_SIZE 2880

/*

 Data ordering for input vectors is (running from slowest to fastest)
 [time][channel][station][polarization][complexity]

 SM: This (^) will ultimately be handled by VCSTOOLS

 Output matrix has ordering
 [channel][station][station][polarization][polarization][complexity]

 SM: This (^) is currently the only place I know of that handles WRITING to
     GPUBox files. Once I understand how this works, I will move to export
     this functionality into libmwa, which already handles READING gpubox
     files

 */

/*--------------------------------------------------------------------------*/

//get command options
typedef struct Options_t {
    int coarse_chan;
    int edge;
    int chan_to_aver;
    time_t starttime;
    char *in_file;
    char *out_file;
    char *obsid;
    int offline;
    int dumps_per_second;
} Options;

void usage()
{
    printf( "offline_correlator: a light-weight correlator for the MWA. "
            "Takes VCS data files and correlates as per the "
            "parameters of the linked xGPU library\n" );
    printf( "usage: offline_correlator -c <coarse_channel> -d <infile> [options]\n" );
    printf( "Options:\n" );
    printf( "   -e EDGES\n" );
    printf( "        Set EDGES channels at both top and bottom of the band to 0 "
            "[default: 0]\n" );
    printf( "   -h\n" );
    printf( "        Display this help and exit\n" );
    printf( "   -n CHAN_AVERAGE\n" );
    printf( "        Average CHAN_AVERAGE adjacent channels in the final output "
            "[default: 4]\n" );
    printf( "   -o OBSID\n" );
    printf( "        The observation ID of the input VCS data "
            "[required]\n" );
    printf( "   -r DUMPS_PER_SECOND\n" );
    printf( "        The number of correlator dumps per second to write to file "
            "[default: 1]\n" );
    printf( "   -s STARTTIME\n" );
    printf( "        The time (in Unix seconds) corresponding to the input file"
            "[required]\n" );
}

void parse_cmdline(int argc, char *argv[], Options *opt)
{
    if (argc == 1)
    {
        usage();
        exit(EXIT_FAILURE);
    }

    int arg = 0;
    while ((arg = getopt(argc, argv, "c:d:e:hn:o:r:s:")) != -1) {

        switch (arg) {
            case 'c':
                opt->coarse_chan = atoi(optarg);
                break;
            case 'd':
                opt->in_file = strdup(optarg);
                break;
            case 'e':
                opt->edge = atoi(optarg);
                break;
            case 'h':
                usage();
                exit(EXIT_SUCCESS);
            case 'n':
                opt->chan_to_aver=atoi(optarg);
                break;
            case 'o':
                opt->obsid = strdup(optarg);
                break;
            case 'r':
                opt->dumps_per_second = atoi(optarg);
                break;
            case 's':
                opt->starttime = (time_t) atol(optarg);
                break;
            default:
                usage();
                exit(EXIT_FAILURE);
                break;
        }
    }

    if (opt->in_file == NULL)
    {
        usage();
        fprintf( stderr, "error: no input file given\n" );
        exit(EXIT_FAILURE);
    }

    if (opt->coarse_chan < 0)
    {
        usage();
        fprintf( stderr, "error: a positive coarse channel must be given\n" );
        exit(EXIT_FAILURE);
    }

    if (opt->out_file == NULL && opt->starttime < 0)
    {
        usage();
        fprintf( stderr, "error: starttime required\n" );
        exit(EXIT_FAILURE);
    }

    if (opt->obsid == NULL)
    {
        usage();
        fprintf(stderr, "error: obsid required\n" );
        exit(EXIT_FAILURE);
    }

    // If no explicit out_file is given, a default one is constructed:
    // [obsid]_[timestamp]_gpubox[chan]_00.fits

    // Convert the obsid and starttime into an output file name

    // Convert the start time to a string
    struct tm current_utctime;
    char file_time[16];
    gmtime_r( &opt->starttime, &current_utctime );
    strftime( file_time, 15, "%Y%m%d%H%M%S", &current_utctime );

    // Write the file name to opt->out_file
    opt->out_file = (char *)malloc( 128 ); // should be more than enough
    sprintf( opt->out_file, "%s_%s_gpubox%02d_00.fits",
            opt->obsid, file_time, opt->coarse_chan );

}


void read_vcs( FILE *fin, size_t ntime, int edge, ComplexInput *out )
/* Read VCS data files and convert from (4+4)-bit complex to
 * (8+8)-bit complex.
 */
{
    // Check that the out pointer has been allocated
    if (out == NULL)
    {
        fprintf( stderr, "error: read_vcs: out pointer is NULL\n" );
        exit(EXIT_FAILURE);
    }

    // SM: TODO: Put this elsewhere?
    build_eight_bit_lookup();

    // Read data in from file
    size_t bytes_to_read = ntime*NSAMPLES_PER_TIMESTEP*NBIT/8;
    char *vcsdata = (char *)malloc( bytes_to_read );
    size_t nread = fread( vcsdata, 1, bytes_to_read, fin );
    if (nread != bytes_to_read)
    {
        fprintf( stderr, "error: read only "
                "%ld of %ld bytes\n", nread, bytes_to_read );
        exit(EXIT_FAILURE);
    }

    // Loop through the VCS data, convert from 4-bit VCS format to
    // int8_t samples, and write result to "out"

    size_t real_samples_per_step = NPOL*NDIM; // = 4
    size_t s; // (s)ample (advances NPOL*NDIM*NBIT = 16 bits = 2 bytes at a time)
    size_t c; // (c)hannel number (for determining if it's an edge channel)
    size_t r; // iterate over (r)eal_samples_per_step
    ComplexInput *out_ptr = out; // Loop through output array via this pointer
    int8_t expanded[4]; // Temp variable for doing the 4-bit expansion
    size_t samps_per_chan = (NSTATION*NPOL*NDIM*NBIT)/8;

    for (s = 0; s < nread; s += (NPOL*NDIM*NBIT)/8)
    {
        // Determine the channel number of this sample
        c = (s / samps_per_chan) % NFREQUENCY;

        if ((c < (size_t)edge ) || c >= (NFREQUENCY - edge))
        {
            // Set edge channels to 0
            for (r = 0; r < real_samples_per_step; r++)
                expanded[r] = 0;
        }
        else
        {
            // Convert 2x pols:
            // ((4+4)*2 = 1*16 bits) --> ((8+8)*2 = 32 = 4*8 bits)
            expand_4bit( (uint16_t *)&vcsdata[s], expanded );
        }

        // Pack the answer into the output array and
        // advance 2 complex samples in the output array
        out_ptr[0].real = (ReImInput)expanded[0];
        out_ptr[0].imag = (ReImInput)expanded[1];
        out_ptr[1].real = (ReImInput)expanded[2];
        out_ptr[1].imag = (ReImInput)expanded[3];
        out_ptr += 2;
    }

    free( vcsdata );
}

int main(int argc, char **argv)
{
    clock_t start = clock();
    printf( "[%9.5lf] Starting offline correlator\n", 0.0 );

    // Set default values for command line options
    Options opt;

    opt.in_file          = NULL;
    opt.obsid            = NULL;
    opt.starttime        = -1;
    opt.chan_to_aver     = 4; // number of channels to combine on output
    opt.edge             = 0;
    opt.coarse_chan      = -1; // only set in the header if this is >= 0
    opt.out_file         = NULL;
    opt.dumps_per_second = 1;

    // Parse the command line
    parse_cmdline( argc, argv, &opt );

    // Prepare structs/variables for xGPU-related info
    XGPUInfo xgpu_info;
    int xgpu_error = 0;

    // Open the input file for reading
    FILE *fin = fopen( opt.in_file, "r" );
    if (fin == NULL)
    {
        fprintf( stderr, "error: unable to open '%s' for reading\n",
                opt.in_file );
        exit(EXIT_FAILURE);
    }

    // Open the out_file for writing
    FILE *fout = fopen( opt.out_file, "w" );
    if (fout == NULL)
    {
        fprintf( stderr, "error: failed to open gpubox file '%s' "
                "for writing\n", opt.out_file );
        exit(EXIT_FAILURE);
    }

    // Get sizing info from library, and make sure that xGPU was compiled with the
    // settings required for processsing VCS data

    xgpuInfo(&xgpu_info);

    if ((xgpu_info.npol       != NPOL      ) ||
        (xgpu_info.nstation   != NSTATION  ) ||
        (xgpu_info.nfrequency != NFREQUENCY) ||
        (xgpu_info.input_type != XGPU_INT8 ))
    {
        fprintf( stderr, "error: XGPU library was compiled with:\n"
                "\tNPOL       = %d\n"
                "\tNSTATION   = %d\n"
                "\tNFREQUENCY = %d\n"
                "\tINPUT_TYPE = %d\n"
                "but offline_correlator requires:\n"
                "\tNPOL       = %lu\n"
                "\tNSTATION   = %lu\n"
                "\tNFREQUENCY = %lu\n"
                "\tINPUT_TYPE = %d\n",
                xgpu_info.npol,
                xgpu_info.nstation,
                xgpu_info.nfrequency,
                xgpu_info.input_type,
                NPOL,
                NSTATION,
                NFREQUENCY,
                0 );
        exit(EXIT_FAILURE);
    }

    // Also, double check that the VCS file sizes are an
    // integral number of xGPU's vecLength
    if (NCMPLX_SAMPLES % xgpu_info.vecLength)
    {
        fprintf( stderr, "error: The number of complex samples in the VCS "
                "input files (%lu) is not an integral number of the vector "
                "size that xGPU is expecting (%llu)\n", NCMPLX_SAMPLES,
                xgpu_info.vecLength );
        exit(EXIT_FAILURE);
    }
    int nintegrations = NCMPLX_SAMPLES / xgpu_info.vecLength / opt.dumps_per_second;
    fprintf( stderr, "nintegrations = %d / %d / %d = %d\n", NCMPLX_SAMPLES, xgpu_info.vecLength, opt.dumps_per_second, nintegrations );

    // Allocate output matrix array
    size_t full_matLength = NFREQUENCY * NSTATION * NSTATION * NPOL * NPOL;
    size_t full_size_bytes = opt.dumps_per_second * full_matLength * sizeof(Complex);
    Complex *full_matrix_h = (Complex *)malloc( full_size_bytes );

    // The PRIMARY header + image + padding
    // SM: This seems to be relevant to the output FITS format gpubox files.
    uint64_t n_visibilities = ((uint64_t) xgpu_info.nbaseline * 4); // SM: Why *4?
    uint64_t blockSize = 0;

    // Work out the blockSize for the output gpubox files
    int hdu_num;
    for (hdu_num = 0; hdu_num < opt.dumps_per_second; hdu_num++) {
        blockSize += HEADER_SIZE;
        blockSize += n_visibilities * (uint64_t)NFREQUENCY *
            sizeof(Complex) / opt.chan_to_aver; // sizeof a data cube

        int remainder = blockSize % HEADER_SIZE; // pad out to the end of the HDU
        blockSize += HEADER_SIZE - remainder; // SM: why an extra HEADER_SIZE?
    }
    if ((blockSize % HEADER_SIZE) != 0)
    {
        fprintf( stderr, "error: Unable to construct correct block size for FITS "
                "output (blockSize = %lu is not divisible by %d)\n",
                blockSize, HEADER_SIZE );
        exit(EXIT_FAILURE);
    }

    // Initialise xGPU
    printf( "[%9.5lf] Initialising xGPU\n", (clock()-start)/(double)CLOCKS_PER_SEC );
    XGPUContext context;
    context.array_h    = NULL;
    context.matrix_h   = NULL;
    context.array_len  = xgpu_info.vecLength;
    context.matrix_len = xgpu_info.matLength;

    xgpu_error = xgpuInit( &context, 0 );
    if(xgpu_error)
    {
        fprintf( stderr, "error: xgpuInit returned error code %d\n", xgpu_error );
        exit(EXIT_FAILURE);
    }

    /* -------------------------- Run xGPU to correlate -------------------- */

    int d; // iterate over (d)umps_per_second
    int i; // iterate over (i)ntegrations
    int ntimesteps = NTIMESTEPS/(nintegrations*opt.dumps_per_second);
    for (d = 0; d < opt.dumps_per_second; d++)
    {
        for (i = 0; i < nintegrations; i++)
        {
            // Read in the next chunk of VCS data
            printf( "[%9.5lf] Reading in next %llu bytes\n",
                    (clock()-start)/(double)CLOCKS_PER_SEC,
                    ntimesteps*NSAMPLES_PER_TIMESTEP*NBIT/8 );
            read_vcs( fin, ntimesteps, opt.edge, context.array_h );

            // Report progress so far
            printf( "[%9.5lf] Running GPU X-Engine (%d/%d)\n",
                    (clock()-start)/(double)CLOCKS_PER_SEC, i, nintegrations );

            // Run the correlator for this integration
            xgpu_error = xgpuCudaXengine( &context, SYNCOP_DUMP );
            if(xgpu_error)
            {
                fprintf( stderr, "\nerror: xgpuCudaXengine returned error "
                        "code %d\n", xgpu_error );
                exit(EXIT_FAILURE);
            }

            fflush(stdout);
        }

        // Copy the result out into the full matrix array
        Complex *ptr = full_matrix_h + (d*context.matrix_len);
        memcpy( ptr, context.matrix_h, (context.matrix_len * sizeof(Complex)) );

        // Clear the xGPU integration array on device
        xgpu_error = xgpuClearDeviceIntegrationBuffer( &context );
        if(xgpu_error)
        {
            fprintf( stderr, "error: xgpuClearDeviceIntegrationBuffer "
                    "returned error code %d\n", xgpu_error );
            exit(EXIT_FAILURE);
        }
    }
    
    // The whole second of VCS data should now have been correlated and packed
    // into full_matrix_h

    /* -------------------------- Write out to gpubox -------------------- */


    printf( "[%9.5lf] Building FITS file\n", (clock()-start)/(double)CLOCKS_PER_SEC );

    manager_t the_manager;
    the_manager.shutdown      = 0;
    the_manager.offline       = opt.offline;
    the_manager.integrate     = 0;
    the_manager.chan_to_aver  = opt.chan_to_aver;
    the_manager.edge          = opt.edge;
    the_manager.nbit          = NBIT;
    the_manager.coarse_chan   = opt.coarse_chan;
    the_manager.nstation      = NSTATION;
    the_manager.nfrequency    = NFREQUENCY;
    the_manager.ndim          = NDIM;
    the_manager.npol          = NPOL;
    the_manager.dumps_per_sec = opt.dumps_per_second;
    the_manager.infile        = 0;
    the_manager.marker        = 0;

    char *FITSbuffer = (char *) malloc(blockSize);

    buildFITSBuffer( xgpu_info,
            full_matrix_h,
            blockSize,
            (void *)FITSbuffer,
            opt.starttime,
            opt.dumps_per_second,
            &the_manager );

    fwrite(FITSbuffer, blockSize, 1, fout);

    printf( "[%9.5lf] Finished. Output written to '%s'\n",
            (clock()-start)/(double)CLOCKS_PER_SEC, opt.out_file );

    /* ------------------ Close files and free memory ------------------- */
    xgpuFree( &context );
    free( full_matrix_h );
    free( FITSbuffer );
    fclose( fin );
    fclose( fout );
    free( opt.out_file );
    if (opt.in_file)
        free( opt.in_file );

    return EXIT_SUCCESS;
}
