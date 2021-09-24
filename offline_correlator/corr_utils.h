#ifndef __CORR_UTILS_H
#define __CORR_UTILS_H

#include "xgpu.h"


typedef struct management {

        int data_line[8];
        int data_line_count;
        char start_obs_id[64];
        char start_obs_UTC[256];
        int new_obs;
        int offline;
        int integrate; // how many seconds to integrate
        int chan_to_aver;     // how many channels to average
        int marker;
        int shutdown;
        int nfrequency;
        int ntime;
        int nstation;
        int ndim;
        int npol;
        int edge;
        int nbit;
        int ncoarse;
        int coarse_chan;
        int dumps_per_sec;
	int infile;

} manager_t;

void printerror(int);
void buildFITSBuffer(XGPUInfo,Complex *,size_t,void *,time_t,int, manager_t *);

#endif
