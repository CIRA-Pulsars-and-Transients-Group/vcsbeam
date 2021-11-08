#!/usr/bin/env python3
# This can be run on the text files output by mwa_plot_calibration

import numpy as np
import matplotlib.pyplot as plt
import sys

def get_tilenames(filename, first_line=14, nlines=16):
    tilenames = []
    with open(filename, 'r') as reader:
        for i in range(first_line-1):
            reader.readline() # Consume the first several uninteresting lines
        for i in range(nlines):
            tilenames.extend(reader.readline().split()[1:])

    return tilenames

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("usage: {} [text file output from mwa_plot_calibration]")
        exit()

    cal_data = np.loadtxt(sys.argv[-1])

    nchans = cal_data.shape[0]
    nants = cal_data.shape[1]//8

    ncols = 16
    nrows = nants // ncols

    # Get tilenames
    tilenames = get_tilenames(sys.argv[-1])

    # Phases
    print("Plotting phases...")

    fig_phase, axs_phase = plt.subplots(nrows, ncols, sharex=True, sharey=True)

    s = 0.1

    for r in range(nrows):
        for c in range(ncols):
            ant = r*ncols + c

            axs_phase[r,c].scatter(np.arange(nchans), cal_data[:,ant*8 + 3], c='cyan', s=s, marker='.') # PQ phase
            axs_phase[r,c].scatter(np.arange(nchans), cal_data[:,ant*8 + 5], c='magenta', s=s, marker='.') # QP phase
            axs_phase[r,c].scatter(np.arange(nchans), cal_data[:,ant*8 + 1], c='b', s=s, marker='.') # PP phase
            axs_phase[r,c].scatter(np.arange(nchans), cal_data[:,ant*8 + 7], c='r', s=s, marker='.') # QQ phase

            axs_phase[r,c].set_xticklabels([])
            axs_phase[r,c].set_yticklabels([])
            axs_phase[r,c].set_ylim([-np.pi, np.pi])

            axs_phase[r,c].set_title(tilenames[ant], y=1.0, pad=-14)

    plt.subplots_adjust(wspace=0, hspace=0)
    fig_phase.suptitle("Phases\n[b c]\n[m r]")

    '''
    # Amps
    print("Plotting amps...")

    fig_amps, axs_amps = plt.subplots(nrows, ncols, sharex=True, sharey=True)

    for r in range(nrows):
        for c in range(ncols):
            ant = r*ncols + c

            axs_amps[r,c].scatter(np.arange(nchans), cal_data[:,ant*8 + 2], c='cyan', s=s, marker='.') # PQ amps
            axs_amps[r,c].scatter(np.arange(nchans), cal_data[:,ant*8 + 4], c='magenta', s=s, marker='.') # QP amps
            axs_amps[r,c].scatter(np.arange(nchans), cal_data[:,ant*8 + 0], c='b', s=s, marker='.') # PP amps
            axs_amps[r,c].scatter(np.arange(nchans), cal_data[:,ant*8 + 6], c='r', s=s, marker='.') # QQ amps

            axs_amps[r,c].set_xticklabels([])
            axs_amps[r,c].set_yticks([0, 1, 2])
            if c == 0:
                axs_amps[r,c].set_yticklabels(["0", "1", "2"])
            else:
                axs_amps[r,c].set_yticklabels([])

            axs_amps[r,c].set_title(tilenames[ant], y=1.0, pad=-14)

    plt.subplots_adjust(wspace=0, hspace=0)
    fig_amps.suptitle("Amps\n[b c]\n[m r]")
    '''

    plt.show()
