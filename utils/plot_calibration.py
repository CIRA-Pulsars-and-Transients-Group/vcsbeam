# This can be run on the text files output by mwa_plot_calibration

import numpy as np
import matplotlib.pyplot as plt
import sys

cal_data = np.loadtxt(sys.argv[-1])

nchans = cal_data.shape[0]
nants = cal_data.shape[1]//8

ncols = 16
nrows = nants // ncols

fig_phase, axs = plt.subplots(nrows, ncols, sharex=True, sharey=True)

s = 0.2

for r in range(nrows):
    for c in range(ncols):
        ant = r*ncols + c

        #axs[r,c].scatter(np.arange(nchans), cal_data[:,ant*8 + 3], c='cyan', s=s, marker='.') # PQ phase
        #axs[r,c].scatter(np.arange(nchans), cal_data[:,ant*8 + 5], c='magenta', s=s, marker='.') # QP phase
        axs[r,c].scatter(np.arange(nchans), cal_data[:,ant*8 + 1], c='b', s=s, marker='.') # PP phase
        axs[r,c].scatter(np.arange(nchans), cal_data[:,ant*8 + 7], c='r', s=s, marker='.') # QQ phase

        axs[r,c].set_xticklabels([])
        axs[r,c].set_yticklabels([])
        axs[r,c].set_ylim([-np.pi, np.pi])

plt.subplots_adjust(wspace=0, hspace=0)
plt.show()
