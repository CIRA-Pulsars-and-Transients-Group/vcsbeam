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

def plot_tile(ax, PP, PQ, QP, QQ, phases_or_amps='phases', show_yticks=False, **kwargs):
    ax.scatter(np.arange(len(PP)), PQ, c='cyan', **kwargs)
    ax.scatter(np.arange(len(PP)), QP, c='magenta', **kwargs)
    ax.scatter(np.arange(len(PP)), PP, c='b', **kwargs)
    ax.scatter(np.arange(len(PP)), QQ, c='r', **kwargs)

    ax.set_xticklabels([])
    if phases_or_amps == 'phases':
        ax.set_ylim([-np.pi, np.pi])
        ax.set_yticks([-np.pi/2, 0, np.pi/2])
        if show_yticks:
            ax.set_yticklabels(["-π/2", "0", "π/2"])
        else:
            ax.set_yticklabels([])
    elif phases_or_amps == 'amps':
        ax.set_yticks([0, 1, 2])
        if show_yticks:
            ax.set_yticklabels(["0", "1", "2"])
        else:
            ax.set_yticklabels([])

def plot_all_tiles(cal_data, phases_or_amps='phases', tilenames=None, ncols=16, **kwargs):

    print("Plotting {}...".format(phases_or_amps))

    if phases_or_amps == 'phases':
        column_offset = 1
    elif phases_or_amps == 'amps':
        column_offset = 0

    nants = cal_data.shape[1]//8
    nrows = nants // ncols

    fig, axs = plt.subplots(nrows, ncols, sharex=True, sharey=True);

    for r in range(nrows):
        for c in range(ncols):
            ant = r*ncols + c

            show_yticks = (c == 0)

            PP = cal_data[:,ant*8 + 0 + column_offset]
            PQ = cal_data[:,ant*8 + 2 + column_offset]
            QP = cal_data[:,ant*8 + 4 + column_offset]
            QQ = cal_data[:,ant*8 + 6 + column_offset]

            plot_tile(axs[r,c], PP, PQ, QP, QQ, phases_or_amps=phases_or_amps, show_yticks=show_yticks, **kwargs)

            axs[r,c].set_title(tilenames[ant], y=1.0, pad=-14)

    plt.subplots_adjust(wspace=0, hspace=0)
    fig.suptitle("{}\n[b c]\n[m r]".format(phases_or_amps))

    return fig, axs


if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("usage: {} [text file output from mwa_plot_calibration]".format(sys.argv[0]))
        exit()

    cal_data = np.loadtxt(sys.argv[-1])

    # Get tilenames
    tilenames = get_tilenames(sys.argv[-1])

    s = 0.1

    plot_all_tiles(cal_data, phases_or_amps='phases', tilenames=tilenames, s=s, marker='.')
    plot_all_tiles(cal_data, phases_or_amps='amps', tilenames=tilenames, s=s, marker='.')

    plt.show()
