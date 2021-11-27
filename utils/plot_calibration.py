#!/usr/bin/env python3
# This can be run on the text files output by mwa_plot_calibration

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import argparse

def get_header_info(filename, header_nlines=30):

    tilenames_first_line = 14
    tilenames_nlines     = 16
    tilenames            = []

    mwa_plot_calibration_line = 4

    with open(filename, 'r') as reader:
        for i in range(header_nlines):
            line = reader.readline() # Consume the next line

            if i+1 == mwa_plot_calibration_line:
                mwa_plot_calibration_cmd = line.strip(' #')
            elif i+1 >= tilenames_first_line and i+1 < tilenames_first_line + tilenames_nlines:
                tilenames.extend(line.split()[1:])

    return tilenames, mwa_plot_calibration_cmd

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

def plot_all_tiles(cal_data, mwa_plot_calibration_cmd, phases_or_amps='phases', tilenames=None, ncols=16, figsize=None, **kwargs):

    if phases_or_amps == 'phases':
        column_offset = 1
    elif phases_or_amps == 'amps':
        column_offset = 0

    nants = cal_data.shape[1]//8
    nrows = nants // ncols

    fig, axs = plt.subplots(nrows, ncols, sharex=True, sharey=True, figsize=figsize);

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

    # Put a legend on the top right
    custom_lines = [Line2D([0], [0], color='w', markerfacecolor="blue",    marker='o', label="Dqq"),
                    Line2D([0], [0], color='w', markerfacecolor="magenta", marker='o', label="Dpq"),
                    Line2D([0], [0], color='w', markerfacecolor="cyan",    marker='o', label="Dqp"),
                    Line2D([0], [0], color='w', markerfacecolor="red",     marker='o', label="Dpp")]
    axs[0,-1].legend(handles=custom_lines,
                     loc='upper center',
                     bbox_to_anchor=(0.5, 0.95),
                     bbox_transform=fig.transFigure,
                     ncol=2,
                     labelspacing=0.5,
                     columnspacing=0.0)

    # Squish the plots together
    plt.subplots_adjust(wspace=0, hspace=0, left=0.01, right=0.99, top=0.89, bottom=0.01)

    # Add a global title
    fig.suptitle("{}\n{}".format(phases_or_amps, mwa_plot_calibration_cmd))

    return fig, axs


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Make phase and amplitude plots of calibration solutions")
    parser.add_argument("calibration_solution", help="file for plotting (output of mwa_plot_calibration)")
    parser.add_argument("--phases_png", help="output plot of phases", default=None)
    parser.add_argument("--amps_png", help="output plot of amps", default=None)
    args = parser.parse_args()

    cal_data = np.loadtxt(args.calibration_solution)

    # Get tilenames
    tilenames, mwa_plot_calibration_cmd = get_header_info(args.calibration_solution)

    s = 0.1
    figsize = [12.8, 9.6]

    if args.phases_png is not None:
        print("Plotting phases...")
        fig, axs = plot_all_tiles(cal_data, mwa_plot_calibration_cmd, phases_or_amps='phases', tilenames=tilenames, s=s, marker='.', figsize=figsize)
        plt.savefig(args.phases_png)

    if args.amps_png is not None:
        print("Plotting amps...")
        fig, axs = plot_all_tiles(cal_data, mwa_plot_calibration_cmd, phases_or_amps='amps', tilenames=tilenames, s=s, marker='.', figsize=figsize)
        plt.savefig(args.amps_png)
