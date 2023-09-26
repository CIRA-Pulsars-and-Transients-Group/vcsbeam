#!/usr/bin/env python3

import sys
import os
import argparse
import numpy as np
from astropy.io import fits

def parse_args():
    parser = argparse.ArgumentParser(usage='Usage: %(prog)s [options]',
        description="Converts the 'flagged_tiles.txt' file employed in the " + \
            "RTS into a list of TileNames that can be used as input to " + \
            "make_mwa_tied_array_beam (see the -F option).")

    parser.add_argument('-m', '--metafits',
        action='store', type=str, dest='metafits_file', default=None,
        help='The calibration metafits file [default: %(default)s]')
    parser.add_argument('-i', '--input_flag_file',
        action='store', type=str, dest='input_flag_file', default=None,
        help='Input file containing a space-separated list of tile indices [default: %(default)s]')
    parser.add_argument('-o', '--output_flag_file',
        action='store', type=str, dest='output_flag_file', default=None,
        help='Where to write TileNames [default: stdout]')
    
    args = parser.parse_args()

    if args.metafits_file is None or args.input_flag_file is None:
        parser.error('--metafits and --input_flag_file are required arguments')

    # Check if input files exist
    if not os.path.exists(args.metafits_file):
        parser.error(f'cannot locate {args.metafits_file}')
    if not os.path.exists(args.input_flag_file):
        parser.error(f'cannot locate {args.input_flag_file}')
    
    # Check read access
    if not os.access(args.metafits_file, os.R_OK):
        parser.error(f'cannot open {args.input_flag_file}')
    if not os.access(args.input_flag_file, os.R_OK):
        parser.error(f'cannot open {args.input_flag_file}')

    # Check write access
    if args.output_flag_file is not None:
        output_dir = os.path.dirname(os.path.abspath(args.output_flag_file))
        if not os.access(output_dir, os.W_OK):
            parser.error(f'cannot write to directory {output_dir}')

    return args

if __name__ == '__main__':
    args = parse_args()

    # Get TileNames from metafits
    hdu = fits.open(args.metafits_file)
    tilenames = hdu[1].data['TileName']

    # Get flagged tiles
    flags = np.loadtxt(args.input_flag_file)
    if flags.ndim == 0:
        flags = np.array([flags])

    if args.output_flag_file is None:
        # Write TileNames to stdout
        for i in flags.astype(int):
            sys.stdout.write(f'{tilenames[i*2]}\n')
    else:
        # Write TileNames to file
        with open(args.output_flag_file, 'w') as outfile:
            for i in flags.astype(int):
                outfile.write(f'{tilenames[i*2]}\n')