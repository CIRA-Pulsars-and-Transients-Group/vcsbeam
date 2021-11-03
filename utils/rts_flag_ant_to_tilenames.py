#!/usr/bin/env python3

import sys
from astropy.io import fits
import numpy as np

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Converts the 'flagged_tiles.txt' file employed in the RTS into a list of TileNames that can be used as input to make_mwa_tied_array_beam (see the -F option).")
        print("usage: {}  [calibration metafits file]  [flagged_tiles.txt]")

    metafits_file = sys.argv[1]
    flagged_file = sys.argv[2]

    hdu = fits.open(metafits_file)

    tilenames = hdu[1].data['TileName']

    flags = np.loadtxt(flagged_file)

    for i in flags.astype(int):
        print(tilenames[i*2])

