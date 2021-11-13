#!/usr/bin/env python3

import sys
from astropy.io import fits
import numpy as np

if len(sys.argv) < 3:
    print("usage: {} [columnname1, columnname2, ...] [metafits file]".format(sys.argv[0]))
    exit(0)

hdu = fits.open(sys.argv[-1])

fields = sys.argv[1:-1]
columns = np.array([hdu[1].data[field] for field in fields])

rows = zip(*columns)

print("# " + " | ".join(fields))
for r in rows:
    strings = ["{}".format(item) for item in list(r)]
    print(" ".join(strings))

