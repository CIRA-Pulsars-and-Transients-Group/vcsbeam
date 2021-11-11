#!/usr/bin/env python3

import sys
from astropy.io import fits
import numpy as np

if len(sys.argv) < 3:
    print("usage: {} [columnname1, columnname2, ...] [metafits file]")
    exit(0)

hdu = fits.open(sys.argv[-1])

fields = sys.argv[1:-1]
columns = np.array([hdu[1].data[field] for field in fields])

rows = zip(*columns)

print("# " + " | ".join(fields))
for r in rows:
    print(" ".join(list(r)))

