import sys
from astropy.io import fits

hdu = fits.open(sys.argv[1])

rows = zip(hdu[1].data["Input"], hdu[1].data["Antenna"], hdu[1].data["Tile"], hdu[1].data["TileName"], hdu[1].data["Pol"], hdu[1].data["VCSOrder"])

print("# Input | Antenna | Tile | TileName | Pol | VCSOrder")
for r in rows:
    print("{:3d} {:3d} {:5d} {:8s} {:1s} {:3d}".format(*r))

