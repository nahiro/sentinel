#!/usr/bin/env python
import os
import sys
import numpy as np
import shapefile
from shapely.geometry import Point,Polygon
from optparse import OptionParser,IndentedHelpFormatter

# Default values
RMAX = 1.0

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.set_usage('Usage: %prog shapefile_for_checking [options]')
parser.add_option('-r','--rmax',default=RMAX,type='float',help='Maximum distance in meter to be detected (%default)')
(opts,args) = parser.parse_args()
if len(args) < 1:
    parser.print_help()
    sys.exit()
shpnam = args[0]

np.warnings.filterwarnings('ignore')

xs = []
ys = []
ps = []
r = shapefile.Reader(shpnam)
#r = shapefile.Reader(shpnam,encoding='cp1252') # cp1252 for Windows
for i,shaperec in enumerate(r.iterShapeRecords()):
    rec = shaperec.record
    shp = shaperec.shape
    # Check if the polygon is valid
    p1 = Polygon(shp.points)
    if not p1.is_valid:
        sys.stderr.write('Error, FID {:6d} is not valid\n'.format(i))
    # Check if the barycenter and the center of bbox is the same -> not needed
    pp = np.array(shp.points)
    if p1.is_empty:
        sys.stderr.write('Error, FID {:6d} is empty\n'.format(i))
        xc = np.nan
        yc = np.nan
    else:
        xc = pp[:,0].mean()
        yc = pp[:,1].mean()
    #x1,y1,x2,y2 = shp.bbox
    #xctr = 0.5*(x1+x2)
    #yctr = 0.5*(y1+y2)
    #l2 = np.square(xc-xctr)+np.square(yc-yctr)
    #if l2 > 1.0e-4:
    #    sys.stderr.write('Error, FID {} l2 = {}\n'.format(i,l2))
    xs.append(xc)
    ys.append(yc)
    ps.append(pp)
xs = np.array(xs)
ys = np.array(ys)

n = np.arange(xs.size)
# Check if the polygon is unique
for i in range(xs.size):
    r = np.sqrt(np.square(xs-xs[i])+np.square(ys-ys[i]))
    indx = n[(~np.isnan(r)) & (r < opts.rmax) & (n != i)]
    if len(indx) != 0:
        sys.stderr.write('Warning, FID={:6d}'.format(i))
        for j in indx:
            sys.stderr.write(', r[{:6d}]={:13.6e}, all_equal={}'.format(j,r[j],np.all(ps[i]==ps[j])))
            xs[j] = np.nan
        sys.stderr.write('\n')
