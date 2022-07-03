#!/usr/bin/env python
import sys
import numpy as np
import shapefile
from shapely.geometry import Polygon
from optparse import OptionParser,IndentedHelpFormatter

# Default values
RAT_MIN = 0.9999

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('-s','--shp_fnam',default=None,help='Shape file name (%default)')
parser.add_option('-o','--outline_fnam',default=None,help='Outline shapefile name (%default)')
parser.add_option('-m','--rat_min',default=RAT_MIN,type='float',help='Mininum ratio (%default)')
parser.add_option('--use_index',default=False,action='store_true',help='Use index instead of OBJECTID (%default)')
(opts,args) = parser.parse_args()

r = shapefile.Reader(opts.shp_fnam)
o = shapefile.Reader(opts.outline_fnam)
polygons = []
for outline in o.iterShapes():
    polygons.append(Polygon(outline.points))

for ii,shaperec in enumerate(r.iterShapeRecords()):
    if ii%1000 == 0:
        sys.stderr.write('{}\n'.format(ii))
        sys.stderr.flush()
    rec = shaperec.record
    shp = shaperec.shape
    if opts.use_index:
        object_id = ii+1
    else:
        object_id = rec.OBJECTID
    if len(shp.points) < 1:
        sys.stderr.write('Warning, len(shp.points)={}, ii={}\n'.format(len(shp.points),ii))
        continue
    poly_r = Polygon(shp.points)
    rat = 0.0
    err = False
    for poly_o in polygons:
        try:
            poly_intersect = poly_r.intersection(poly_o)
        except Exception:
            sys.stderr.write('Warning, error occured at ii={}\n'.format(ix,iy,ii))
            err = True
            break
        rat += poly_intersect.area/poly_r.area
    if err:
        continue
    if rat < opts.rat_min:
        sys.stdout.write('{}\n'.format(object_id))
