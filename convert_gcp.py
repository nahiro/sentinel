#!/usr/bin/env python
import os
import sys
import re
import gdal
import numpy as np
from optparse import OptionParser,IndentedHelpFormatter

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('-i','--src_fnam',default=None,help='Source file name (%default)')
parser.add_option('-o','--dst_fnam',default=None,help='Destination file name (%default)')
parser.add_option('-I','--src_geotiff',default=None,help='Source GeoTIFF name (%default)')
parser.add_option('-e','--exp',default=False,action='store_true',help='Output in exp format (%default)')
parser.add_option('--long',default=False,action='store_true',help='Output in long format (%default)')
(opts,args) = parser.parse_args()

ds = gdal.Open(opts.src_geotiff)
src_trans = ds.GetGeoTransform()
if src_trans[2] != 0.0 or src_trans[4] != 0.0:
    raise ValueError('Error, src_trans={}'.format(src_trans))
src_xmin = src_trans[0]
src_xstp = src_trans[1]
src_ymax = src_trans[3]
src_ystp = src_trans[5]
ds = None

src_xi = []
src_yi = []
src_line = []
with open(opts.src_fnam,'r') as fp:
    for line in fp:
        #660.5    120.5 751655.112266 9243156.957321  -0.009634  -0.003979    2.10494      0.824779      0.016902   3555      0.789201   1440
        #720.5    120.5 751667.515378 9243156.775077   0.393478  -0.186223    2.37634      0.832486      0.014921   5446      0.797028   1845
        m = re.search('^\s*(\S+)\s+(\S+)\s+(\S.*)$',line)
        if not m:
            raise ValueError('Error in reading '+opts.src_fnam)
        src_xi.append(float(m.group(1)))
        src_yi.append(float(m.group(2)))
        src_line.append(m.group(3))
src_xi = np.array(src_xi)
src_yi = np.array(src_yi)

dst_xi = src_xmin+src_xstp*src_xi
dst_yi = src_ymax+src_ystp*src_yi
with open(opts.dst_fnam,'w') as fp:
    for xi,yi,line in zip(dst_xi,dst_yi,src_line):
        if opts.exp:
            fp.write('{:15.8e} {:15.8e} {}\n'.format(xi,yi,line))
        elif opts.long:
            fp.write('{:12.6f} {:12.6f} {}\n'.format(xi,yi,line))
        else:
            fp.write('{:8.2f} {:8.2f} {}\n'.format(xi,yi,line))
