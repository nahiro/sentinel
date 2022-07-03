#!/usr/bin/env python
import os
import sys
import shutil
try:
    import gdal
except Exception:
    from osgeo import gdal
import numpy as np
from datetime import datetime
import shapefile
from matplotlib.dates import date2num,num2date
from statsmodels.stats.weightstats import DescrStatsW
from optparse import OptionParser,IndentedHelpFormatter

# Default values
OUT_FNAM = 'output.npz'
AREA_FNAM = 'pixel_area_block.dat'
BAND = 0

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('-i','--inp_fnam',default=None,help='Input file name (%default)')
parser.add_option('-o','--out_fnam',default=OUT_FNAM,help='Output npz file name (%default)')
parser.add_option('--area_fnam',default=AREA_FNAM,help='Pixel area file name (%default)')
parser.add_option('-b','--band',default=BAND,type='int',help='Target band# (%default)')
parser.add_option('-n','--no_block',default=False,action='store_true',help='No block in pixel area file (%default)')
parser.add_option('-d','--debug',default=False,action='store_true',help='Debug mode (%default)')
(opts,args) = parser.parse_args()

ds = gdal.Open(opts.inp_fnam)
if ds.RasterCount < 2:
    data = ds.ReadAsArray().flatten()
else:
    band = ds.GetRasterBand(opts.band+1)
    data = band.ReadAsArray().flatten()
ds = None

object_ids = []
blocks = []
inds = []
areas = []
if opts.no_block:
    n0 = 2
else:
    n0 = 3
with open(opts.area_fnam,'r') as fp:
    for line in fp:
        item = line.split()
        nitem = len(item)
        if nitem < n0 or item[0] == '#':
            continue
        n = int(item[n0-1])
        if nitem != n*2+n0:
            raise ValueError('Error, nitem={}, n0={}, n={}'.format(nitem,n0,n))
        object_ids.append(int(item[0]))
        if not opts.no_block:
            blocks.append(item[1])
        inds.append([])
        areas.append([])
        for nn in range(n0,nitem,2):
            inds[-1].append(int(item[nn]))
            areas[-1].append(float(item[nn+1]))
        inds[-1] = np.array(inds[-1])
        areas[-1] = np.array(areas[-1])
object_ids = np.array(object_ids)
blocks = np.array(blocks)
inds = np.array(inds,dtype='object')
areas = np.array(areas,dtype='object')

data_avg = []
data_std = []
for i in range(object_ids.size):
    object_id = object_ids[i]
    block = blocks[i]
    data_value = data[inds[i]]
    data_weight = areas[i]
    cnd = ~np.isnan(data_value)
    if cnd.sum() <= 1:
        data_avg.append(data_value[cnd].mean())
        data_std.append(0.0)
    else:
        data_weighted_stats = DescrStatsW(data_value[cnd],weights=data_weight[cnd],ddof=0)
        data_avg.append(data_weighted_stats.mean)
        data_std.append(data_weighted_stats.std)
data_avg = np.array(data_avg)
data_std = np.array(data_std)

np.savez(opts.out_fnam,object_ids=object_ids,data_avg=data_avg,data_std=data_std)
