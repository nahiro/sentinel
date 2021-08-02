#!/usr/bin/env python
import os
import sys
import shutil
import gdal
import numpy as np
from datetime import datetime
import shapefile
from matplotlib.dates import date2num,num2date
from statsmodels.stats.weightstats import DescrStatsW
from optparse import OptionParser,IndentedHelpFormatter

# Default values
HOME = os.environ.get('HOME')
if HOME is None:
    HOME = os.environ.get('HOMEPATH')
INPNAM = os.path.join(HOME,'Work','SATREPS','Shapefile','field_GIS','cihea','New_Test_Sites')
OUTNAM = os.path.join('.','transplanting_date')
DATA_FILE = 'output.tif'
AREA_FILE = 'pixel_area_block.dat'

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('--data_file',default=DATA_FILE,help='Estimation data file (%default)')
parser.add_option('--area_file',default=AREA_FILE,help='Area data file (%default)')
parser.add_option('--inpnam',default=INPNAM,help='Input shapefile (%default)')
parser.add_option('--outnam',default=OUTNAM,help='Output shapefile (%default)')
parser.add_option('-d','--debug',default=False,action='store_true',help='Debug mode (%default)')
(opts,args) = parser.parse_args()

ds = gdal.Open(opts.data_file)
data = ds.ReadAsArray()
ds = None
data = data.reshape(len(data),-1)

object_ids = []
blocks = []
inds = []
areas = []
with open(opts.area_file,'r') as fp:
    for line in fp:
        item = line.split()
        if len(item) < 5 or item[0] == '#':
            continue
        object_ids.append(int(item[0]))
        blocks.append(item[1])
        inds.append([])
        areas.append([])
        n = int(item[2])
        for nn in range(3,n*2+3,2):
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
    data_value = data[0,inds[i]]
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

np.save(opts.out_fnam,object_ids=object_ids,data_avg=data_avg,data_std=data_std)
