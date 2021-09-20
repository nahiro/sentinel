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
parser.add_option('--true_file',default=None,help='True data file (%default)')
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

scheduled = {}
if opts.true_file is not None:
    with open(opts.true_file,'r') as fp:
        for line in fp:
            item = line.split()
            scheduled.update({item[0]:date2num(datetime.strptime(item[1],'%Y/%m/%d'))})

data_dict = {}
for i in range(object_ids.size):
    object_id = object_ids[i]
    block = blocks[i]
    if block in scheduled:
        data_scheduled = scheduled[block]
    else:
        data_scheduled = np.nan
    data_value = data[0,inds[i]]
    signal_value = data[1,inds[i]]
    data_weight = signal_value
    #signal_weight = areas[i]
    cnd = ~np.isnan(data_value)
    if cnd.sum() <= 1:
        data_avg = data_value[cnd].mean()
        signal_avg = signal_value[cnd].mean()
        data_std = 0.0
        signal_std = 0.0
    else:
        data_weighted_stats = DescrStatsW(data_value[cnd],weights=data_weight[cnd],ddof=0)
        data_avg = data_weighted_stats.mean
        data_std = data_weighted_stats.std
        signal_avg = signal_value[cnd].mean()
        signal_std = signal_value[cnd].std()
    data_dict.update({object_id:[block,data_scheduled,data_avg,signal_avg]})

r = shapefile.Reader(opts.inpnam)
w = shapefile.Writer(opts.outnam)
w.shapeType = shapefile.POLYGON
w.fields = r.fields[1:] # skip first deletion field
w.field('Block','C',6,0)
w.field('TANAM','F',13,6)
w.field('TANAM_text','C',10,0)
w.field('trans_d','F',13,6)
w.field('trans_t','C',10,0)
w.field('trans_s','F',13,6)
w.field('trans_n','F',13,6)
w.field('bsc_min','F',13,6)
w.field('post_avg','F',13,6)
w.field('post_min','F',13,6)
w.field('post_max','F',13,6)
w.field('risetime','F',13,6)
for shaperec in r.iterShapeRecords():
    rec = shaperec.record
    shp = shaperec.shape
    data_list = data_dict[rec.OBJECTID].copy()
    data_list.insert(2,'N/A' if np.isnan(data_dict[rec.OBJECTID][1]) else num2date(np.round(data_dict[rec.OBJECTID][1])+0.1).strftime('%Y/%m/%d'))
    data_list.insert(4,'N/A' if np.isnan(data_dict[rec.OBJECTID][2]) else num2date(np.round(data_dict[rec.OBJECTID][2])+0.1).strftime('%Y/%m/%d'))
    rec.extend(data_list)
    w.shape(shp)
    w.record(*rec)
w.close()
shutil.copy2(opts.inpnam+'.prj',opts.outnam+'.prj')
