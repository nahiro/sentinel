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
INPNAM = os.path.join(HOME,'Work','SATREPS','Shapefile','field_GIS','All_area_polygon_20210914','All_area_polygon_20210914')
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
parser.add_option('-n','--no_block',default=False,action='store_true',help='No block in area_file (%default)')
parser.add_option('--use_index',default=False,action='store_true',help='Use index instead of OBJECTID (%default)')
parser.add_option('-d','--debug',default=False,action='store_true',help='Debug mode (%default)')
(opts,args) = parser.parse_args()

ds = gdal.Open(opts.data_file)
data = ds.ReadAsArray()
ds = None
data = data.reshape(len(data),-1).astype(np.float64)

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
nobject = object_ids.size

scheduled = {}
if opts.true_file is not None:
    with open(opts.true_file,'r') as fp:
        for line in fp:
            item = line.split()
            scheduled.update({item[0]:date2num(datetime.strptime(item[1],'%Y/%m/%d'))})

data_dict = {}
for i in range(object_ids.size):
    if opts.use_index:
        object_id = i+1
    else:
        object_id = object_ids[i]
    data_list = [None]*12
    # block, tanam_d
    if not opts.no_block:
        block = blocks[i]
        if block in scheduled:
            tanam_d = scheduled[block]
        else:
            tanam_d = np.nan
    else:
        block = 'None'
        tanam_d = np.nan
    data_list[0] = block
    data_list[1] = tanam_d
    # tanam_t
    tanam_t = 'N/A' if np.isnan(tanam_d) else num2date(np.round(tanam_d)+0.1).strftime('%Y/%m/%d')
    data_list[2] = tanam_t
    # trans_d
    data_value = data[0,inds[i]]
    #data_weight = areas[i]
    data_weight = data[1,inds[i]]
    cnd = ~np.isnan(data_value)
    dcnd = data_value[cnd].copy()
    wcnd = data_weight[cnd].copy()
    wsum = wcnd.sum()
    trans_d = (dcnd*wcnd).sum()/wsum
    data_list[3] = trans_d
    # trans_t
    trans_t = 'N/A' if np.isnan(trans_d) else num2date(np.round(trans_d)+0.1).strftime('%Y/%m/%d')
    data_list[4] = trans_t
    # trans_s
    trans_s = wsum/wcnd.size
    data_list[5] = trans_s
    # trans_n
    data_value = data[2,inds[i]]
    cnd = ~np.isnan(data_value)
    dcnd = data_value[cnd].copy()
    wcnd = data_weight[cnd].copy()
    wsum = wcnd.sum()
    trans_n = (dcnd*wcnd).sum()/wsum
    data_list[6] = trans_n
    # bsc_min
    data_value = data[3,inds[i]]
    dcnd = data_value[cnd].copy()
    bsc_min = (dcnd*wcnd).sum()/wsum
    data_list[7] = bsc_min
    # post_avg
    data_value = data[4,inds[i]]
    cnd = ~np.isnan(data_value)
    dcnd = data_value[cnd].copy()
    wcnd = data_weight[cnd].copy()
    wsum = wcnd.sum()
    post_avg = (dcnd*wcnd).sum()/wsum
    data_list[8] = post_avg
    # post_min
    data_value = data[5,inds[i]]
    dcnd = data_value[cnd].copy()
    post_min = (dcnd*wcnd).sum()/wsum
    data_list[9] = post_min
    # post_max
    data_value = data[6,inds[i]]
    dcnd = data_value[cnd].copy()
    post_max = (dcnd*wcnd).sum()/wsum
    data_list[10] = post_max
    # risetime
    data_value = data[7,inds[i]]
    cnd = ~np.isnan(data_value)
    dcnd = data_value[cnd].copy()
    wcnd = data_weight[cnd].copy()
    wsum = wcnd.sum()
    risetime = (dcnd*wcnd).sum()/wsum
    data_list[11] = risetime
    data_dict.update({object_id:data_list})

r = shapefile.Reader(opts.inpnam)
if len(r) != nobject:
    raise ValueError('Error, len(r)={}, nobject={}'.format(len(r),nobject))
w = shapefile.Writer(opts.outnam)
w.shapeType = shapefile.POLYGON
w.fields = r.fields[1:] # skip first deletion field
w.field('block','C',6,0)
w.field('tanam_d','F',13,6)
w.field('tanam_t','C',10,0)
w.field('trans_d','F',13,6)
w.field('trans_t','C',10,0)
w.field('trans_s','F',13,6)
w.field('trans_n','F',13,6)
w.field('bsc_min','F',13,6)
w.field('post_avg','F',13,6)
w.field('post_min','F',13,6)
w.field('post_max','F',13,6)
w.field('risetime','F',13,6)
for iobj,shaperec in enumerate(r.iterShapeRecords()):
    rec = shaperec.record
    shp = shaperec.shape
    if opts.use_index:
        object_id = iobj+1
    else:
        object_id = getattr(rec,'OBJECTID')
    data_list = data_dict[object_id]
    rec.extend(data_list)
    w.shape(shp)
    w.record(*rec)
w.close()
shutil.copy2(opts.inpnam+'.prj',opts.outnam+'.prj')
