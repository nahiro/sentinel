#!/usr/bin/env python
import os
import sys
import re
import warnings
import numpy as np
import gdal
from statsmodels.stats.weightstats import DescrStatsW
from optparse import OptionParser,IndentedHelpFormatter

# Default values
BAND = 'NDVI'
BAND_COL = 1
AREA_FNAM = 'pixel_area_block.dat'

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.set_usage('Usage: %prog input_fnam [options]')
parser.add_option('--band',default=BAND,help='Target band (%default)')
parser.add_option('-B','--band_fnam',default=None,help='Band file name (%default)')
parser.add_option('--band_col',default=BAND_COL,help='Band column number (%default)')
parser.add_option('--area_fnam',default=AREA_FNAM,help='Pixel area file name (%default)')
parser.add_option('--factor_fnam',default=None,help='Factor file name (%default)')
parser.add_option('--offset_fnam',default=None,help='Offset file name (%default)')
parser.add_option('--debug',default=False,action='store_true',help='Debug mode (%default)')
(opts,args) = parser.parse_args()
if len(args) < 1:
    parser.print_help()
    sys.exit(0)
input_fnam = args[0]
m = re.search('^('+'\d'*8+')_',os.path.basename(input_fnam))
if not m:
    raise ValueError('Error in file name >>> '+input_fnam)
dstr = m.group(1)
if not opts.debug:
    warnings.simplefilter('ignore')
if opts.band.upper() == 'NDVI':
    band_s = 'ndvi'
else:
    band_s = 'band'+opts.band
if opts.factor_fnam is None:
    opts.factor_fnam = 'atcor_{}_factor_{}.npy'.format(band_s,dstr)
if not os.path.exists(opts.factor_fnam):
    raise IOError('Error, no such file >>> '+opts.factor_fnam)
else:
    factor = np.load(opts.factor_fnam)
offset = np.load(os.path.join(datdir,'atcor_{}_offset_{}.npy'.format(band_s,dstr)))

ds = gdal.Open(input_fnam)
data = ds.ReadAsArray()

band_list = []
if opts.band_fnam is not None:
    with open(opts.band_fnam,'r') as fp:
        for line in fp:
            item = line.split()
            if len(item) <= opts.band_col or item[0][0]=='#':
                continue
            band_list.append(item[opts.band_col])
    if len(data) != len(band_list):
        raise ValueError('Error, len(data)={}, len(band_list)={} >>> {}'.format(len(data),len(band_list),input_fnam))
else:
    for i in range(ds.RasterCount):
        band = ds.GetRasterBand(i+1)
        band_list.append(band.GetDescription())
ds = None
if opts.band.upper() == 'NDVI':
    band4_index = band_list.index('B4')
    band8_index = band_list.index('B8')
else:
    band_index = band_list.index('B{}'.format(opts.band))
if opts.band.upper() == 'NDVI':
    b4_img = data[band4_index].flatten()
    b8_img = data[band8_index].flatten()
    data_img = (b8_img-b4_img)/(b8_img+b4_img)
else:
    data_img = data[band_index].flatten()*1.0e-4

object_ids = []
blocks = []
inds = []
areas = []
with open(opts.area_fnam,'r') as fp:
    for line in fp:
        item = line.split()
        if len(item) < 3 or item[0] == '#':
            continue
        object_ids.append(int(item[0]))
        blocks.append(item[1])
        inds.append([])
        areas.append([])
        n = int(item[2])
        if len(item) < 5 and n != 0:
            raise ValueError('Error, len(item)={}, n={}, expected n=0.'.format(len(item),n))
        for nn in range(3,n*2+3,2):
            inds[-1].append(int(item[nn]))
            areas[-1].append(float(item[nn+1]))
        inds[-1] = np.array(inds[-1])
        areas[-1] = np.array(areas[-1])
object_ids = np.array(object_ids)
blocks = np.array(blocks)
nobject = object_ids.size

data_org = []
for iobj in range(nobject):
    if inds[iobj].size < 1:
        data_org.append(np.nan)
        continue
    data_value = data_img[inds[iobj]]
    data_weight = areas[iobj]
    cnd = ~np.isnan(data_value)
    if cnd.sum() <= 1:
        data_org.append(data_value[cnd].mean())
    else:
        data_weighted_stats = DescrStatsW(data_value[cnd],weights=data_weight[cnd],ddof=0)
        data_org.append(data_weighted_stats.mean)
data_org = np.array(data_org)
data_cor = data_org*factor+offset
