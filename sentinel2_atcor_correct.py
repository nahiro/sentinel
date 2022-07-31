#!/usr/bin/env python
import os
import sys
import re
import warnings
import numpy as np
try:
    import gdal
except Exception:
    from osgeo import gdal
from statsmodels.stats.weightstats import DescrStatsW
from optparse import OptionParser,IndentedHelpFormatter

# Default values
BAND = 'NDVI'
BAND_COL = 1
BAND4_MAX = 0.35
R_MIN = 0.3
AREA_FNAM = 'pixel_area_block.dat'

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.set_usage('Usage: %prog input_fnam [options]')
parser.add_option('--band',default=BAND,help='Target band (%default)')
parser.add_option('-B','--band_fnam',default=None,help='Band file name (%default)')
parser.add_option('--band_col',default=BAND_COL,help='Band column number (%default)')
parser.add_option('--band4_max',default=BAND4_MAX,type='float',help='Band4 threshold (%default)')
parser.add_option('--r_min',default=R_MIN,type='float',help='R threshold (%default)')
parser.add_option('--area_fnam',default=AREA_FNAM,help='Pixel area file name (%default)')
parser.add_option('--param_fnam',default=None,help='Atcor parameter file name (%default)')
parser.add_option('-o','--output_fnam',default=None,help='Output NPZ name (%default)')
parser.add_option('--ignore_band4',default=False,action='store_true',help='Ignore exceeding the band4 threshold (%default)')
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
    band_l = 'ndvi'
else:
    band_l = 'band'+opts.band
if opts.param_fnam is None:
    opts.param_fnam = 'atcor_param_{}_{}.npz'.format(band_l,dstr)
if not os.path.exists(opts.param_fnam):
    raise IOError('Error, no such file >>> '+opts.param_fnam)
if opts.output_fnam is None:
    opts.output_fnam = 'atcor_data_{}_{}.npz'.format(band_l,dstr)
param = np.load(opts.param_fnam)
corcoef = param['corcoef']
factor = param['factor']
offset = param['offset']
cnd = (corcoef < opts.r_min)
factor[cnd] = np.nan
offset[cnd] = np.nan
npar = factor.size

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
        band_name = band.GetDescription()
        if not band_name:
            raise ValueError('Error, faild to read band name >>> {}'.format(input_fnam))
        band_list.append(band_name)
ds = None
if opts.band.upper() == 'NDVI':
    band_name = 'B4'
    if not band_name in band_list:
        raise ValueError('Error, faild to search index for {}'.format(band_name))
    band4_index = band_list.index(band_name)
    b4_img = data[band4_index].astype(np.float64).flatten()
    band_name = 'B8'
    if not band_name in band_list:
        raise ValueError('Error, faild to search index for {}'.format(band_name))
    band8_index = band_list.index(band_name)
    b8_img = data[band8_index].astype(np.float64).flatten()
    data_img = (b8_img-b4_img)/(b8_img+b4_img)
    if not opts.ignore_band4:
        b4_img *= 1.0e-4
else:
    band_name = 'B{}'.format(opts.band)
    if not band_name in band_list:
        raise ValueError('Error, faild to search index for {}'.format(band_name))
    band_index = band_list.index(band_name)
    data_img = data[band_index].astype(np.float64).flatten()*1.0e-4
    if not opts.ignore_band4:
        if opts.band == 4:
            b4_img = data_img.copy()
        else:
            band_name = 'B4'
            if not band_name in band_list:
                raise ValueError('Error, faild to search index for {}'.format(band_name))
            band4_index = band_list.index(band_name)
            b4_img = data[band4_index].astype(np.float64).flatten()*1.0e-4

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
if nobject != npar:
    raise ValueError('Error, nobject={}, npar={}'.format(nobject,npar))

data_org = []
for iobj in range(nobject):
    if inds[iobj].size < 1:
        data_org.append(np.nan)
        continue
    data_value = data_img[inds[iobj]]
    data_weight = areas[iobj]
    if opts.ignore_band4:
        cnd = ~np.isnan(data_value)
    else:
        data_band4 = b4_img[inds[iobj]]
        cnd = (~np.isnan(data_value)) & (~np.isnan(data_band4)) & (data_band4 < opts.band4_max)
    if cnd.sum() <= 1:
        data_org.append(data_value[cnd].mean())
    else:
        data_weighted_stats = DescrStatsW(data_value[cnd],weights=data_weight[cnd],ddof=0)
        data_org.append(data_weighted_stats.mean)
data_org = np.array(data_org)
data_cor = data_org*factor+offset
np.savez(opts.output_fnam,data_org=data_org,data_cor=data_cor)
