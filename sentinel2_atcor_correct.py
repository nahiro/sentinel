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
BAND = '4'
AREA_FNAM = 'pixel_area_block.dat'

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('--band',default=BAND,help='Target band (%default)')
parser.add_option('--area_fnam',default=AREA_FNAM,help='Pixel area file name (%default)')
parser.add_option('-B','--band_fnam',default=None,help='Band file name (%default)')
parser.add_option('--debug',default=False,action='store_true',help='Debug mode (%default)')
(opts,args) = parser.parse_args()
if not opts.debug:
    warnings.simplefilter('ignore')

if opts.band.upper() == 'NDVI':
    band_s = 'nv'
else:
    band_s = 'b'+opts.band

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

band_list = []
with open(opts.band_fnam,'r') as fp:
    for line in fp:
        item = line.split()
        if len(item) != 2:
            raise ValueError('Error, len(item)={}'.format(len(item)))
        band_list.append(item[1])
if opts.band.upper() == 'NDVI':
    band4_index = band_list.index('B4')
    band8_index = band_list.index('B8')
else:
    band_index = band_list.index('B{}'.format(opts.band))

datdir = 'atcor_'+band_s
data_org = []
data_cor = []
for f in sorted(os.listdir(datdir)):
    m = re.search('atcor_{}_factor_(\d+)\.npy'.format(band_s),f)
    if not m:
        continue
    dstr = m.group(1)
    sys.stderr.write(dstr+'\n')

    fnam = os.path.join('/home/naohiro/Work/Sentinel-2/L2A/Bojongsoang/resample',dstr+'_geocor_resample.tif')
    ds = gdal.Open(fnam)
    data = ds.ReadAsArray()
    if len(data) != len(band_list):
        raise ValueError('Error, len(data)={}, len(band_list)={} >>> {}'.format(len(data),len(band_list),fnam))
    ds = None
    if opts.band.upper() == 'NDVI':
        b4_img = data[band4_index].flatten()
        b8_img = data[band8_index].flatten()
        data_img = (b8_img-b4_img)/(b8_img+b4_img)
    else:
        data_img = data[band_index].flatten()*1.0e-4

    factor = np.load(os.path.join(datdir,'atcor_{}_factor_{}.npy'.format(band_s,dstr)))
    offset = np.load(os.path.join(datdir,'atcor_{}_offset_{}.npy'.format(band_s,dstr)))

    data_avg = []
    for iobj in range(nobject):
        if inds[iobj].size < 1:
            data_avg.append(np.nan)
            continue
        data_value = data_img[inds[iobj]]
        data_weight = areas[iobj]
        cnd = ~np.isnan(data_value)
        if cnd.sum() <= 1:
            data_avg.append(data_value[cnd].mean())
        else:
            data_weighted_stats = DescrStatsW(data_value[cnd],weights=data_weight[cnd],ddof=0)
            data_avg.append(data_weighted_stats.mean)
    data_avg = np.array(data_avg)
    data_org.append(data_avg)
    data_cor.append(data_avg*factor+offset)
data_org = np.array(data_org)
data_cor = np.array(data_cor)
np.save('atcor_{}_original.npy'.format(band_s),data_org)
np.save('atcor_{}_corrected.npy'.format(band_s),data_cor)
