#!/usr/bin/env python
import os
import sys
import re
from glob import glob
from datetime import datetime
try:
    import gdal
except Exception:
    from osgeo import gdal
try:
    import osr
except Exception:
    from osgeo import osr
import numpy as np
from matplotlib.dates import date2num,num2date
from csaps import UnivariateCubicSmoothingSpline
from statsmodels.stats.weightstats import DescrStatsW
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from optparse import OptionParser,IndentedHelpFormatter

# Default values
TMIN = '20190315'
TMAX = '20190615'
TMGN = 15.0
TSTP = 0.1
TSTR_1 = -10.0
TEND_1 = 10.0
TSTR_2 = 20.0
TEND_2 = 80.0
SMOOTH = 0.01
SEN1_DISTANCE = 10
SEN1_PROMINENCE = 0.1
XSGM = 4.0 # day
LSGM = 30.0 # m
N_NEAREST = 36
DATDIR = '.'
INCIDENCE_ANGLE = 'incidence_angle.dat'
NEAR_FNAM = 'find_nearest.npz'
AREA_FNAM = 'pixel_area_block.dat'
NPY_FNAM = 'output.npy'

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('-s','--tmin',default=TMIN,help='Min date of transplanting in the format YYYYMMDD (%default)')
parser.add_option('-e','--tmax',default=TMAX,help='Max date of transplanting in the format YYYYMMDD (%default)')
parser.add_option('--data_tmin',default=None,help='Min date of input data in the format YYYYMMDD (%default)')
parser.add_option('--data_tmax',default=None,help='Max date of input data in the format YYYYMMDD (%default)')
parser.add_option('--tmgn',default=TMGN,type='float',help='Margin of input data in day (%default)')
parser.add_option('--tstp',default=TSTP,type='float',help='Precision of transplanting date in day (%default)')
parser.add_option('--tstr_1',default=TSTR_1,type='float',help='Start day of transplanting period seen from the min. peak (%default)')
parser.add_option('--tend_1',default=TEND_1,type='float',help='End day of transplanting period seen from the min. peak (%default)')
parser.add_option('--tstr_2',default=TSTR_2,type='float',help='Start day of heading period seen from the min. peak (%default)')
parser.add_option('--tend_2',default=TEND_2,type='float',help='End day of heading period seen from the min. peak (%default)')
parser.add_option('-I','--incidence_angle',default=INCIDENCE_ANGLE,help='Incidence angle file, format: date(%Y%m%d) angle(deg) (%default)')
parser.add_option('-l','--incidence_list',default=None,help='Incidence angle list, format: flag(0|1=baseline) pol(VH|VV) angle(deg) filename (%default)')
parser.add_option('-S','--smooth',default=SMOOTH,type='float',help='Smoothing factor from 0 to 1 (%default)')
parser.add_option('--sen1_distance',default=SEN1_DISTANCE,type='int',help='Minimum peak distance in day for Sentinel-1 (%default)')
parser.add_option('--sen1_prominence',default=SEN1_PROMINENCE,type='float',help='Minimum prominence in dB for Sentinel-1 (%default)')
parser.add_option('-w','--xsgm',default=XSGM,type='float',help='Standard deviation of gaussian in day (%default)')
parser.add_option('-W','--lsgm',default=LSGM,type='float',help='Standard deviation of gaussian in m (%default)')
parser.add_option('--n_nearest',default=N_NEAREST,type='int',help='Number of nearest pixels to be considered (%default)')
parser.add_option('--output_epsg',default=None,type='int',help='Output EPSG (guessed from input data)')
parser.add_option('--near_fnam',default=NEAR_FNAM,help='Nearby index file name (%default)')
parser.add_option('--area_fnam',default=AREA_FNAM,help='Pixel area file name (%default)')
parser.add_option('--npy_fnam',default=NPY_FNAM,help='Output npy file name (%default)')
parser.add_option('-D','--datdir',default=DATDIR,help='Input data directory, not used if input_fnam is given (%default)')
parser.add_option('--search_key',default=None,help='Search key for input data, not used if input_fnam is given (%default)')
parser.add_option('-i','--inp_fnam',default=None,help='Input GeoTIFF name (%default)')
(opts,args) = parser.parse_args()

nmin = date2num(datetime.strptime(opts.tmin,'%Y%m%d'))
nmax = date2num(datetime.strptime(opts.tmax,'%Y%m%d'))
if opts.data_tmin is not None:
    dmin = datetime.strptime(opts.data_tmin,'%Y%m%d')
else:
    dmin = num2date(nmin+opts.tstr_1-opts.tmgn).replace(tzinfo=None)
if opts.data_tmax is not None:
    dmax = datetime.strptime(opts.data_tmax,'%Y%m%d')
else:
    dmax = num2date(nmax+opts.tend_2+opts.tmgn).replace(tzinfo=None)

# read nearby indices
data = np.load(opts.near_fnam)
sid_0 = data['sid_0']
for nn in range(1,opts.n_nearest+1):
    exec('sid_{} = data["sid_{}"]'.format(nn,nn))
    exec('leng_{} = data["leng_{}"]'.format(nn,nn))

object_ids = []
blocks = []
inds = []
areas = []
with open(opts.area_fnam,'r') as fp:
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
inds = np.array(inds)
areas = np.array(areas)
nobject = object_ids.size

if opts.inp_fnam is not None:
    ds = gdal.Open(opts.inp_fnam)
    if opts.output_epsg is None:
        prj = ds.GetProjection()
        srs = osr.SpatialReference(wkt=prj)
        epsg = srs.GetAttrValue('AUTHORITY',1)
        if re.search('\D',epsg):
            raise ValueError('Error in EPSG >>> '+epsg)
        output_epsg = int(epsg)
    else:
        output_epsg = opts.output_epsg
    data = ds.ReadAsArray()
    data_shape = data[0].shape
    data_trans = ds.GetGeoTransform()
    band_list = []
    for i in range(len(data)):
        band_list.append(ds.GetRasterBand(i+1).GetDescription())
    band_list = np.array(band_list)
    ds = None
else:
    output_epsg = None
    data_shape = None
    data_trans = None
    data = []
    band_list = []
    fs = sorted(glob(os.path.join(opts.datdir,'*'+'[0-9]'*8+'*.tif')))
    for fnam in fs:
        f = os.path.basename(fnam)
        if opts.search_key is not None and not re.search(opts.search_key,f):
            continue
        m = re.search('\D('+'\d'*8+')\D',f)
        if not m:
            m = re.search('^('+'\d'*8+')\D',f)
            if not m:
                raise ValueError('Error in finding date >>> '+f)
        dstr = m.group(1)
        d = datetime.strptime(dstr,'%Y%m%d')
        if d < dmin or d > dmax:
            continue
        sys.stderr.write(f+' '+dstr+'\n')
        ds = gdal.Open(fnam)
        if output_epsg is None:
            if opts.output_epsg is None:
                prj = ds.GetProjection()
                srs = osr.SpatialReference(wkt=prj)
                epsg = srs.GetAttrValue('AUTHORITY',1)
                if re.search('\D',epsg):
                    raise ValueError('Error in EPSG >>> '+epsg)
                output_epsg = int(epsg)
            else:
                output_epsg = opts.output_epsg
        dtmp = ds.ReadAsArray()
        if data_shape is None:
            data_shape = dtmp[0].shape
        elif dtmp[0].shape != data_shape:
            raise ValueError('Error, dtmp[0].shape={}, data_shape={}'.format(dtmp[0].shape,data_shape))
        trans = ds.GetGeoTransform()
        if data_trans is None:
            data_trans = trans
        elif trans != data_trans:
            raise ValueError('Error, trans={}, data_trans={}'.format(trans,data_trans))
        for i in range(len(dtmp)):
            data.append(dtmp[i])
            band_list.append(ds.GetRasterBand(i+1).GetDescription())
        ds = None
    data = np.array(data)
    band_list = np.array(band_list)

if opts.incidence_list is not None:
    incidence_flag = []
    incidence_angle = []
    incidence_fnam = []
    incidence_select = []
    with open (opts.incidence_list,'r') as fp:
        for line in fp:
            item = line.split()
            if len(item) < 4:
                continue
            if item[0][0] == '#':
                continue
            pol = item[1].upper()
            if pol != 'VH':
                continue
            incidence_flag.append(int(item[0]))
            incidence_angle.append(float(item[2]))
            incidence_fnam.append(item[3])
            if len(item) > 4:
                if item[4][0].lower() == 'f':
                    incidence_select.append(False)
                else:
                    incidence_select.append(True)
            else:
                incidence_select.append(True)
    incidence_flag = np.array(incidence_flag)
    incidence_angle = np.array(incidence_angle)
    incidence_select = np.array(incidence_select)
    nangle = incidence_angle.size
    cnd = (incidence_flag == 1)
    if cnd.sum() != 1:
        raise ValueError('Error in incidence_flag, cnd.sum()={}'.format(cnd.sum()))
    baseline_indx = np.arange(nangle)[cnd][0]
    signal_avg = []
    for fnam in incidence_fnam:
        avg = np.load(fnam)
        if avg.shape != data_shape:
            raise ValueError('Error, avg.shape={}, data_shape={}, fnam={}'.format(avg.shape,data_shape,fnam))
        signal_avg.append(avg)
    signal_dif = []
    for i in range(nangle):
        if i == baseline_indx:
            signal_dif.append(0.0)
        else:
            signal_dif.append(signal_avg[baseline_indx]-signal_avg[i])
    incidence_indx = {}
    with open(opts.incidence_angle,'r') as fp:
        for line in fp:
            item = line.split()
            if len(item) < 2:
                continue
            if item[0][0] == '#':
                continue
            ang = float(item[1])
            dif = np.abs(incidence_angle-ang)
            indx = np.argmin(dif)
            if dif[indx] > 0.1:
                raise ValueError('Error, dif={}'.format(dif[indx]))
            incidence_indx.update({item[0]:indx})

vh_dtim = []
vh_data = []
for i,band in enumerate(band_list):
    if not re.search('VH',band):
        continue
    m = re.search('_('+'\d'*8+')$',band)
    if not m:
        raise ValueError('Error in finding date >>> '+band)
    dstr = m.group(1)
    if opts.incidence_list is not None:
        if not incidence_select[incidence_indx[dstr]]:
            continue
        dtmp = (data[i]+signal_dif[incidence_indx[dstr]]).flatten()
    else:
        dtmp = data[i].flatten()
    data_avg = []
    for i in range(nobject):
        data_value = dtmp[inds[i]]
        data_weight = areas[i]
        cnd = ~np.isnan(data_value)
        if cnd.sum() <= 1:
            data_avg.append(data_value[cnd].mean())
        else:
            data_weighted_stats = DescrStatsW(data_value[cnd],weights=data_weight[cnd],ddof=0)
            data_avg.append(data_weighted_stats.mean)
    data_avg = np.array(data_avg)
    vh_dtim.append(datetime.strptime(dstr,'%Y%m%d'))
    vh_data.append(data_avg)
vh_dtim = np.array(vh_dtim)
vh_data = np.array(vh_data)
vh_ntim = date2num(vh_dtim)

k1_offset = int(opts.tstr_1/opts.tstp+(-0.1 if opts.tstr_1 < 0.0 else 0.1))
k2_offset = int(opts.tend_1/opts.tstp+(-0.1 if opts.tend_1 < 0.0 else 0.1))+1
k3_offset = int(opts.tstr_2/opts.tstp+(-0.1 if opts.tstr_2 < 0.0 else 0.1))
k4_offset = int(opts.tend_2/opts.tstp+(-0.1 if opts.tend_2 < 0.0 else 0.1))+1
xx = np.arange(np.floor(vh_ntim.min()),np.ceil(vh_ntim.max())+1.0,opts.tstp)
xpek_sid = [[] for i in range(nobject)]
ypek_sid = [[] for i in range(nobject)]
for i in range(nobject):
    yi = vh_data[:,i] # VH
    sp = UnivariateCubicSmoothingSpline(vh_ntim,yi,smooth=opts.smooth)
    yy = sp(xx)
    min_peaks,properties = find_peaks(-yy,distance=opts.sen1_distance,prominence=opts.sen1_prominence)
    if len(min_peaks) > 0:
        for k in min_peaks:
            if xx[k] < nmin or xx[k] > nmax:
                continue
            k1 = max(k+k1_offset,0)
            k2 = min(k+k2_offset,xx.size)
            vmin = yy[k1:k2].mean()
            k3 = min(k+k3_offset,xx.size-1)
            k4 = min(k+k4_offset,xx.size)
            vmax = yy[k3:k4].mean()
            if vmax > vmin:
                xpek_sid[i].append(xx[k])
                ypek_sid[i].append(vmax-vmin)
for i in range(nobject):
    xpek_sid[i] = np.array(xpek_sid[i])
    ypek_sid[i] = np.array(ypek_sid[i])

nb = 2
output_data = np.full((nb,nobject),np.nan)
for i in range(nobject):
    if i != sid_0[i]:
        raise ValueError('Error, i={}, sid_0={}'.format(i,sid_0[i]))
    indx = []
    lengs = []
    for nn in range(1,opts.n_nearest+1):
        exec('indx.append(sid_{}[i])'.format(nn))
        exec('lengs.append(leng_{}[i])'.format(nn))
    indx = np.array(indx)
    lengs = np.array(lengs)
    yy = np.zeros_like(xx)
    for xi,yi in zip(xpek_sid[i],ypek_sid[i]):
        ytmp = yi*np.exp(-0.5*np.square((xx-xi)/opts.xsgm))
        yy += ytmp
    for j,leng in zip(indx,lengs):
        fact = np.exp(-0.5*np.square(leng/opts.lsgm))
        for xj,yj in zip(xpek_sid[j],ypek_sid[j]):
            ytmp = yj*fact*np.exp(-0.5*np.square((xx-xj)/opts.xsgm))
            yy += ytmp
    k = np.argmax(yy)
    if xx[k] < nmin or xx[k] > nmax:
        continue
    output_data[0,i] = xx[k]
    output_data[1,i] = yy[k]
np.save(opts.npy_fnam,output_data)
