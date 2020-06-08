#!/usr/bin/env python
import re
from glob import glob
from datetime import datetime
import gdal
import osr
import numpy as np
from matplotlib.dates import date2num
from csaps import UnivariateCubicSmoothingSpline
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from optparse import OptionParser,IndentedHelpFormatter

# Default values
TMIN = '20190315'
TMAX = '20190615'
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
DATDIR = '.'
INCIDENCE_ANGLE = 'incidence_angle.dat'
NEAR_FNAM = 'find_nearest.npz'
OUT_FNAM = 'output.tif'

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('-s','--tmin',default=TMIN,help='Min date in the format YYYYMMDD (%default)')
parser.add_option('-e','--tmax',default=TMAX,help='Max date in the format YYYYMMDD (%default)')
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
parser.add_option('--output_epsg',default=None,type='int',help='Output EPSG (guessed from input data)')
parser.add_option('--near_fnam',default=NEAR_FNAM,help='Nearby index file name (%default)')
parser.add_option('--npy_fnam',default=None,help='Output npy file name (%default)')
parser.add_option('-D','--datdir',default=DATDIR,help='Input data directory (%default)')
parser.add_option('-i','--inp_fnam',default=None,help='Input GeoTIFF name (%default)')
parser.add_option('-o','--out_fnam',default=OUT_FNAM,help='Output GeoTIFF name (%default)')
(opts,args) = parser.parse_args()

nmin = date2num(datetime.strptime(opts.tmin,'%Y%m%d'))
nmax = date2num(datetime.strptime(opts.tmax,'%Y%m%d'))

# read nearby indices
data = np.load(opts.near_fnam)
sid_0 = data['sid_0']
sid_1 = data['sid_1']
sid_2 = data['sid_2']
sid_3 = data['sid_3']
sid_4 = data['sid_4']
sid_5 = data['sid_5']
sid_6 = data['sid_6']
sid_7 = data['sid_7']
sid_8 = data['sid_8']
sid_9 = data['sid_9']
sid_a = data['sid_a']
sid_b = data['sid_b']
sid_c = data['sid_c']
leng_1 = data['leng_1']
leng_2 = data['leng_2']
leng_3 = data['leng_3']
leng_4 = data['leng_4']
leng_5 = data['leng_5']
leng_6 = data['leng_6']
leng_7 = data['leng_7']
leng_8 = data['leng_8']
leng_9 = data['leng_9']
leng_a = data['leng_a']
leng_b = data['leng_b']
leng_c = data['leng_c']

if opts.inp_fnam is not None:
    ds = gdal.Open(opts.inp_fnam)
    prj = ds.GetProjection()
    srs = osr.SpatialReference(wkt=prj)
    if opts.output_epsg is None:
        epsg = srs.GetAttrValue('AUTHORITY',1)
        if re.search('\D',epsg):
            raise ValueError('Error in EPSG >>> '+epsg)
        output_epsg = int(epsg)
    else:
        output_epsg = opts.output_epsg
    data = ds.ReadAsArray()
    trans = ds.GetGeoTransform()
    band_list = []
    for i in range(len(data)):
        band_list.append(ds.GetRasterBand(i+1).GetDescription())
    band_list = np.array(band_list)
    ds = None
else:

data_shape = data[0].shape
nx = data_shape[1]
ny = data_shape[0]
ngrd = nx*ny

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
        vh_data.append(data[i]+signal_dif[incidence_indx[dstr]])
    else:
        vh_data.append(data[i])
    vh_dtim.append(datetime.strptime(dstr,'%Y%m%d'))
vh_dtim = np.array(vh_dtim)
vh_data = np.array(vh_data)
vh_ntim = date2num(vh_dtim)

k1_offset = int(opts.tstr_1/opts.tstp-0.1) # assuming negative value
k2_offset = int(opts.tend_1/opts.tstp+0.1)+1
k3_offset = int(opts.tstr_2/opts.tstp+0.1)
k4_offset = int(opts.tend_2/opts.tstp+0.1)+1
xx = np.arange(np.floor(vh_ntim.min()),np.ceil(vh_ntim.max())+1.0,opts.tstp)
xpek_sid = [[] for i in range(ngrd)]
ypek_sid = [[] for i in range(ngrd)]
for i in range(ny):
#for i in range(810,871):
    if i%100 == 0:
        print(i)
    for j in range(nx):
#    for j in range(810,841):
        yi = vh_data[:,i,j] # VH
        sp = UnivariateCubicSmoothingSpline(vh_ntim,yi,smooth=opts.smooth)
        yy = sp(xx)
        min_peaks,properties = find_peaks(-yy,distance=opts.sen1_distance,prominence=opts.sen1_prominence)
        if len(min_peaks) > 0:
            sid = np.ravel_multi_index((i,j),data_shape)
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
                    xpek_sid[sid].append(xx[k])
                    ypek_sid[sid].append(vmax-vmin)
for i in range(ngrd):
    xpek_sid[i] = np.array(xpek_sid[i])
    ypek_sid[i] = np.array(ypek_sid[i])

nb = 2
output_data = np.full((nb,ny,nx),np.nan)
for i in range(ngrd):
    if i != sid_0[i]:
        raise ValueError('Error, i={}, sid_0={}'.format(i,sid_0[i]))
    indx = np.array([sid_1[i],sid_2[i],sid_3[i],sid_4[i],sid_5[i],sid_6[i],sid_7[i],sid_8[i],sid_9[i],sid_a[i],sid_b[i],sid_c[i]])
    lengs = np.array([leng_1[i],leng_2[i],leng_3[i],leng_4[i],leng_5[i],leng_6[i],leng_7[i],leng_8[i],leng_9[i],leng_a[i],leng_b[i],leng_c[i]])
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
    indy,indx = np.unravel_index(i,data_shape)
    output_data[0,indy,indx] = xx[k]
    output_data[1,indy,indx] = yy[k]
if opts.npy_fnam is not None:
    np.save(opts.npy_fnam,output_data)

drv = gdal.GetDriverByName('GTiff')
ds = drv.Create(opts.out_fnam,nx,ny,nb,gdal.GDT_Float32)
ds.SetGeoTransform(trans)
srs = osr.SpatialReference()
srs.ImportFromEPSG(output_epsg)
ds.SetProjection(srs.ExportToWkt())
band_name = ['xpek','ypek']
for i in range(nb):
    band = ds.GetRasterBand(i+1)
    band.WriteArray(output_data[i])
    band.SetDescription(band_name[i])
band.SetNoDataValue(np.nan) # The TIFFTAG_GDAL_NODATA only support one value per dataset
ds.FlushCache()
ds = None # close dataset
