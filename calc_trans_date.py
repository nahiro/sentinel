#!/usr/bin/env python
import os
import sys
import re
from glob import glob
from datetime import datetime
import gdal
import osr
import json
from collections import OrderedDict
import numpy as np
from matplotlib.dates import date2num,num2date
from csaps import csaps
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from optparse import OptionParser,IndentedHelpFormatter

# Default values
TMIN = '20190315'
TMAX = '20190615'
TMGN = 30.0 # day
TSTP = 0.1 # day
TSTR = -20.0 # day
TEND = 20.0 # day
SMOOTH = 0.01
SEN1_DISTANCE = 10
SEN1_PROMINENCE = 0.1
SIG_DISTANCE = 10
SIG_PROMINENCE = 0.1
VTHR = -13.0 # dB
XSGM = 6.0 # day
LSGM = 30.0 # m
N_NEAREST = 120 # pixel
OFFSET = 0.0 # day
DATDIR = '.'
NEAR_FNAM = 'find_nearest.npz'
JSON_FNAM = 'output.json'
TEMP_FNAM = 'temporary.bin'
OUT_FNAM = 'output.tif'

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('-s','--tmin',default=TMIN,help='Min date of transplanting in the format YYYYMMDD (%default)')
parser.add_option('-e','--tmax',default=TMAX,help='Max date of transplanting in the format YYYYMMDD (%default)')
parser.add_option('--data_tmin',default=None,help='Min date of input data in the format YYYYMMDD (%default)')
parser.add_option('--data_tmax',default=None,help='Max date of input data in the format YYYYMMDD (%default)')
parser.add_option('--tmgn',default=TMGN,type='float',help='Margin of input data in day (%default)')
parser.add_option('--tstp',default=TSTP,type='float',help='Precision of transplanting date in day (%default)')
parser.add_option('--tstr',default=TSTR,type='float',help='Start day of transplanting period seen from the min. peak (%default)')
parser.add_option('--tend',default=TEND,type='float',help='End day of transplanting period seen from the min. peak (%default)')
parser.add_option('-x','--xmin',default=None,type='int',help='Min X index (inclusive, %default)')
parser.add_option('-X','--xmax',default=None,type='int',help='Max X index (exclusive, %default)')
parser.add_option('-y','--ymin',default=None,type='int',help='Min Y index (inclusive, %default)')
parser.add_option('-Y','--ymax',default=None,type='int',help='Max Y index (exclusive, %default)')
parser.add_option('-S','--smooth',default=SMOOTH,type='float',help='Smoothing factor from 0 to 1 (%default)')
parser.add_option('--sen1_distance',default=SEN1_DISTANCE,type='int',help='Minimum peak distance in tstp for Sentinel-1 (%default)')
parser.add_option('--sen1_prominence',default=SEN1_PROMINENCE,type='float',help='Minimum prominence in dB for Sentinel-1 (%default)')
parser.add_option('--sig_distance',default=SIG_DISTANCE,type='int',help='Minimum peak distance in tstp for synthesized signal (%default)')
parser.add_option('--sig_prominence',default=SIG_PROMINENCE,type='float',help='Minimum prominence in dB for synthesized signal (%default)')
parser.add_option('-v','--vthr',default=VTHR,type='float',help='Threshold (max value) of minimum VH in dB (%default)')
parser.add_option('-w','--xsgm',default=XSGM,type='float',help='Standard deviation of gaussian in day (%default)')
parser.add_option('-W','--lsgm',default=LSGM,type='float',help='Standard deviation of gaussian in m (%default)')
parser.add_option('--n_nearest',default=N_NEAREST,type='int',help='Number of nearest pixels to be considered (%default)')
parser.add_option('--offset',default=OFFSET,type='float',help='Offset of transplanting date in day (%default)')
parser.add_option('--output_epsg',default=None,type='int',help='Output EPSG (guessed from input data)')
parser.add_option('-I','--incidence_angle',default=None,help='Incidence angle file, format: date(%Y%m%d) angle(deg) (%default)')
parser.add_option('-l','--incidence_list',default=None,help='Incidence angle list, format: flag(0|1=baseline) pol(VH|VV) angle(deg) filename (%default)')
parser.add_option('-D','--datdir',default=DATDIR,help='Input data directory, not used if input_fnam is given (%default)')
parser.add_option('--search_key',default=None,help='Search key for input data, not used if input_fnam is given (%default)')
parser.add_option('--near_fnam',default=NEAR_FNAM,help='Nearby index file name (%default)')
parser.add_option('--mask_fnam',default=None,help='Mask file name (%default)')
parser.add_option('-i','--inp_fnam',default=None,help='Input GeoTIFF name (%default)')
parser.add_option('-j','--json_fnam',default=JSON_FNAM,help='Output JSON name (%default)')
parser.add_option('--temp_fnam',default=TEMP_FNAM,help='Temporary file name (%default)')
parser.add_option('-o','--out_fnam',default=OUT_FNAM,help='Output GeoTIFF name (%default)')
parser.add_option('--npy_fnam',default=None,help='Output npy file name (%default)')
parser.add_option('--early',default=False,action='store_true',help='Early estimation mode (%default)')
(opts,args) = parser.parse_args()

data_info = OrderedDict()
data_info['tmin'] = opts.tmin
data_info['tmax'] = opts.tmax
data_info['data_tmin'] = opts.data_tmin
data_info['data_tmax'] = opts.data_tmax
data_info['tmgn'] = opts.tmgn
data_info['tstp'] = opts.tstp
data_info['tstr'] = opts.tstr
data_info['tend'] = opts.tend
data_info['xmin'] = opts.xmin
data_info['xmax'] = opts.xmax
data_info['ymin'] = opts.ymin
data_info['ymax'] = opts.ymax
data_info['smooth'] = opts.smooth
data_info['sen1_distance'] = opts.sen1_distance
data_info['sen1_prominence'] = opts.sen1_prominence
data_info['sig_distance'] = opts.sig_distance
data_info['sig_prominence'] = opts.sig_prominence
data_info['vthr'] = opts.vthr
data_info['xsgm'] = opts.xsgm
data_info['lsgm'] = opts.lsgm
data_info['n_nearest'] = opts.n_nearest
data_info['offset'] = opts.offset
data_info['early'] = opts.early
data_info['incidence_angle'] = opts.incidence_angle if opts.incidence_angle is None else os.path.basename(opts.incidence_angle)
data_info['incidence_list'] = opts.incidence_list if opts.incidence_list is None else os.path.basename(opts.incidence_list)
data_info['near_fnam'] = opts.near_fnam if opts.near_fnam is None else os.path.basename(opts.near_fnam)
data_info['mask_fnam'] = opts.mask_fnam if opts.mask_fnam is None else os.path.basename(opts.mask_fnam)

nmin = date2num(datetime.strptime(opts.tmin,'%Y%m%d'))
nmax = date2num(datetime.strptime(opts.tmax,'%Y%m%d'))
if opts.data_tmin is not None:
    dmin = datetime.strptime(opts.data_tmin,'%Y%m%d')
else:
    dmin = num2date(nmin+opts.tstr-opts.tmgn).replace(tzinfo=None)
if opts.data_tmax is not None:
    dmax = datetime.strptime(opts.data_tmax,'%Y%m%d')
else:
    dmax = num2date(nmax+opts.tend+opts.tmgn).replace(tzinfo=None)
data_nmin = date2num(dmin)
data_nmax = date2num(dmax)

# read nearby indices
data = np.load(opts.near_fnam)
sid_0 = data['sid_0']
for nn in range(1,opts.n_nearest+1):
    exec('sid_{} = data["sid_{}"]'.format(nn,nn))
    exec('leng_{} = data["leng_{}"]'.format(nn,nn))

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
    iangle = {}
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
        meta = ds.GetMetadata()
        if 'iangle' in meta and 'PASS' in meta:
            ang = float(meta['iangle'])
            pas = meta['PASS'].upper()
            if pas == 'DESCENDING':
                ang *= -1.0
            elif pas != 'ASCENDING':
                raise ValueError('Error, PASS='+meta['PASS'])
            iangle.update({dstr:ang})
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
nx = data_shape[1]
ny = data_shape[0]
ngrd = nx*ny

if opts.mask_fnam is not None:
    ds = gdal.Open(opts.mask_fnam)
    mask = ds.ReadAsArray()
    ds = None
    if mask.shape != data_shape:
        raise ValueError('Error, mask.shape={}, data_shape={}'.format(mask.shape,data_shape))
    mask_flat = mask.flatten()

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
    if opts.incidence_angle is not None:
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
                    raise ValueError('Error, dif={} ()'.format(dif[indx],item[0]))
                incidence_indx.update({item[0]:indx})
    else:
        for dstr in iangle:
            ang = iangle[dstr]
            dif = np.abs(incidence_angle-ang)
            indx = np.argmin(dif)
            if dif[indx] > 0.1:
                raise ValueError('Error, dif={} ({})'.format(dif[indx],dstr))
            incidence_indx.update({dstr:indx})

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
        if not dstr in incidence_indx:
            raise ValueError('Error, no incidence angle for '+dstr)
        elif not incidence_select[incidence_indx[dstr]]:
            continue
        dtmp = data[i]+signal_dif[incidence_indx[dstr]]
    else:
        dtmp = data[i]
    #if opts.mask_fnam is not None:
    #    dtmp[mask < 0.5] = np.nan
    vh_data.append(dtmp)
    vh_dtim.append(datetime.strptime(dstr,'%Y%m%d'))
vh_dtim = np.array(vh_dtim)
vh_data = np.array(vh_data)
vh_ntim = date2num(vh_dtim)
data_info['dtim'] = ','.join([d.strftime('%Y%m%d') for d in vh_dtim])
with open(opts.json_fnam,'w') as json_file:
    json.dump(data_info,json_file,indent=4)
    json_file.write('\n') # Add newline cause PyJSON does not

k1_offset = int(opts.tstr/opts.tstp+(-0.1 if opts.tstr < 0.0 else 0.1))
k2_offset = int(opts.tend/opts.tstp+(-0.1 if opts.tend < 0.0 else 0.1))+1
xx = np.arange(data_nmin,data_nmax+0.1*opts.tstp,opts.tstp)
cnd = (xx > data_nmin-1.0e-6) & (xx < data_nmax+1.0e-6)
xc_indx = np.where(cnd)[0]
if xc_indx.size < 3:
    raise ValueError('Error, not enough data, xc_indx.size={}'.format(xc_indx.size))
xc_ind0 = xc_indx[0]
xc_ind1 = xc_indx[1]
xc_ind_1 = xc_indx[-1]
xc_ind_2 = xc_indx[-2]
xpek_sid = [[] for indp in range(ngrd)]
ypek_sid = [[] for indp in range(ngrd)]
min_peaks_list = [None for indp in range(ngrd)]

if opts.xmin is None:
    opts.xmin = 0
if opts.xmax is None:
    opts.xmax = nx
if opts.ymin is None:
    opts.ymin = 0
if opts.ymax is None:
    opts.ymax = ny
with open(opts.temp_fnam,'wb') as fp:
    for indy in range(opts.ymin,opts.ymax):
        if indy%100 == 0:
            sys.stderr.write('{}\n'.format(indy))
        for indx in range(opts.xmin,opts.xmax):
            indp = np.ravel_multi_index((indy,indx),data_shape)
            yi = vh_data[:,indy,indx] # VH
            yy = csaps(vh_ntim,yi,xx,smooth=opts.smooth)
            yy.tofile(fp)
            min_peaks,properties = find_peaks(-yy,distance=opts.sen1_distance,prominence=opts.sen1_prominence)
            if opts.early:
                if not xc_ind_1 in min_peaks:
                    if (yy[xc_ind_1] < opts.vthr) & (yy[xc_ind_1] < yy[xc_ind_2]):
                        min_peaks = np.append(min_peaks,xc_ind_1)
            min_peaks_list[indp] = min_peaks.astype(np.int16)
            if len(min_peaks) > 0:
                for k in min_peaks:
                    if k == xc_ind_1:
                        vmin = yy[k]
                        if vmin < opts.vthr:
                            xpek_sid[indp].append(xx[k])
                            ypek_sid[indp].append(opts.vthr-vmin)
                    else:
                        k1 = max(k+k1_offset,0)
                        k2 = min(k+k2_offset,xx.size)
                        vmin = yy[k1:k2].mean()
                        if vmin < opts.vthr:
                            xpek_sid[indp].append(xx[k]+opts.offset)
                            ypek_sid[indp].append(opts.vthr-vmin)
for indp in range(ngrd):
    xpek_sid[indp] = np.array(xpek_sid[indp])
    ypek_sid[indp] = np.array(ypek_sid[indp])

nb = 8
output_data = np.full((nb,ny,nx),np.nan)
with open(opts.temp_fnam,'rb') as fp:
    for indy in range(opts.ymin,opts.ymax):
        for indx in range(opts.xmin,opts.xmax):
            indp = np.ravel_multi_index((indy,indx),data_shape)
            if indp != sid_0[indp]:
                raise ValueError('Error, indp={}, sid_0={}'.format(indp,sid_0[indp]))
            inds = []
            lengs = []
            for nn in range(1,opts.n_nearest+1):
                exec('inds.append(sid_{}[indp])'.format(nn))
                exec('lengs.append(leng_{}[indp])'.format(nn))
            inds = np.array(inds)
            lengs = np.array(lengs)
            ys = np.zeros_like(xx)
            ws = 1.0
            for xj,yj in zip(xpek_sid[indp],ypek_sid[indp]):
                ytmp = yj*np.exp(-0.5*np.square((xx-xj)/opts.xsgm))
                ys += ytmp
            if opts.mask_fnam is not None:
                for indq,leng in zip(inds,lengs):
                    if mask_flat[indq] < 0.5:
                        continue
                    fact = np.exp(-0.5*np.square(leng/opts.lsgm))
                    ws += fact
                    for xj,yj in zip(xpek_sid[indq],ypek_sid[indq]):
                        ytmp = fact*yj*np.exp(-0.5*np.square((xx-xj)/opts.xsgm))
                        ys += ytmp
            else:
                for indq,leng in zip(inds,lengs):
                    fact = np.exp(-0.5*np.square(leng/opts.lsgm))
                    ws += fact
                    for xj,yj in zip(xpek_sid[indq],ypek_sid[indq]):
                        ytmp = fact*yj*np.exp(-0.5*np.square((xx-xj)/opts.xsgm))
                        ys += ytmp
            ys *= 1.0/ws
            max_peaks,properties = find_peaks(ys,distance=opts.sig_distance,prominence=opts.sig_prominence)
            if opts.early:
                max_peaks = np.append(max_peaks,xc_ind_1)
            if len(max_peaks) < 1:
                continue
            xp = xx[max_peaks]
            yp = ys[max_peaks]
            yp[(xp < nmin-1.0e-6) | (xp > nmax+1.0e-6)] = -1.0e10
            k = np.argmax(yp)
            if yp[k] < -1.0e2:
                continue
            output_data[0,indy,indx] = xp[k]
            output_data[1,indy,indx] = yp[k]
            yy = np.fromfile(fp,count=xx.size)
            min_peaks = min_peaks_list[indp]
if os.path.exists(opts.temp_fnam):
    os.remove(opts.temp_fnam)
if opts.npy_fnam is not None:
    np.save(opts.npy_fnam,output_data)

drv = gdal.GetDriverByName('GTiff')
ds = drv.Create(opts.out_fnam,nx,ny,nb,gdal.GDT_Float32)
ds.SetGeoTransform(data_trans)
srs = osr.SpatialReference()
srs.ImportFromEPSG(output_epsg)
ds.SetProjection(srs.ExportToWkt())
ds.SetMetadata({'offset':'{:.4f}'.format(opts.offset)})
band_name = ['trans_d','trans_s','trans_n','bsc_min','post_avg','post_min','post_max','risetime']
for i in range(nb):
    band = ds.GetRasterBand(i+1)
    band.WriteArray(output_data[i])
    band.SetDescription(band_name[i])
band.SetNoDataValue(np.nan) # The TIFFTAG_GDAL_NODATA only support one value per dataset
ds.FlushCache()
ds = None # close dataset
