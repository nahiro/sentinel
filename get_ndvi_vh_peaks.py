#!/usr/bin/env python
import os
import sys
import shutil
import re
from datetime import datetime,timedelta
import gdal
import osr
import shapefile
import numpy as np
from scipy.signal import find_peaks
from csaps import UnivariateCubicSmoothingSpline
from matplotlib.dates import date2num
from matplotlib.path import Path
from optparse import OptionParser,IndentedHelpFormatter

# Default values
SCL_MIN = 3.9
SCL_MAX = 7.1
SEN1_DISTANCE = 30
SEN2_DISTANCE = 30
SEN1_PROMINENCE = 0.5
SEN2_PROMINENCE = 0.1
SEN1_THRESHOLD = -5.0 # dB
SEN1_SEN2_DIF = 30.0 # day
POLARIZATION = 'VH'
MAXDIS = 100.0 # m
VINT = 100
SHPNAM = os.path.join('New_Test_Sites','New_Test_Sites')
OUTNAM = 'peak_data'
NPYNAM = 'peak_data.npy'
INCIDENCE_ANGLE = 'incidence_angle.dat'

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.set_usage('Usage: %prog sen1_collocated_geotiff_file sen2_collocated_geotiff_file [options]')
parser.add_option('-w','--wmin',default=None,type='float',help='Minimum growing time in day (%default)')
parser.add_option('-W','--wmax',default=None,type='float',help='Maximum growing time in day (%default)')
parser.add_option('-p','--pmin',default=None,help='Minimum planting date in the format YYYYMMDD (%default)')
parser.add_option('-P','--pmax',default=None,help='Maximum planting date in the format YYYYMMDD (%default)')
parser.add_option('-l','--scl_min',default=SCL_MIN,type='float',help='Minimum scene classification value (%default)')
parser.add_option('-L','--scl_max',default=SCL_MAX,type='float',help='Maximum scene classification value (%default)')
parser.add_option('--sen1_distance',default=SEN1_DISTANCE,type='int',help='Minimum peak distance in day for Sentinel-1 (%default)')
parser.add_option('--sen2_distance',default=SEN2_DISTANCE,type='int',help='Minimum peak distance in day for Sentinel-2 (%default)')
parser.add_option('--sen1_prominence',default=SEN1_PROMINENCE,type='float',help='Minimum prominence in dB for Sentinel-1 (%default)')
parser.add_option('--sen2_prominence',default=SEN2_PROMINENCE,type='float',help='Minimum prominence in reflectance for Sentinel-2 (%default)')
parser.add_option('--sen1_threshold',default=SEN1_THRESHOLD,type='float',help='Minimum signal in dB for Sentinel-1 (%default)')
parser.add_option('--sen1_sen2_dif',default=SEN1_SEN2_DIF,type='float',help='Maximum transplanting-date difference between Sentinel-1 and Sentinel-2 in day (%default)')
parser.add_option('--polarization',default=POLARIZATION,help='Polarization (%default)')
parser.add_option('-I','--incidence_angle',default=INCIDENCE_ANGLE,help='Incidence angle file, format: date(%Y%m%d) angle(deg) (%default)')
parser.add_option('-i','--incidence_list',default=None,help='Incidence angle list, format: flag(0|1=baseline) pol(VH|VV) angle(deg) filename (%default)')
parser.add_option('-m','--maxdis',default=MAXDIS,type='float',help='Max distance in m (%default)')
parser.add_option('-s','--shpnam',default=SHPNAM,help='Input shapefile name (%default)')
parser.add_option('-S','--sind',default=None,type='int',action='append',help='Selected indices (%default)')
parser.add_option('-o','--npynam',default=NPYNAM,help='Output NPY name (%default)')
parser.add_option('-O','--outnam',default=OUTNAM,help='Output shapefile name (%default)')
parser.add_option('--vint',default=VINT,type='int',help='Verbose output interval (%default)')
parser.add_option('-v','--verbose',default=False,action='store_true',help='Verbose mode (%default)')
(opts,args) = parser.parse_args()
if len(args) < 2:
    parser.print_help()
    sys.exit(0)
sen1_fnam = args[0]
sen2_fnam = args[1]
if opts.npynam.upper() == 'NONE':
    opts.npynam = None
if opts.outnam.upper() == 'NONE':
    opts.outnam = None

xstp = 10.0
ystp = -10.0
xmin,xmax,ymin,ymax = (743800.0,756800.0,9236000.0,9251800.0)
xg,yg = np.meshgrid(np.arange(xmin,xmax+0.1*xstp,xstp),np.arange(ymax,ymin-0.1*ystp,ystp))
ngrd = xg.size
ny,nx = xg.shape

if opts.incidence_list is not None:
    incidence_flag = []
    incidence_angle = []
    incidence_fnam = []
    with open (opts.incidence_list,'r') as fp:
        for line in fp:
            item = line.split()
            if len(item) < 4:
                continue
            if item[0][0] == '#':
                continue
            pol = item[1].upper()
            if pol != opts.polarization.upper():
                continue
            incidence_flag.append(int(item[0]))
            incidence_angle.append(float(item[2]))
            incidence_fnam.append(item[3])
    incidence_flag = np.array(incidence_flag)
    incidence_angle = np.array(incidence_angle)
    nangle = incidence_angle.size
    cnd = (incidence_flag == 1)
    if cnd.sum() != 1:
        raise ValueError('Error in incidence_flag, cnd.sum()={}'.format(cnd.sum()))
    baseline_indx = np.arange(nangle)[cnd][0]
    signal_avg = []
    for fnam in incidence_fnam:
        avg = np.load(fnam)
        if avg.shape != xg.shape:
            raise ValueError('Error, avg.shape={}, xg.shape={}, fnam={}'.format(avg.shape,xg.shape,fnam))
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

ibands = [4,8,17]
nband = len(ibands)

# Read Sentinel-1 data
ds = gdal.Open(sen1_fnam)
sen1_data = ds.ReadAsArray()
if sen1_data[0].shape != xg.shape:
    raise ValueError('Error, sen1_data[0].shape={}, xg.shape={}'.format(sen1_data[0].shape,xg.shape))
sen1_band_list = []
sen1_band_indx = []
sen1_dtim = []
for i in range(len(sen1_data)):
    band = ds.GetRasterBand(i+1).GetDescription()
    if not opts.polarization in band.upper():
        continue
    sen1_band_list.append(band)
    #sys.stderr.write(band+'\n')
    m = re.search('.*_(\d+)$',band)
    if not m:
        raise ValueError('Error in finding date >>> '+band)
    dstr = m.group(1)
    dtim = datetime.strptime(dstr,'%Y%m%d')
    sen1_band_indx.append(i)
    sen1_dtim.append(dtim)
    if opts.incidence_list is not None:
        sen1_data[i] += signal_dif[incidence_indx[dstr]]
ds = None
sen1_band_list = np.array(sen1_band_list)
sen1_band_indx = np.array(sen1_band_indx)
sen1_dtim = np.array(sen1_dtim)
sen1_ntim = date2num(sen1_dtim)

# Read Sentinel-2 data
ds = gdal.Open(sen2_fnam)
sen2_data = ds.ReadAsArray()
if sen2_data[0].shape != xg.shape:
    raise ValueError('Error, sen2_data[0].shape={}, xg.shape={}'.format(sen2_data[0].shape,xg.shape))
sen2_band_list = []
sen2_band_indx = [[] for i in ibands]
sen2_dtim_list = [[] for i in ibands]
for i in range(len(sen2_data)):
    band = ds.GetRasterBand(i+1).GetDescription()
    sen2_band_list.append(band)
    #sys.stderr.write(band+'\n')
    m = re.search('band_(\d+)_(\d+)$',band)
    if not m:
        raise ValueError('Error in finding date >>> '+band)
    band_num = int(m.group(1))
    dtim = datetime.strptime(m.group(2),'%Y%m%d')
    try:
        j = ibands.index(band_num)
        sen2_band_indx[j].append(i)
        sen2_dtim_list[j].append(dtim)
    except Exception:
        pass
ds = None
sen2_band_list = np.array(sen2_band_list)
sen2_band_indx = np.array(sen2_band_indx)
sen2_dtim_list = np.array(sen2_dtim_list)
for i in range(1,nband):
    if not np.all(sen2_dtim_list[i] == sen2_dtim_list[0]):
        raise ValueError('Error, different time {}'.format(i))
sen2_dtim = sen2_dtim_list[0].copy()
sen2_ntim = date2num(sen2_dtim)
if opts.pmin is not None:
    dmin = datetime.strptime(opts.pmin,'%Y%m%d')
else:
    dmin = sen2_dtim.min()
if opts.pmax is not None:
    dmax = datetime.strptime(opts.pmax,'%Y%m%d')
else:
    dmax = sen2_dtim.max()
pmin = date2num(dmin)
pmax = date2num(dmax)

all_bands = []
for i in range(nband):
    all_bands.append(sen2_data[sen2_band_indx[i]])
all_bands = np.array(all_bands)
all_bands[:nband-1] *= 1.0e-4

# Apply Flag and calculate NDVI -----------------------------------------------#
scl_data = all_bands[ibands.index(17)]
cnd = (scl_data < opts.scl_min) | (scl_data > opts.scl_max)
for i in range(0,nband-1):
    all_bands[i][cnd] = np.nan
b04_data = all_bands[ibands.index(4)]
b08_data = all_bands[ibands.index(8)]

ndvi = (b08_data-b04_data)/(b08_data+b04_data)

r = shapefile.Reader(opts.shpnam)
nr = len(r)
if opts.sind is not None:
    indi = opts.sind
else:
    indi = range(nr)
# Create output file
peak_data = np.full((8,nr),np.nan,dtype=np.float32)
xx = np.arange(np.floor(sen2_ntim.min()),np.ceil(sen2_ntim.max())+0.1,1.0)
for inum,ishp in enumerate(indi):
    if opts.verbose and inum%opts.vint == 0:
        sys.stderr.write('{} {}\n'.format(inum,len(indi)))
    shp = r.shape(ishp)
    p = Path(shp.points)
    flags = p.contains_points(np.hstack((xg.flatten()[:,np.newaxis],yg.flatten()[:,np.newaxis]))).reshape(xg.shape)
    ndat = flags.sum()
    sen1_yi = None
    sen2_yi = None
    if ndat > 0:
        dtmp = sen1_data[0][flags]
        cnd = ~np.isnan(dtmp)
        if cnd.sum() > 0:
            sen1_yi = sen1_data[:,flags].reshape(sen1_dtim.size,-1).mean(axis=1)
            sen2_yi = ndvi[:,flags].reshape(sen2_dtim.size,-1).mean(axis=1)
    else:
        pp = np.array(shp.points)
        xc = pp[:,0].mean()
        yc = pp[:,1].mean()
        if p.contains_point((xc,yc)):
            dp = np.square(xg-xc)+np.square(yg-yc)
            indx_y,indx_x = np.unravel_index(np.argmin(dp),xg.shape)
            dp_min = np.sqrt(dp[indx_y,indx_x])
            if dp_min < opts.maxdis:
                sen1_yi = sen1_data[:,indx_y,indx_x]
                sen2_yi = ndvi[:,indx_y,indx_x]
                ndat = -1
            else:
                sys.stderr.write('Case A, x={}, y={}, dp_min={}\n'.format(indx_x,indx_y,dp_min))
        else:
            dp_min = 1.0e10
            indx_y = None
            indx_x = None
            for ip in range(len(pp)):
                xc = pp[ip,0]
                yc = pp[ip,1]
                dp = np.square(xg-xc)+np.square(yg-yc)
                iy,ix = np.unravel_index(np.argmin(dp),xg.shape)
                dp_tmp = dp[iy,ix]
                if dp_tmp < dp_min:
                    indx_y = iy
                    indx_x = ix
                    dp_min = dp_tmp
            dp_min = np.sqrt(dp_min)
            if dp_min < opts.maxdis:
                sen1_yi = sen1_data[:,indx_y,indx_x]
                sen2_yi = ndvi[:,indx_y,indx_x]
                ndat = -2
            else:
                sys.stderr.write('Case B, x={}, y={}, dp_min={}\n'.format(indx_x,indx_y,dp_min))
    if sen1_yi is not None and sen2_yi is not None:
        # Smoothing
        cnd = ~np.isnan(sen1_yi)
        if cnd.sum() < 5:
            continue
        sen1_xc = sen1_ntim[cnd]
        sen1_yc = sen1_yi[cnd]
        sp = UnivariateCubicSmoothingSpline(sen1_xc,sen1_yc,smooth=1.0e-2)
        sen1_yy = sp(xx)

        # Peak search
        min_peaks,properties = find_peaks(-sen1_yy,distance=opts.sen1_distance,prominence=opts.sen1_prominence)
        if min_peaks.size < 1:
            if sen1_yy[0] < sen1_yy[-1]:
                min_peaks = np.append(0,min_peaks)
            else:
                continue
        sen1_x1 = xx[min_peaks]
        sen1_y1 = sen1_yy[min_peaks]

        # Smoothing
        cnd = ~np.isnan(sen2_yi)
        if cnd.sum() < 5:
            continue
        sen2_xc = sen2_ntim[cnd]
        sen2_yc = sen2_yi[cnd]
        sp = UnivariateCubicSmoothingSpline(sen2_xc,sen2_yc,smooth=2.0e-3)
        sen2_yy = sp(xx)

        # Peak search
        min_peaks,properties = find_peaks(-sen2_yy,distance=opts.sen2_distance,prominence=opts.sen2_prominence)
        max_peaks,properties = find_peaks(+sen2_yy,distance=opts.sen2_distance,prominence=opts.sen2_prominence)
        if min_peaks.size < 1 and max_peaks.size < 1:
            if sen2_yy[0] < sen2_yy[-1]:
                min_peaks = np.append(0,min_peaks)
                max_peaks = np.append(max_peaks,xx.size-1)
            else:
                continue
        elif min_peaks.size < 1:
            if max_peaks.size > 1:
                sys.stderr.write('Warning, more than one ({}) max peaks have been found, ishp={}.\n'.format(max_peaks.size,ishp))
            min_peaks = np.append(0,min_peaks)
        elif max_peaks.size < 1:
            if min_peaks.size > 1:
                sys.stderr.write('Warning, more than one ({}) min peaks have been found, ishp={}.\n'.format(min_peaks.size,ishp))
            max_peaks = np.append(max_peaks,xx.size-1)
        else:
            if max_peaks[0] < min_peaks[0]: # upslope
                min_peaks = np.append(0,min_peaks)
            if min_peaks[-1] > max_peaks[-1]: # upslope
                max_peaks = np.append(max_peaks,xx.size-1)
        if min_peaks.size != max_peaks.size:
            sys.stderr.write('Warning, min_peaks.size={}, max_peaks.size={}, ishp={}\n'.format(min_peaks.size,max_peaks.size,ishp))
            continue
        sen2_x1 = xx[min_peaks]
        sen2_x2 = xx[max_peaks]
        sen2_y1 = sen2_yy[min_peaks]
        sen2_y2 = sen2_yy[max_peaks]
        sen2_width = sen2_x2-sen2_x1
        sen2_height = sen2_y2-sen2_y1

        # Find planting/heading stage based on peak points
        cnd = (sen2_x1 >= pmin) & (sen2_x1 <= pmax) & (sen2_width >= opts.wmin) & (sen2_width <= opts.wmax)
        if cnd.sum() >= 1:
            indx = np.argmax(sen2_height[cnd])
            peak_data[0,ishp] = sen2_x1[cnd][indx]  # Put date of minNDVI (this is defined as planting stage)
            peak_data[1,ishp] = sen2_y1[cnd][indx]  # Put minNDVI value
            peak_data[2,ishp] = sen2_x2[cnd][indx]  # Put date of maxNDVI (this is defined as heading stage)
            peak_data[3,ishp] = sen2_y2[cnd][indx]  # Put maxNDVI value
            x1 = sen2_x1[cnd][indx]
        else:
            cnd = (sen2_x1 >= pmin) & (sen2_x1 <= pmax)
            if cnd.sum() >= 1:
                indx = np.argmax(sen2_height[cnd])
                peak_data[0,ishp] = sen2_x1[cnd][indx]  # Put date of minNDVI (this is defined as planting stage)
                peak_data[1,ishp] = sen2_y1[cnd][indx]  # Put minNDVI value
                peak_data[2,ishp] = sen2_x2[cnd][indx]  # Put date of maxNDVI (this is defined as heading stage)
                peak_data[3,ishp] = sen2_y2[cnd][indx]  # Put maxNDVI value
                x1 = sen2_x1[cnd][indx]
        cnd = (sen1_y1 < opts.sen1_threshold) & (np.abs(sen1_x1-x1)<opts.sen1_sen2_dif)
        if cnd.sum() >= 1:
            indx = np.argmin(sen1_y1[cnd])
            peak_data[4,ishp] = sen1_x1[cnd][indx]  # Put date of minVH (this is defined as planting stage)
            peak_data[5,ishp] = sen1_y1[cnd][indx]  # Put minVH value
        cnd = (sen1_y1 < opts.sen1_threshold)
        if cnd.sum() >= 1:
            indx = np.argmin(np.abs(sen1_x1[cnd]-x1))
            peak_data[6,ishp] = sen1_x1[cnd][indx]  # Put date of minVH (this is defined as planting stage)
            peak_data[7,ishp] = sen1_y1[cnd][indx]  # Put minVH value

if opts.npynam is not None:
    np.save(opts.npynam,peak_data)

# Output results
if opts.outnam is not None:
    w = shapefile.Writer(opts.outnam)
    w.shapeType = shapefile.POLYGON
    w.fields = r.fields[1:] # skip first deletion field
    band_name = ['planting_date','planting_ndvi','heading_date','heading_ndvi','sen1_min_date','sen1_min_peak','sen1_near_value','sen1_near_peak']
    for band in band_name:
        w.field(band,'F',13,6)
    for i,shaperec in enumerate(r.iterShapeRecords()):
        rec = shaperec.record
        shp = shaperec.shape
        for j in range(len(band_name)):
            rec.append(peak_data[j,i])
        w.shape(shp)
        w.record(*rec)
    w.close()
    shutil.copy2(opts.shpnam+'.prj',opts.outnam+'.prj')
