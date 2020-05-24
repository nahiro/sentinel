import re
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
SMOOTH = 0.01
SEN1_DISTANCE = 10
SEN1_PROMINENCE = 0.1
OUTNAM = 'output.tif'

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('-s','--tmin',default=TMIN,help='Min date in the format YYYYMMDD (%default)')
parser.add_option('-e','--tmax',default=TMAX,help='Max date in the format YYYYMMDD (%default)')
parser.add_option('-S','--smooth',default=SMOOTH,type='float',help='Smoothing factor from 0 to 1 (%default)')
parser.add_option('--sen1_distance',default=SEN1_DISTANCE,type='int',help='Minimum peak distance in day for Sentinel-1 (%default)')
parser.add_option('--sen1_prominence',default=SEN1_PROMINENCE,type='float',help='Minimum prominence in dB for Sentinel-1 (%default)')
parser.add_option('-o','--outnam',default=OUTNAM,help='Output GeoTIFF name (%default)')
(opts,args) = parser.parse_args()

nmin = date2num(datetime.strptime(opts.tmin,'%Y%m%d'))
nmax = date2num(datetime.strptime(opts.tmax,'%Y%m%d'))

# read nearby indices
data = np.load('find_nearest.npz')
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

fnam = 'collocate_all_resample.tif'
ds = gdal.Open(fnam)
data = ds.ReadAsArray()
trans = ds.GetGeoTransform()
band_list = []
for i in range(len(data)):
    band_list.append(ds.GetRasterBand(i+1).GetDescription())
band_list = np.array(band_list)
ds = None
nx = data.shape[2]
ny = data.shape[1]
ngrd = nx*ny
dtim = []
vh_indx = []
vv_indx = []
for i,band in enumerate(band_list):
    m = re.search('_('+'\d'*8+')$',band)
    if m:
        dtim.append(datetime.strptime(m.group(1),'%Y%m%d'))
    else:
        raise ValueError('Error in finding date >>> '+band)
    if re.search('VH',band):
        vh_indx.append(i)
    elif re.search('VV',band):
        vv_indx.append(i)
    else:
        raise ValueError('Error in band >>> '+band)
dtim = np.array(dtim)
ntim = date2num(dtim)
vh_indx = np.array(vh_indx)
vv_indx = np.array(vv_indx)

vh_data = data[vh_indx]
vh_band_list = band_list[vh_indx]
vh_dtim = dtim[vh_indx]
vh_ntim = ntim[vh_indx]

#vv_data = data[vv_indx]
#vv_band_list = band_list[vv_indx]
#vv_dtim = dtim[vv_indx]
#vv_ntim = ntim[vv_indx]

xx = np.arange(np.floor(vh_ntim.min()),np.ceil(vh_ntim.max())+1.0,0.1)
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
            sid = np.ravel_multi_index((i,j),(ny,nx))
            for k in min_peaks:
                if xx[k] < nmin or xx[k] > nmax:
                    continue
                k1 = max(k-100,0)
                k2 = min(k+101,xx.size)
                vmin = yy[k1:k2].mean()
                k1 = min(k+500,xx.size-1)
                k2 = min(k+701,xx.size)
                vmax = yy[k1:k2].mean()
                if vmax > vmin:
                    xpek_sid[sid].append(xx[k])
                    ypek_sid[sid].append(vmax-vmin)
for i in range(ngrd):
    xpek_sid[i] = np.array(xpek_sid[i])
    ypek_sid[i] = np.array(ypek_sid[i])

nb = 4
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
    output_data[0,indy,indx] = xx[k]
    output_data[1,indy,indx] = yy[k]
    #k = np.argmax(zz)
    #output_data[2,indy,indx] = xx[k]
    #output_data[3,indy,indx] = zz[k]
np.save('output_data.npy',output_data)

drv = gdal.GetDriverByName('GTiff')
ds = drv.Create(opts.outnam,nx,ny,nb,gdal.GDT_Float32)
ds.SetGeoTransform(trans)
srs = osr.SpatialReference()
srs.ImportFromEPSG(32748)
ds.SetProjection(srs.ExportToWkt())
band_name = ['xmin','ymin','vmin','xmax','ymax','vmax']
for i in range(nb):
    band = ds.GetRasterBand(i+1)
    band.WriteArray(output_data[i])
    band.SetDescription(band_name[i])
band.SetNoDataValue(np.nan) # The TIFFTAG_GDAL_NODATA only support one value per dataset
ds.FlushCache()
ds = None # close dataset
