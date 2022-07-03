#!/usr/bin/env python
#from shapely.geometry.polygon import Polygon
import os
import sys
import warnings
from datetime import datetime
import numpy as np
from scipy.interpolate import griddata
import shapefile
import shapely
try:
    import gdal
except Exception:
    from osgeo import gdal
try:
    import osr
except Exception:
    from osgeo import osr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
from optparse import OptionParser,IndentedHelpFormatter

# Default values
HOME = os.environ.get('HOME')
if HOME is None:
    HOME = os.environ.get('HOMEPATH')
TMIN = '20190415'
TMAX = '20190601'
SMIN = 0.0
SMAX = 6.0
COORDS_COLOR = '#aaaaaa'
BLOCK_FNAM = os.path.join(HOME,'Work','SATREPS','Shapefile','studyarea','studyarea.shp')
TRANS_FNAM = os.path.join('.','output.tif')
MASK_FNAM = os.path.join('.','paddy_mask4_resample.tif')
OUTPUT_FNAM = 'trans_date_cihea.pdf'

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('-s','--tmin',default=TMIN,help='Min date in the format YYYYMMDD (%default)')
parser.add_option('-e','--tmax',default=TMAX,help='Max date in the format YYYYMMDD (%default)')
parser.add_option('--smin',default=SMIN,type='float',help='Min trans_s in dB (%default)')
parser.add_option('--smax',default=SMAX,type='float',help='Max trans_s in dB (%default)')
parser.add_option('-t','--title',default=None,help='Figure title (%default)')
parser.add_option('--block_fnam',default=BLOCK_FNAM,help='Block shape file (%default)')
parser.add_option('--trans_fnam',default=TRANS_FNAM,help='Transplanting tiff file (%default)')
parser.add_option('--mask_fnam',default=MASK_FNAM,help='Mask file (%default)')
parser.add_option('--output_fnam',default=OUTPUT_FNAM,help='Output figure name (%default)')
parser.add_option('--add_tmin',default=False,action='store_true',help='Add tmin in colorbar (%default)')
parser.add_option('--add_tmax',default=False,action='store_true',help='Add tmax in colorbar (%default)')
parser.add_option('--add_coords',default=False,action='store_true',help='Add geographical coordinates (%default)')
parser.add_option('--coords_color',default=COORDS_COLOR,help='Color of geographical coordinates (%default)')
parser.add_option('--early',default=False,action='store_true',help='Early estimation mode (%default)')
parser.add_option('-b','--batch',default=False,action='store_true',help='Batch mode (%default)')
parser.add_option('--debug',default=False,action='store_true',help='Debug mode (%default)')
(opts,args) = parser.parse_args()
if not opts.debug:
    warnings.simplefilter('ignore')

if opts.batch:
    import matplotlib
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap,LinearSegmentedColormap,to_rgba
from matplotlib.dates import date2num,num2date
from matplotlib.path import Path

def transform_wgs84_to_utm(longitude,latitude):
    utm_zone = (int(1+(longitude.mean()+180.0)/6.0))
    is_northern = (1 if latitude.mean() > 0 else 0)
    utm_coordinate_system = osr.SpatialReference()
    utm_coordinate_system.SetWellKnownGeogCS('WGS84') # Set geographic coordinate system to handle lat/lon
    utm_coordinate_system.SetUTM(utm_zone,is_northern)
    wgs84_coordinate_system = utm_coordinate_system.CloneGeogCS() # Clone ONLY the geographic coordinate system
    wgs84_to_utm_geo_transform = osr.CoordinateTransformation(wgs84_coordinate_system,utm_coordinate_system) # create transform component
    xyz = np.array(wgs84_to_utm_geo_transform.TransformPoints(np.dstack((longitude,latitude)).reshape((-1,2)))).reshape(longitude.shape[0],longitude.shape[1],3)
    return xyz[:,:,0],xyz[:,:,1],xyz[:,:,2] # returns easting, northing, altitude

if opts.add_coords:
    center_x = 107.268
    center_y = -6.839
    lon = np.arange(107+10/60,107+23/60,2.0/60.0)
    lat = np.arange(-6-56/60,-6.756,2.0/60.0)
    xg,yg = np.meshgrid(lon,lat)
    x,y,z = transform_wgs84_to_utm(xg,yg)
    ind_x = np.argmin(np.abs(lon-center_x))
    ind_y = np.argmin(np.abs(lat-center_y))
    center_x_utm = x[ind_y,:]
    center_y_utm = y[:,ind_x]
    x_labels = ['{:.0f}'.format(int(x))+'$^{\circ}$'+'{:.0f}'.format((x-int(x))*60.0)+'$^{\prime}$E' for x in lon]
    y_labels = ['{:.0f}'.format(int(y))+'$^{\circ}$'+'{:.0f}'.format((y-int(y))*60.0)+'$^{\prime}$S' for y in -lat]
    x_labels = ['{:d}'.format(int(x))+'$^{\circ}$'+'{:02d}'.format(int((x-int(x))*60.0+0.1))+'$^{\prime}$E' for x in lon]
    y_labels = ['{:d}'.format(int(y))+'$^{\circ}$'+'{:02d}'.format(int((y-int(y))*60.0+0.1))+'$^{\prime}$S' for y in -lat]

color = cm.hsv(np.linspace(0.0,1.0,365))
colors = np.vstack((color,color,color,color,color,color))
mymap = LinearSegmentedColormap.from_list('my_colormap',colors,N=len(colors)*2)

prj = ccrs.UTM(zone=48,southern_hemisphere=True)

block_shp = list(shpreader.Reader(opts.block_fnam).geometries())
block_rec = list(shpreader.Reader(opts.block_fnam).records())

ds = gdal.Open(opts.mask_fnam)
mask = ds.ReadAsArray()
mask_shape = mask.shape
ds = None

ds = gdal.Open(opts.trans_fnam)
data = ds.ReadAsArray()
data_trans = ds.GetGeoTransform()
data_shape = data[0].shape
ds = None
if data_shape != mask_shape:
    raise ValueError('Error, data_shape={}, mask_shape={}'.format(data_shape,mask_shape))
data[:,mask<0.5] = np.nan
#data[0] -= 9.0 # offset corrected
#data[0] += date2num(np.datetime64('0000-12-31')) # Matplotlib>3.3
xmin = data_trans[0]
xstp = data_trans[1]
xmax = xmin+xstp*data_shape[1]
ymax = data_trans[3]
ystp = data_trans[5]
ymin = ymax+ystp*data_shape[0]

sys.stderr.write('tmin: {}\n'.format(num2date(np.nanmin(data[0])).strftime('%Y%m%d')))
sys.stderr.write('tmax: {}\n'.format(num2date(np.nanmax(data[0])).strftime('%Y%m%d')))
sys.stderr.write('smin: {}\n'.format(np.nanmin(data[1])))
sys.stderr.write('smax: {}\n'.format(np.nanmax(data[1])))
if opts.tmin is not None:
    tmin = date2num(datetime.strptime(opts.tmin,'%Y%m%d'))
if opts.tmax is not None:
    tmax = date2num(datetime.strptime(opts.tmax,'%Y%m%d'))
if opts.smin is not None:
    smin = opts.smin
if opts.smax is not None:
    smax = opts.smax
tdif = tmax-tmin

values = []
labels = []
ticks = []
ds = tdif/365
for y in range(num2date(tmin).year,num2date(tmax).year+1):
    if ds > 2.0:
        for m in range(1,13,3):
            d = datetime(y,m,1)
            values.append(date2num(d))
            labels.append(d.strftime('%Y-%m'))
        for m in range(1,13,1):
            d = datetime(y,m,1)
            ticks.append(date2num(d))
    elif ds > 1.0:
        for m in range(1,13,2):
            d = datetime(y,m,1)
            values.append(date2num(d))
            labels.append(d.strftime('%Y-%m'))
        for m in range(1,13,1):
            d = datetime(y,m,1)
            ticks.append(date2num(d))
    else:
        for m in range(1,13,1):
            for day in [1,15]:
                d = datetime(y,m,day)
                values.append(date2num(d))
                labels.append(d.strftime('%m/%d'))
            for day in [5,10,20,25]:
                d = datetime(y,m,day)
                ticks.append(date2num(d))
dmin = num2date(tmin)
dmax = num2date(tmax)
if opts.add_tmin:
    if not tmin in values:
        if ds > 1.0:
            values.append(tmin)
            labels.append(dmin.strftime('%Y-%m'))
        else:
            values.append(tmin)
            labels.append(dmin.strftime('%m/%d'))
if opts.add_tmax:
    if not tmax in values:
        if ds > 1.0:
            values.append(tmax)
            labels.append(dmax.strftime('%Y-%m'))
        else:
            values.append(tmax)
            labels.append(dmax.strftime('%m/%d'))
torg = date2num(datetime(dmin.year,1,1))
twid = 365.0*2.0
newcolors = mymap(np.linspace((tmin-torg)/twid,(tmax-torg)/twid,mymap.N))
if opts.early:
    indx = int(mymap.N*0.995+0.5)
    newcolors[indx:,:] = to_rgba('maroon')
mymap2 = ListedColormap(newcolors)

if not opts.batch:
    plt.interactive(True)
fig = plt.figure(1,facecolor='w',figsize=(8.3,5.8))
plt.subplots_adjust(top=0.97,bottom=0.01,left=0.026,right=0.963,wspace=0.085,hspace=0.08)
fig.clear()

ax1 = plt.subplot(121,projection=prj)
ax2 = plt.subplot(122,projection=prj)

im1 = ax1.imshow(data[0],extent=(xmin,xmax,ymax,ymin),vmin=tmin,vmax=tmax,cmap=mymap2,interpolation='none')
ax12 = plt.colorbar(im1,ax=ax1,orientation='horizontal',shrink=1.0,pad=0.01).ax
ax12.xaxis.set_major_locator(plt.FixedLocator(values))
ax12.xaxis.set_major_formatter(plt.FixedFormatter(labels))
ax12.xaxis.set_minor_locator(plt.FixedLocator(ticks))
#ax1.set_title('(a) Transplanting date')
for l in ax12.xaxis.get_ticklabels():
    l.set_rotation(30)
ax12.set_xlabel('Estimated transplanting date (MM/DD)')
ax12.xaxis.set_label_coords(0.5,-2.8)
ax1.add_geometries(block_shp,prj,edgecolor='k',facecolor='none')

im2 = ax2.imshow(data[1],extent=(xmin,xmax,ymax,ymin),vmin=smin,vmax=smax,cmap=cm.jet,interpolation='none')
ax22 = plt.colorbar(im2,ax=ax2,orientation='horizontal',shrink=1.0,pad=0.01).ax
#ax22.xaxis.set_major_locator(plt.FixedLocator(values))
#ax22.xaxis.set_major_formatter(plt.FixedFormatter(labels))
#ax22.xaxis.set_minor_locator(plt.FixedLocator(ticks))
#ax12.set_xlabel('Transplanting date')
#ax2.set_title('(b) Signal')
#for l in ax22.xaxis.get_ticklabels():
#    l.set_rotation(30)
ax22.set_xlabel('Signal (dB)')
ax22.xaxis.set_label_coords(0.5,-2.6)
ax2.add_geometries(block_shp,prj,edgecolor='k',facecolor='none')

if opts.add_coords:
    ax1.xaxis.set_tick_params(labelsize=6,direction='in',pad=2)
    ax1.yaxis.set_tick_params(labelsize=6,direction='in',pad=2)
    ax1.xaxis.set_label_position('top')
    ax1.xaxis.tick_top()
    ax1.yaxis.set_label_position('right')
    ax1.yaxis.tick_right()
    plt.setp(ax1.get_xticklabels(),color=opts.coords_color)
    plt.setp(ax1.get_yticklabels(),color=opts.coords_color)
    ax1.set_xticks(np.arange(100.0,120.0,0.1))
    ax1.set_yticks(np.arange(-7.5,-5.5,0.1))
    ax1.xaxis.set_major_locator(plt.FixedLocator(center_x_utm))
    ax1.yaxis.set_major_locator(plt.FixedLocator(center_y_utm))
    ax1.xaxis.set_major_formatter(plt.FixedFormatter(x_labels))
    ax1.yaxis.set_major_formatter(plt.FixedFormatter(y_labels))

    ax2.xaxis.set_tick_params(labelsize=6,direction='in',pad=2)
    ax2.yaxis.set_tick_params(labelsize=6,direction='in',pad=2)
    ax2.xaxis.set_label_position('top')
    ax2.xaxis.tick_top()
    ax2.yaxis.set_label_position('right')
    ax2.yaxis.tick_right()
    plt.setp(ax2.get_xticklabels(),color=opts.coords_color)
    plt.setp(ax2.get_yticklabels(),color=opts.coords_color)
    ax2.set_xticks(np.arange(100.0,120.0,0.1))
    ax2.set_yticks(np.arange(-7.5,-5.5,0.1))
    ax2.xaxis.set_major_locator(plt.FixedLocator(center_x_utm))
    ax2.yaxis.set_major_locator(plt.FixedLocator(center_y_utm))
    ax2.xaxis.set_major_formatter(plt.FixedFormatter(x_labels))
    ax2.yaxis.set_major_formatter(plt.FixedFormatter(y_labels))

ax1.set_xlim(xmin,xmax)
ax1.set_ylim(ymin,ymax)
ax2.set_xlim(xmin,xmax)
ax2.set_ylim(ymin,ymax)
#ax1.background_patch.set_facecolor('#eeeeee')
#ax3.background_patch.set_facecolor('#eeeeee')
if opts.title is not None:
    plt.suptitle(opts.title)

plt.savefig(opts.output_fnam,dpi=1200)
if not opts.batch:
    plt.draw()
