#!/usr/bin/env python
#from shapely.geometry.polygon import Polygon
import os
import sys
import warnings
from datetime import datetime
import numpy as np
import osr
import shapefile
import shapely
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
from optparse import OptionParser,IndentedHelpFormatter

# Default values
HOME = os.environ.get('HOME')
if HOME is None:
    HOME = os.environ.get('HOMEPATH')
TMIN = '20190501'
TMAX = '20190915'
PMIN = 0.0
PMAX = 0.8
MMIN = 0
MMAX = 5
NMIN = 0
NMAX = 12
XMGN = 100.0
YMGN = 100.0
NCAN = 1
COORDS_COLOR = '#aaaaaa'
TRANS_FNAM = os.path.join('.','transplanting_date.shp')
OUTPUT_FNAM = 'trans_date_testsite.pdf'

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('--site',default=None,help='Site name (%default)')
parser.add_option('-s','--tmin',default=TMIN,help='Min date in the format YYYYMMDD (%default)')
parser.add_option('-e','--tmax',default=TMAX,help='Max date in the format YYYYMMDD (%default)')
parser.add_option('-p','--pmin',default=PMIN,type='float',help='Min NDVI (%default)')
parser.add_option('-P','--pmax',default=PMAX,type='float',help='Max NDVI (%default)')
parser.add_option('--mmin',default=MMIN,type='int',help='Min #data within 5 days (%default)')
parser.add_option('--mmax',default=MMAX,type='int',help='Max #data within 5 days (%default)')
parser.add_option('--nmin',default=NMIN,type='int',help='Min #data within 15 days (%default)')
parser.add_option('--nmax',default=NMAX,type='int',help='Max #data within 15 days (%default)')
parser.add_option('--xmgn',default=XMGN,type='float',help='X margin in m (%default)')
parser.add_option('--ymgn',default=YMGN,type='float',help='Y margin in m (%default)')
parser.add_option('-N','--ncan',default=NCAN,type='int',help='Candidate number between 1 and 3 (%default)')
parser.add_option('-t','--title',default=None,help='Figure title (%default)')
parser.add_option('--trans_fnam',default=TRANS_FNAM,help='Transplanting shape file (%default)')
parser.add_option('--outline_fnam',default=None,help='Outline shapefile name (%default)')
parser.add_option('--output_fnam',default=OUTPUT_FNAM,help='Output figure name (%default)')
parser.add_option('--add_tmin',default=False,action='store_true',help='Add tmin in colorbar (%default)')
parser.add_option('--add_tmax',default=False,action='store_true',help='Add tmax in colorbar (%default)')
parser.add_option('--add_coords',default=False,action='store_true',help='Add geographical coordinates (%default)')
parser.add_option('--coords_color',default=COORDS_COLOR,help='Color of geographical coordinates (%default)')
parser.add_option('--early',default=False,action='store_true',help='Early estimation mode (%default)')
parser.add_option('-b','--batch',default=False,action='store_true',help='Batch mode (%default)')
parser.add_option('--debug',default=False,action='store_true',help='Debug mode (%default)')
(opts,args) = parser.parse_args()
if opts.ncan < 1 or opts.ncan > 3:
    raise ValueError('Error, ncan={}'.format(opts.ncan))
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
    center_x = 107.67225
    center_y = -6.98795
    lon = np.arange(107+36/60,107+44/60,1.0/60.0)
    lat = np.arange(-7-2/60,-6-56/60,1.0/60.0)
    xg,yg = np.meshgrid(lon,lat)
    x,y,z = transform_wgs84_to_utm(xg,yg)
    ind_x = np.argmin(np.abs(lon-center_x))
    ind_y = np.argmin(np.abs(lat-center_y))
    center_x_utm = x[ind_y,:]
    center_y_utm = y[:,ind_x]
    #x_labels = ['{:d}'.format(int(x))+'$^{\circ}$'+'{:02d}'.format(int((x-int(x))*60.0+0.1))+'$^{\prime}$'+'{:02d}'.format(int((x*60.0-int(x*60.0))*60.0+0.1))+'$^{\prime\prime}$E' for x in lon]
    #y_labels = ['{:d}'.format(int(y))+'$^{\circ}$'+'{:02d}'.format(int((y-int(y))*60.0+0.1))+'$^{\prime}$'+'{:02d}'.format(int((y*60.0-int(y*60.0))*60.0+0.1))+'$^{\prime\prime}$S' for y in -lat]
    x_labels = ['{:d}'.format(int(x))+'$^{\circ}$'+'{:02d}'.format(int((x-int(x))*60.0+0.1))+'$^{\prime}$E' for x in lon]
    y_labels = ['{:d}'.format(int(y))+'$^{\circ}$'+'{:02d}'.format(int((y-int(y))*60.0+0.1))+'$^{\prime}$S' for y in -lat]

color = cm.hsv(np.linspace(0.0,1.0,365))
colors = np.vstack((color,color,color,color,color,color))
mymap = LinearSegmentedColormap.from_list('my_colormap',colors,N=len(colors)*2)
cmap5 = cm.get_cmap('jet',opts.mmax-opts.mmin+1)
cmap15 = cm.get_cmap('jet',opts.nmax-opts.nmin+1)

prj = ccrs.UTM(zone=48,southern_hemisphere=True)

shapes = list(shpreader.Reader(opts.trans_fnam).geometries())
records = list(shpreader.Reader(opts.trans_fnam).records())
if opts.outline_fnam is not None:
    outline_shapes = list(shpreader.Reader(opts.outline_fnam).geometries())
    pp = outline_shapes
else:
    pp = shapes
xmin = 1.0e10
xmax = -1.0e10
ymin = 1.0e10
ymax = -1.0e10
for p in pp:
    x1,y1,x2,y2 = p.bounds
    if x1 < xmin:
        xmin = x1
    if x2 > xmax:
        xmax = x2
    if y1 < ymin:
        ymin = y1
    if y2 > ymax:
        ymax = y2
xmin -= opts.xmgn
xmax += opts.xmgn
ymin -= opts.ymgn
ymax += opts.ymgn

tmin = 1.0e10
tmax = -1.0e10
pmin = 1.0e10
pmax = -1.0e10
mmin = 1.0e10
mmax = -1.0e10
nmin = 1.0e10
nmax = -1.0e10
for rec in records:
    t = rec.attributes['trans_d{:d}'.format(opts.ncan)]#+date2num(np.datetime64('0000-12-31'))
    p = rec.attributes['signal_{:d}'.format(opts.ncan)]
    m = rec.attributes['ndat5_{:d}'.format(opts.ncan)]
    n = rec.attributes['ndat15_{:d}'.format(opts.ncan)]
    if t > 1000.0:
        if t < tmin:
            tmin = t
        if t > tmax:
            tmax = t
        if p < pmin:
            pmin = p
        if p > pmax:
            pmax = p
        if m < mmin:
            mmin = m
        if m > mmax:
            mmax = m
        if n < nmin:
            nmin = n
        if n > nmax:
            nmax = n
sys.stderr.write('tmin: {}\n'.format(num2date(tmin).strftime('%Y%m%d')))
sys.stderr.write('tmax: {}\n'.format(num2date(tmax).strftime('%Y%m%d')))
sys.stderr.write('pmin: {}\n'.format(pmin))
sys.stderr.write('pmax: {}\n'.format(pmax))
sys.stderr.write('mmin: {}\n'.format(mmin))
sys.stderr.write('mmax: {}\n'.format(mmax))
sys.stderr.write('nmin: {}\n'.format(nmin))
sys.stderr.write('nmax: {}\n'.format(nmax))
if opts.tmin is not None:
    tmin = date2num(datetime.strptime(opts.tmin,'%Y%m%d'))
if opts.tmax is not None:
    tmax = date2num(datetime.strptime(opts.tmax,'%Y%m%d'))
if opts.pmin is not None:
    pmin = opts.pmin
if opts.pmax is not None:
    pmax = opts.pmax
tdif = tmax-tmin
pdif = pmax-pmin
qmin = 0.0
qmax = 0.5
qdif = qmax-qmin
smin = -0.3
smax = 0.3
sdif = smax-smin
mmin = opts.mmin-0.5
mmax = opts.mmax+0.5
mdif = mmax-mmin
nmin = opts.nmin-0.5
nmax = opts.nmax+0.5
ndif = nmax-nmin

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
            for day in [1]:
                d = datetime(y,m,day)
                values.append(date2num(d))
                labels.append(d.strftime('%m/%d'))
            for day in [5,10,15,20,25]:
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

site_low = opts.site.lower()
if site_low == 'cihea':
    fig = plt.figure(1,facecolor='w',figsize=(7.6,14.0))
    plt.subplots_adjust(top=0.965,bottom=0.005,left=0.030,right=0.96,wspace=0.12,hspace=0.043)
elif site_low == 'bojongsoang':
    fig = plt.figure(1,facecolor='w',figsize=(7.6,8.6))
    plt.subplots_adjust(top=0.94,bottom=0.06,left=0.030,right=0.96,wspace=0.12,hspace=0.25)
else:
    raise ValueError('Error in site >>> '+opts.site)
fig.clear()

ax1 = plt.subplot(321,projection=prj)
ax2 = plt.subplot(322,projection=prj)
ax3 = plt.subplot(323,projection=prj)
ax4 = plt.subplot(324,projection=prj)
ax5 = plt.subplot(325,projection=prj)
ax6 = plt.subplot(326,projection=prj)

for shp,rec in zip(shapes,records):
    t = rec.attributes['trans_d{:d}'.format(opts.ncan)]#-9.0#+date2num(np.datetime64('0000-12-31')) # offset corrected
    p = rec.attributes['ndvi_{:d}'.format(opts.ncan)]
    q = rec.attributes['ndvimax{:d}'.format(opts.ncan)]
    s = rec.attributes['signal_{:d}'.format(opts.ncan)]
    m = rec.attributes['ndat5_{:d}'.format(opts.ncan)]
    n = rec.attributes['ndat15_{:d}'.format(opts.ncan)]
    if not np.isnan(t):
        ax1.add_geometries(shp,prj,edgecolor='none',facecolor=mymap2((t-tmin)/tdif))
    if not np.isnan(p):
        ax2.add_geometries(shp,prj,edgecolor='none',facecolor=cm.jet((p-pmin)/pdif))
    if not np.isnan(q):
        ax3.add_geometries(shp,prj,edgecolor='none',facecolor=cm.jet((q-qmin)/qdif))
    if not np.isnan(s):
        ax4.add_geometries(shp,prj,edgecolor='none',facecolor=cm.jet((s-smin)/sdif))
    if not np.isnan(t):
        ax5.add_geometries(shp,prj,edgecolor='none',facecolor=cmap5((m-mmin)/mdif))
    if not np.isnan(t):
        ax6.add_geometries(shp,prj,edgecolor='none',facecolor=cmap15((n-nmin)/ndif))
im1 = ax1.imshow(np.arange(4).reshape(2,2),extent=(-2,-1,-2,-1),vmin=tmin,vmax=tmax,cmap=mymap2)
ax12 = plt.colorbar(im1,ax=ax1,orientation='horizontal',shrink=1.00,pad=0.01).ax
ax12.xaxis.set_major_locator(plt.FixedLocator(values))
ax12.xaxis.set_major_formatter(plt.FixedFormatter(labels))
ax12.xaxis.set_minor_locator(plt.FixedLocator(ticks))
#ax1.set_title('(a)')
for l in ax12.xaxis.get_ticklabels():
    l.set_rotation(30)
ax12.xaxis.set_label_coords(0.5,-3.0)
ax12.set_xlabel('Estimated transplanting date (MM/DD)')
if opts.outline_fnam is not None:
    ax1.add_geometries(outline_shapes,prj,edgecolor='k',facecolor='none')

im2 = ax2.imshow(np.arange(4).reshape(2,2),extent=(-2,-1,-2,-1),vmin=pmin,vmax=pmax,cmap=cm.jet)
ax22 = plt.colorbar(im2,ax=ax2,orientation='horizontal',shrink=1.00,pad=0.01).ax
#ax22.xaxis.set_major_locator(plt.MultipleLocator(5.0))
ax22.minorticks_on()
#ax2.set_title('(b)')
ax22.set_xlabel('NDVI at transplanting')
ax22.xaxis.set_label_coords(0.5,-3.0)
if opts.outline_fnam is not None:
    ax2.add_geometries(outline_shapes,prj,edgecolor='k',facecolor='none')

im3 = ax3.imshow(np.arange(4).reshape(2,2),extent=(-2,-1,-2,-1),vmin=qmin,vmax=qmax,cmap=cm.jet)
ax32 = plt.colorbar(im3,ax=ax3,orientation='horizontal',shrink=1.00,pad=0.01).ax
ax32.minorticks_on()
#ax3.set_title('(c)')
ax32.set_xlabel('NDVI increase after transplanting')
ax32.xaxis.set_label_coords(0.5,-2.2)
if opts.outline_fnam is not None:
    ax3.add_geometries(outline_shapes,prj,edgecolor='k',facecolor='none')

im4 = ax4.imshow(np.arange(4).reshape(2,2),extent=(-2,-1,-2,-1),vmin=smin,vmax=smax,cmap=cm.jet)
ax42 = plt.colorbar(im4,ax=ax4,orientation='horizontal',shrink=1.00,pad=0.01).ax
ax42.minorticks_on()
#ax4.set_title('(d)')
ax42.set_xlabel('Signal at transplanting')
ax42.xaxis.set_label_coords(0.5,-2.2)
if opts.outline_fnam is not None:
    ax4.add_geometries(outline_shapes,prj,edgecolor='k',facecolor='none')

im5 = ax5.imshow(np.arange(4).reshape(2,2),extent=(-2,-1,-2,-1),vmin=mmin,vmax=mmax,cmap=cmap5)
ax52 = plt.colorbar(im5,ax=ax5,orientation='horizontal',ticks=np.arange(opts.mmin,opts.mmax+1,1),shrink=1.00,pad=0.01).ax
#ax52.minorticks_on()
#ax5.set_title('(e)')
ax52.set_xlabel('#Data within $\pm$5 days')
ax52.xaxis.set_label_coords(0.5,-2.2)
if opts.outline_fnam is not None:
    ax5.add_geometries(outline_shapes,prj,edgecolor='k',facecolor='none')

im6 = ax6.imshow(np.arange(4).reshape(2,2),extent=(-2,-1,-2,-1),vmin=nmin,vmax=nmax,cmap=cmap15)
ax62 = plt.colorbar(im6,ax=ax6,orientation='horizontal',ticks=np.arange(opts.nmin,opts.nmax+1,1),shrink=1.00,pad=0.01).ax
#ax62.minorticks_on()
#ax6.set_title('(f)')
ax62.set_xlabel('#Data within $\pm$15 days')
ax62.xaxis.set_label_coords(0.5,-2.2)
if opts.outline_fnam is not None:
    ax6.add_geometries(outline_shapes,prj,edgecolor='k',facecolor='none')

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

    ax3.xaxis.set_tick_params(labelsize=6,direction='in',pad=2)
    ax3.yaxis.set_tick_params(labelsize=6,direction='in',pad=2)
    ax3.xaxis.set_label_position('top')
    ax3.xaxis.tick_top()
    ax3.yaxis.set_label_position('right')
    ax3.yaxis.tick_right()
    plt.setp(ax3.get_xticklabels(),color=opts.coords_color)
    plt.setp(ax3.get_yticklabels(),color=opts.coords_color)
    ax3.set_xticks(np.arange(100.0,120.0,0.1))
    ax3.set_yticks(np.arange(-7.5,-5.5,0.1))
    ax3.xaxis.set_major_locator(plt.FixedLocator(center_x_utm))
    ax3.yaxis.set_major_locator(plt.FixedLocator(center_y_utm))
    ax3.xaxis.set_major_formatter(plt.FixedFormatter(x_labels))
    ax3.yaxis.set_major_formatter(plt.FixedFormatter(y_labels))

    ax4.xaxis.set_tick_params(labelsize=6,direction='in',pad=2)
    ax4.yaxis.set_tick_params(labelsize=6,direction='in',pad=2)
    ax4.xaxis.set_label_position('top')
    ax4.xaxis.tick_top()
    ax4.yaxis.set_label_position('right')
    ax4.yaxis.tick_right()
    plt.setp(ax4.get_xticklabels(),color=opts.coords_color)
    plt.setp(ax4.get_yticklabels(),color=opts.coords_color)
    ax4.set_xticks(np.arange(100.0,120.0,0.1))
    ax4.set_yticks(np.arange(-7.5,-5.5,0.1))
    ax4.xaxis.set_major_locator(plt.FixedLocator(center_x_utm))
    ax4.yaxis.set_major_locator(plt.FixedLocator(center_y_utm))
    ax4.xaxis.set_major_formatter(plt.FixedFormatter(x_labels))
    ax4.yaxis.set_major_formatter(plt.FixedFormatter(y_labels))

    ax5.xaxis.set_tick_params(labelsize=6,direction='in',pad=2)
    ax5.yaxis.set_tick_params(labelsize=6,direction='in',pad=2)
    ax5.xaxis.set_label_position('top')
    ax5.xaxis.tick_top()
    ax5.yaxis.set_label_position('right')
    ax5.yaxis.tick_right()
    plt.setp(ax5.get_xticklabels(),color=opts.coords_color)
    plt.setp(ax5.get_yticklabels(),color=opts.coords_color)
    ax5.set_xticks(np.arange(100.0,120.0,0.1))
    ax5.set_yticks(np.arange(-7.5,-5.5,0.1))
    ax5.xaxis.set_major_locator(plt.FixedLocator(center_x_utm))
    ax5.yaxis.set_major_locator(plt.FixedLocator(center_y_utm))
    ax5.xaxis.set_major_formatter(plt.FixedFormatter(x_labels))
    ax5.yaxis.set_major_formatter(plt.FixedFormatter(y_labels))

    ax6.xaxis.set_tick_params(labelsize=6,direction='in',pad=2)
    ax6.yaxis.set_tick_params(labelsize=6,direction='in',pad=2)
    ax6.xaxis.set_label_position('top')
    ax6.xaxis.tick_top()
    ax6.yaxis.set_label_position('right')
    ax6.yaxis.tick_right()
    plt.setp(ax6.get_xticklabels(),color=opts.coords_color)
    plt.setp(ax6.get_yticklabels(),color=opts.coords_color)
    ax6.set_xticks(np.arange(100.0,120.0,0.1))
    ax6.set_yticks(np.arange(-7.5,-5.5,0.1))
    ax6.xaxis.set_major_locator(plt.FixedLocator(center_x_utm))
    ax6.yaxis.set_major_locator(plt.FixedLocator(center_y_utm))
    ax6.xaxis.set_major_formatter(plt.FixedFormatter(x_labels))
    ax6.yaxis.set_major_formatter(plt.FixedFormatter(y_labels))

ax1.set_xlim(xmin,xmax)
ax1.set_ylim(ymin,ymax)
ax2.set_xlim(xmin,xmax)
ax2.set_ylim(ymin,ymax)
ax3.set_xlim(xmin,xmax)
ax3.set_ylim(ymin,ymax)
ax4.set_xlim(xmin,xmax)
ax4.set_ylim(ymin,ymax)
ax5.set_xlim(xmin,xmax)
ax5.set_ylim(ymin,ymax)
ax6.set_xlim(xmin,xmax)
ax6.set_ylim(ymin,ymax)
if opts.title is not None:
    plt.suptitle(opts.title)
plt.savefig(opts.output_fnam)
if not opts.batch:
    plt.show()
