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
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap,LinearSegmentedColormap,to_rgba
from matplotlib.dates import date2num,num2date
from matplotlib.path import Path
from optparse import OptionParser,IndentedHelpFormatter

# Default values
HOME = os.environ.get('HOME')
if HOME is None:
    HOME = os.environ.get('HOMEPATH')
TMIN = '20190415'
TMAX = '20190601'
PMIN = 0.0
PMAX = 300.0
COORDS_COLOR = '#aaaaaa'
BLOCK_FNAM = os.path.join(HOME,'Work','SATREPS','Shapefile','cihea_testsite_200819','cihea_testsite_200819.shp')
TRANS_FNAM = os.path.join('.','transplanting_date.shp')
OUTPUT_FNAM = 'trans_date_testsite.pdf'

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('-s','--tmin',default=TMIN,help='Min date in the format YYYYMMDD (%default)')
parser.add_option('-e','--tmax',default=TMAX,help='Max date in the format YYYYMMDD (%default)')
parser.add_option('-p','--pmin',default=PMIN,type='float',help='Min signal in dB (%default)')
parser.add_option('-P','--pmax',default=PMAX,type='float',help='Max signal in dB (%default)')
parser.add_option('-t','--title',default=None,help='Figure title (%default)')
parser.add_option('--block_fnam',default=BLOCK_FNAM,help='Block shape file (%default)')
parser.add_option('--trans_fnam',default=TRANS_FNAM,help='Transplanting shape file (%default)')
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
    lon = np.arange(107+10/60,107+23/60,2.0/360.0)
    lat = np.arange(-6-56/60,-6.756,2.0/360.0)
    xg,yg = np.meshgrid(lon,lat)
    x,y,z = transform_wgs84_to_utm(xg,yg)
    ind_x = np.argmin(np.abs(lon-center_x))
    ind_y = np.argmin(np.abs(lat-center_y))
    center_x_utm = x[ind_y,:]
    center_y_utm = y[:,ind_x]
    x_labels = ['{:d}'.format(int(x))+'$^{\circ}$'+'{:02d}'.format(int((x-int(x))*60.0+0.1))+'$^{\prime}$'+'{:02d}'.format(int((x*60.0-int(x*60.0))*60.0+0.1))+'$^{\prime\prime}$E' for x in lon]
    y_labels = ['{:d}'.format(int(y))+'$^{\circ}$'+'{:02d}'.format(int((y-int(y))*60.0+0.1))+'$^{\prime}$'+'{:02d}'.format(int((y*60.0-int(y*60.0))*60.0+0.1))+'$^{\prime\prime}$S' for y in -lat]

color = cm.hsv(np.linspace(0.0,1.0,365))
colors = np.vstack((color,color,color,color,color,color))
mymap = LinearSegmentedColormap.from_list('my_colormap',colors,N=len(colors)*2)

prj = ccrs.UTM(zone=48,southern_hemisphere=True)

block_shp = list(shpreader.Reader(opts.block_fnam).geometries())
block_rec = list(shpreader.Reader(opts.block_fnam).records())

shapes = list(shpreader.Reader(opts.trans_fnam).geometries())
records = list(shpreader.Reader(opts.trans_fnam).records())
xmin = 1.0e10
xmax = -1.0e10
ymin = 1.0e10
ymax = -1.0e10
for p in block_shp:
    x1,y1,x2,y2 = p.bounds
    if x1 < xmin:
        xmin = x1
    if x2 > xmax:
        xmax = x2
    if y1 < ymin:
        ymin = y1
    if y2 > ymax:
        ymax = y2
xmin -= 10.0
xmax += 10.0
ymin -= 10.0
ymax += 10.0

xp = np.arange(xmin,xmax,10.0)
yp = np.arange(ymin,ymax,10.0)
xg,yg = np.meshgrid(xp,yp)

tmin = 1.0e10
tmax = -1.0e10
pmin = 1.0e10
pmax = -1.0e10
for rec in records:
    t = rec.attributes['trans_date']#+date2num(np.datetime64('0000-12-31'))
    p = rec.attributes['peak_value']
    if t > 1000.0:
        if t < tmin:
            tmin = t
        if t > tmax:
            tmax = t
        if p < pmin:
            pmin = p
        if p > pmax:
            pmax = p
sys.stderr.write('tmin: {}\n'.format(num2date(tmin).strftime('%Y%m%d')))
sys.stderr.write('tmax: {}\n'.format(num2date(tmax).strftime('%Y%m%d')))
sys.stderr.write('pmin: {}\n'.format(pmin))
sys.stderr.write('pmax: {}\n'.format(pmax))
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
fig = plt.figure(1,facecolor='w',figsize=(7.4,6.0))
fig.clear()
plt.subplots_adjust(top=0.93,bottom=0.01,left=0.030,right=0.940,wspace=0.145,hspace=0.0)

ax1 = plt.subplot(121,projection=prj)
ax2 = plt.subplot(122,projection=prj)

for shp,rec in zip(shapes,records):
    t = rec.attributes['trans_date']#-9.0#+date2num(np.datetime64('0000-12-31')) # offset corrected
    p = rec.attributes['peak_value']
    if not np.isnan(t):
        ax1.add_geometries(shp,prj,edgecolor='none',facecolor=mymap2((t-tmin)/tdif))
    if not np.isnan(p):
        ax2.add_geometries(shp,prj,edgecolor='none',facecolor=cm.jet((p-pmin)/pdif))
im1 = ax1.imshow(np.arange(4).reshape(2,2),extent=(-2,-1,-2,-1),vmin=tmin,vmax=tmax,cmap=mymap2)
ax12 = plt.colorbar(im1,ax=ax1,orientation='horizontal',shrink=1.0,pad=0.01).ax
ax12.xaxis.set_major_locator(plt.FixedLocator(values))
ax12.xaxis.set_major_formatter(plt.FixedFormatter(labels))
ax12.xaxis.set_minor_locator(plt.FixedLocator(ticks))
#ax1.set_title('(a)')
for l in ax12.xaxis.get_ticklabels():
    l.set_rotation(30)
ax12.xaxis.set_label_coords(0.5,-3.2)
ax12.set_xlabel('Estimated transplanting date (MM/DD)')
ax1.add_geometries(block_shp,prj,edgecolor='k',facecolor='none')

im2 = ax2.imshow(np.arange(4).reshape(2,2),extent=(-2,-1,-2,-1),vmin=pmin,vmax=pmax,cmap=cm.jet)
ax22 = plt.colorbar(im2,ax=ax2,orientation='horizontal',shrink=1.0,pad=0.01).ax
ax22.minorticks_on()
#ax2.set_title('(b)')
ax22.set_xlabel('Signal (dB)')
ax22.xaxis.set_label_coords(0.5,-3.2)
ax2.add_geometries(block_shp,prj,edgecolor='k',facecolor='none')

for bs,br in zip(block_shp,block_rec):
    for n in range(len(bs)):
        #xc = bs[0].boundary.centroid.x
        #yc = bs[0].boundary.centroid.y
        points = list(zip(*bs[n].exterior.coords.xy))
        p = Path(points)
        flags = p.contains_points(np.hstack((xg.flatten()[:,np.newaxis],yg.flatten()[:,np.newaxis]))).reshape(xg.shape)
        distance = np.empty(xg.shape,dtype=np.float64)
        for i,x_pt in enumerate(xp):
            for j,y_pt in enumerate(yp):
                if flags[j,i]:
                    distance[j,i] = bs[n].exterior.distance(shapely.geometry.Point(x_pt,y_pt))
                else:
                    distance[j,i] = -1.0e10
        indy,indx = np.unravel_index(np.argmax(distance),distance.shape)
        if indx > 0 and indy > 0:
            xt = xp[indx]
            yt = yp[indy]
            t = br.attributes['Blok']
            ax1.text(xt,yt,t,transform=prj,ha='center',va='center')
            ax2.text(xt,yt,t,transform=prj,ha='center',va='center')
        #break

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
if opts.title is not None:
    plt.suptitle(opts.title,y=0.99)
plt.savefig(opts.output_fnam)
if not opts.batch:
    plt.draw()
