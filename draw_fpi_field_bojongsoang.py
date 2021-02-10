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
COORDS_COLOR = '#aaaaaa'
TRANS_FNAM = os.path.join('.','transplanting_date.shp')
OUTPUT_FNAM = 'fpi_bojongsoang.pdf'

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('-t','--title',default=None,help='Figure title (%default)')
parser.add_option('--trans_fnam',default=TRANS_FNAM,help='Transplanting shape file (%default)')
parser.add_option('--output_fnam',default=OUTPUT_FNAM,help='Output figure name (%default)')
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

prj = ccrs.UTM(zone=48,southern_hemisphere=True)

shapes = list(shpreader.Reader(opts.trans_fnam).geometries())
records = list(shpreader.Reader(opts.trans_fnam).records())
xmin = 1.0e10
xmax = -1.0e10
ymin = 1.0e10
ymax = -1.0e10
for p in shapes:
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

fmin = 0.0
fmax = 1.0
fdif = fmax-fmin

fig = plt.figure(1,facecolor='w',figsize=(4.0,3.0))
fig.clear()
plt.subplots_adjust(top=0.97,bottom=0.09,left=0.034,right=0.96,wspace=0.12,hspace=0.25)

ax1 = plt.subplot(111,projection=prj)

for shp,rec in zip(shapes,records):
    f = rec.attributes['fpi_e']
    ax1.add_geometries(shp,prj,edgecolor='none',facecolor=cm.jet((f-fmin)/fdif))

im4 = ax1.imshow(np.arange(4).reshape(2,2),extent=(-2,-1,-2,-1),vmin=fmin,vmax=fmax,cmap=cm.jet)
ax12 = plt.colorbar(im4,ax=ax1,orientation='horizontal',shrink=1.0,pad=0.01).ax
ax12.minorticks_on()
#ax2.set_title('(b)')
ax12.set_xlabel('Fishpond Index')
ax12.xaxis.set_label_coords(0.5,-1.8)
#ax2.add_geometries(shapes,prj,edgecolor='k',facecolor='none')

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

ax1.set_xlim(xmin,xmax)
ax1.set_ylim(ymin,ymax)
if opts.title is not None:
    plt.suptitle(opts.title)
plt.savefig(opts.output_fnam)
if not opts.batch:
    plt.show()
