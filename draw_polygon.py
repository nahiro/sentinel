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
FIELD_NAME = 'value'
FIELD_TYPE = 'float'
HOME = os.environ.get('HOME')
if HOME is None:
    HOME = os.environ.get('HOMEPATH')
SHP_FNAM = os.path.join(HOME,'Work','SATREPS','Shapefile','field_GIS','Bojongsoang','Bojongsoang')
OUTPUT_FNAM = 'output.pdf'
COORDS_COLOR = '#aaaaaa'

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('--site',default=None,help='Site name (%default)')
parser.add_option('-a','--field_name',default=FIELD_NAME,help='Field name to draw (%default)')
parser.add_option('-t','--field_type',default=FIELD_TYPE,help='Field type, time, float, or number (%default)')
parser.add_option('-z','--zmin',default=None,help='Min value (%default)')
parser.add_option('-Z','--zmax',default=None,help='Max value (%default)')
parser.add_option('-i','--inp_fnam',default=None,help='Input file name (%default)')
parser.add_option('--shp_fnam',default=SHP_FNAM,help='Input shapefile name (%default)')
parser.add_option('--outline_fnam',default=None,help='Outline shapefile name (%default)')
parser.add_option('--output_fnam',default=OUTPUT_FNAM,help='Output figure name (%default)')
parser.add_option('-T','--title',default=None,help='Figure title (%default)')
parser.add_option('--zlabel',default=None,help='Colorbar label (%default)')
parser.add_option('--add_tmin',default=False,action='store_true',help='Add tmin in colorbar (%default)')
parser.add_option('--add_tmax',default=False,action='store_true',help='Add tmax in colorbar (%default)')
parser.add_option('--add_coords',default=False,action='store_true',help='Add geographical coordinates (%default)')
parser.add_option('--coords_color',default=COORDS_COLOR,help='Color of geographical coordinates (%default)')
parser.add_option('-b','--batch',default=False,action='store_true',help='Batch mode (%default)')
parser.add_option('--debug',default=False,action='store_true',help='Debug mode (%default)')
(opts,args) = parser.parse_args()
if not opts.debug:
    warnings.simplefilter('ignore')
field_type = opts.field_type[0].upper()
if not field_type.upper() in ['T','F','N']:
    raise ValueError('Error in field type >>> '+opts.field_type)

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

shapes = list(shpreader.Reader(opts.shp_fnam).geometries())
records = list(shpreader.Reader(opts.shp_fnam).records())
nobject = len(records)

if opts.outline_fnam is not None:
    outline_shapes = list(shpreader.Reader(opts.outline_fnam).geometries())

ext = os.path.splitext(opts.inp_fnam)[1].lower()
if ext == '.npz':
    data = np.load(opts.inp_fnam)[opts.field_name]
elif ext == '.npy':
    data = np.load(opts.inp_fnam)
else:
    data = []
    for rec in list(shpreader.Reader(opts.inp_fnam).records()):
        data.append(rec.attributes[opts.field_name])
if len(data) != nobject:
    raise ValueError('Error, len(data)={}, nobject={}'.format(len(data),nobject))

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
if opts.outline_fnam is not None:
    for p in outline_shapes:
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

zmin = 1.0e10
zmax = -1.0e10
for z in data:
    if z < zmin:
        zmin = z
    if z > zmax:
        zmax = z
if field_type == 'T':
    sys.stderr.write('zmin: {}\n'.format(num2date(zmin).strftime('%Y%m%d')))
    sys.stderr.write('zmax: {}\n'.format(num2date(zmax).strftime('%Y%m%d')))
else:
    sys.stderr.write('zmin: {}\n'.format(zmin))
    sys.stderr.write('zmax: {}\n'.format(zmax))
if opts.zmin is not None:
    if field_type == 'T':
        zmin = date2num(datetime.strptime(opts.zmin,'%Y%m%d'))
    elif field_type == 'N':
        zmin = int(opts.zmin)
    else:
        zmin = float(opts.zmin)
if opts.zmax is not None:
    if field_type == 'T':
        zmax = date2num(datetime.strptime(opts.zmax,'%Y%m%d'))
    elif field_type == 'N':
        zmax = int(opts.zmax)
    else:
        zmax = float(opts.zmax)
zdif = zmax-zmin

if field_type == 'T':
    dmin = num2date(zmin)
    dmax = num2date(zmax)
    color = cm.hsv(np.linspace(0.0,1.0,365))
    colors = np.vstack((color,color,color,color,color,color))
    mymap = LinearSegmentedColormap.from_list('my_colormap',colors,N=len(colors)*2)
    values = []
    labels = []
    ticks = []
    ds = tdif/365
    for y in range(dmin.year,dmax.year+1):
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
    if opts.add_tmin:
        if not zmin in values:
            if ds > 1.0:
                values.append(zmin)
                labels.append(dmin.strftime('%Y-%m'))
            else:
                values.append(zmin)
                labels.append(dmin.strftime('%m/%d'))
    if opts.add_tmax:
        if not zmax in values:
            if ds > 1.0:
                values.append(zmax)
                labels.append(dmax.strftime('%Y-%m'))
            else:
                values.append(zmax)
                labels.append(dmax.strftime('%m/%d'))
    torg = date2num(datetime(dmin.year,1,1))
    twid = 365.0*2.0
    newcolors = mymap(np.linspace((zmin-torg)/twid,(zmax-torg)/twid,mymap.N))
    mycmap = ListedColormap(newcolors)
else:
    mycmap = cm.jet

site_low = opts.site.lower()
if site_low == 'cihea':
    fig = plt.figure(1,facecolor='w',figsize=(8.0,5.3))
    plt.subplots_adjust(top=0.99,bottom=0.00,left=0.044,right=0.95,wspace=0.12,hspace=0.25)
elif site_low == 'bojongsoang':
    fig = plt.figure(1,facecolor='w',figsize=(5.0,3.8))
    plt.subplots_adjust(top=0.97,bottom=0.09,left=0.034,right=0.96,wspace=0.12,hspace=0.25)
else:
    raise ValueError('Error in site >>> '+opts.site)
fig.clear()

prj = ccrs.UTM(zone=48,southern_hemisphere=True)
ax1 = plt.subplot(111,projection=prj)

for i,shp in enumerate(shapes):
    z = data[i]
    if not np.isnan(z):
        ax1.add_geometries(shp,prj,edgecolor='none',facecolor=mycmap((z-zmin)/zdif))
im1 = ax1.imshow(np.arange(4).reshape(2,2),extent=(-2,-1,-2,-1),vmin=zmin,vmax=zmax,cmap=mycmap)
ax12 = plt.colorbar(im1,ax=ax1,orientation='horizontal',shrink=1.00,pad=0.01).ax
if field_type == 'T':
    ax12.xaxis.set_major_locator(plt.FixedLocator(values))
    ax12.xaxis.set_major_formatter(plt.FixedFormatter(labels))
    ax12.xaxis.set_minor_locator(plt.FixedLocator(ticks))
    for l in ax12.xaxis.get_ticklabels():
        l.set_rotation(30)
if opts.zlabel is not None:
    ax12.set_xlabel(opts.zlabel)
    ax12.xaxis.set_label_coords(0.5,-1.8)
if opts.outline_fnam is not None:
    ax1.add_geometries(outline_shapes,prj,edgecolor='k',facecolor='none')

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
