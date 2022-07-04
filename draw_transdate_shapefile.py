#from shapely.geometry.polygon import Polygon
import sys
from datetime import datetime
import numpy as np
import shapefile
import shapely
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.dates import date2num,num2date
from matplotlib.path import Path
from optparse import OptionParser,IndentedHelpFormatter

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('-s','--tmin',default=None,help='Min date in the format YYYYMMDD (%default)')
parser.add_option('-e','--tmax',default=None,help='Max date in the format YYYYMMDD (%default)')
parser.add_option('-p','--pmin',default=None,type='float',help='Min signal in dB (%default)')
parser.add_option('-P','--pmax',default=None,type='float',help='Max signal in dB (%default)')
parser.add_option('-t','--title',default=None,help='Figure title (%default)')
(opts,args) = parser.parse_args()

prj = ccrs.UTM(zone=48,southern_hemisphere=True)
utc_offset = 7.0/24.0 # day

fnam = '/home/naohiro/Work/SATREPS/Shapefile/cihea_testsite_mod/cihea_testsite_mod.shp'
block_shp = list(shpreader.Reader(fnam).geometries())
block_rec = list(shpreader.Reader(fnam).records())

fnam = 'transplanting_date.shp'
shapes = list(shpreader.Reader(fnam).geometries())
records = list(shpreader.Reader(fnam).records())
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

xp = np.arange(xmin,xmax,10.0)
yp = np.arange(ymin,ymax,10.0)
xg,yg = np.meshgrid(xp,yp)

tmin = 1.0e10
tmax = -1.0e10
pmin = 1.0e10
pmax = -1.0e10
for rec in records:
    t = rec.attributes['trans_date']+utc_offset
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
                labels.append(d.strftime('%Y/%m/%d'))
            for day in [5,10,20,25]:
                d = datetime(y,m,day)
                ticks.append(date2num(d))

plt.interactive(True)
fig = plt.figure(1,facecolor='w',figsize=(12,6))
fig.clear()
plt.subplots_adjust(top=0.93,bottom=0.03,left=0.02,right=0.95,wspace=0.4)

ax1 = plt.subplot(121,projection=prj)
ax2 = plt.subplot(122,projection=prj)
for shp,rec in zip(shapes,records):
    t = rec.attributes['trans_date']+utc_offset
    p = rec.attributes['peak_value']
    if t > 1000.0:
        ax1.add_geometries([shp],prj,edgecolor='none',facecolor=cm.jet((t-tmin)/tdif))
        ax2.add_geometries([shp],prj,edgecolor='none',facecolor=cm.jet((p-pmin)/pdif))
im1 = ax1.imshow(np.arange(4).reshape(2,2),extent=(-2,-1,-2,-1),vmin=tmin,vmax=tmax,cmap=cm.jet)
ax12 = plt.colorbar(im1,ax=ax1).ax
ax12.yaxis.set_major_locator(plt.FixedLocator(values))
ax12.yaxis.set_major_formatter(plt.FixedFormatter(labels))
ax12.yaxis.set_minor_locator(plt.FixedLocator(ticks))
ax12.set_ylabel('Transplanting date')
ax12.yaxis.set_label_coords(5.5,0.5)
ax1.add_geometries(block_shp,prj,edgecolor='k',facecolor='none')
im2 = ax2.imshow(np.arange(4).reshape(2,2),extent=(-2,-1,-2,-1),vmin=pmin,vmax=pmax,cmap=cm.jet)
ax22 = plt.colorbar(im2,ax=ax2).ax
ax22.minorticks_on()
#ax22.yaxis.set_major_locator(plt.FixedLocator(values))
#ax22.yaxis.set_major_formatter(plt.FixedFormatter(labels))
#ax22.yaxis.set_minor_locator(plt.FixedLocator(ticks))
ax22.set_ylabel('Signal (dB)')
ax22.yaxis.set_label_coords(3.5,0.5)
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
        #break
ax1.set_xlim(xmin,xmax)
ax1.set_ylim(ymin,ymax)
ax2.set_xlim(xmin,xmax)
ax2.set_ylim(ymin,ymax)
if opts.title is not None:
    plt.suptitle(opts.title)
plt.draw()
plt.savefig('transplanting_date.pdf')
