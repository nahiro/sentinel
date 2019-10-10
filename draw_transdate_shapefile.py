#from shapely.geometry.polygon import Polygon
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

prj = ccrs.UTM(zone=48,southern_hemisphere=True)

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
for rec in records:
    t = rec.attributes['trans_date']
    if t > 1000 and t < tmin:
        tmin = t
    if t > tmax:
        tmax = t
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
                labels.append(d.strftime('%Y/%m/%d'))
            for day in [5,10,20,25]:
                d = datetime(y,m,day)
                ticks.append(date2num(d))

plt.interactive(True)
fig = plt.figure(1,facecolor='w',figsize=(7,6))
fig.clear()
plt.subplots_adjust(top=0.95,bottom=0.05,left=0.10,right=0.80)
ax1 = plt.subplot(111,projection=prj)
for shp,rec in zip(shapes,records):
    t = rec.attributes['trans_date']
    if t > 1000.0:
        ax1.add_geometries(shp,prj,edgecolor='none',facecolor=cm.jet((t-tmin)/tdif))
im = ax1.imshow(np.arange(4).reshape(2,2),extent=(-2,-1,-2,-1),vmin=tmin,vmax=tmax,cmap=cm.jet)
ax1.add_geometries(block_shp,prj,edgecolor='k',facecolor='none')
for bs,br in zip(block_shp,block_rec):
    for n in range(len(bs)):
        points = list(zip(*bs[n].exterior.coords.xy))
        p = Path(points)
        flags = p.contains_points(np.hstack((xg.flatten()[:,np.newaxis],yg.flatten()[:,np.newaxis]))).reshape(xg.shape)
        #xc = bs.boundary.centroid.x
        #yc = bs.boundary.centroid.y
        distance = np.empty(xg.shape,dtype=np.float64)
        for i,x_pt in enumerate(xp):
            for j,y_pt in enumerate(yp):
                if flags[j,i]:
                    distance[j,i] = bs[n].exterior.distance(shapely.geometry.Point(x_pt,y_pt))
                else:
                    distance[j,i] = -1.0e10
        indy,indx = np.unravel_index(np.argmax(distance),distance.shape)
        xc = xp[indx]
        yc = yp[indy]
        tc = br.attributes['Blok']
        ax1.text(xc,yc,tc,transform=prj,ha='center',va='center')
        #break
ax2 = plt.colorbar(im).ax
ax2.yaxis.set_major_locator(plt.FixedLocator(values))
ax2.yaxis.set_major_formatter(plt.FixedFormatter(labels))
ax2.yaxis.set_minor_locator(plt.FixedLocator(ticks))
ax1.set_xlim(xmin,xmax)
ax1.set_ylim(ymin,ymax)
plt.draw()
plt.savefig('transplanting_date.pdf')
