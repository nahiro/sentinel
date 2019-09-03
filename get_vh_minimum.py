import os
os.environ['_JAVA_OPTIONS'] = '-Xmx100240m'     # Save memory for JAVA
os.system('export _JAVA_OPTIONS=-Xmx100240m')   # Save memory for JAVA
import shutil
import sys
import re
from datetime import datetime
import numpy as np
from scipy.interpolate import splrep,splev
from csaps import UnivariateCubicSmoothingSpline
import gdal
import osr
import shapefile
from snappy import ProductIO
import matplotlib.pyplot as plt
from matplotlib.dates import date2num
from matplotlib.path import Path
from matplotlib.backends.backend_pdf import PdfPages

def transform_utm_to_wgs84(easting,northing,utm_zone):
    is_northern = (1 if northing.mean() > 0 else 0)
    utm_coordinate_system = osr.SpatialReference()
    utm_coordinate_system.SetWellKnownGeogCS('WGS84') # Set geographic coordinate system to handle lat/lon
    utm_coordinate_system.SetUTM(utm_zone,is_northern)
    wgs84_coordinate_system = utm_coordinate_system.CloneGeogCS() # Clone ONLY the geographic coordinate system
    utm_to_wgs84_geo_transform = osr.CoordinateTransformation(utm_coordinate_system,wgs84_coordinate_system) # create transform component
    xyz = np.array(utm_to_wgs84_geo_transform.TransformPoints(np.dstack((easting,northing)).reshape((-1,2)))).reshape(easting.shape[0],easting.shape[1],3)
    return xyz[:,:,0],xyz[:,:,1],xyz[:,:,2] # returns lon, lat, altitude

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

debug = True

if debug:
    plt.interactive(False)
    pdf = PdfPages('transplanting_date.pdf')
    fig = plt.figure(1,facecolor='w',figsize=(6,3.5))
    plt.subplots_adjust(top=0.85,bottom=0.20,left=0.15,right=0.90)
    plt.draw()

fnam = '../../190823/vh/collocation_VH.tif'
ds = gdal.Open(fnam)
data = ds.ReadAsArray()
trans = ds.GetGeoTransform()
indy,indx = np.indices(data[0].shape)
lon = trans[0]+indx*trans[1]+indy*trans[2]
lat = trans[3]+indx*trans[4]+indy*trans[5]
xp,yp,zp = transform_wgs84_to_utm(lon,lat)
ds = None # close dataset

vh_prod = ProductIO.readProduct(fnam)
vh_list = list(vh_prod.getBandNames())
nh = len(vh_list)
dates = []
for i in range(nh):
    m = re.search('_(\d+)$',vh_list[i])
    if not m:
        raise ValueError('Error in finding date >>> '+vh_list[i])
    dh = datetime.strptime(m.group(1),'%Y%m%d')
    dates.append(dh)
vh = vh_prod.getBand(vh_list[0])
x0 = 0
y0 = 0
d0 = datetime(2019,4,10)
d1 = datetime(2019,6,10)
w0 = vh.getRasterWidth()
h0 = vh.getRasterHeight()
x1 = w0
y1 = h0
w = x1-x0
h = y1-y0
dset = []
dtim = []
for i in range(nh):
    if dates[i] < d0:
        continue
    if dates[i] > d1:
        break
    data = np.zeros((h,w),dtype=np.float32)
    sys.stderr.write(vh_list[i]+'\n')
    vh = vh_prod.getBand(vh_list[i])
    wi = vh.getRasterWidth()
    hi = vh.getRasterHeight()
    if wi != w0 or hi != h0:
        raise ValueError('Error, different size >>> ({},{}), expected ({},{})'.format(wi,hi,w0,h0))
    vh.readPixels(x0,y0,w,h,data)
    dset.append(data)
    dtim.append(dates[i])
vh_prod.dispose()
dset = np.array(dset)
dtim = np.array(dtim)
ntim = date2num(dtim)
t0 = ntim.min()
t1 = ntim.max()
nt = ntim.size
xx = np.arange(t0,t1,0.01)
nx = xx.size
inds = np.arange(nx)

r = shapefile.Reader('../../190830/New_Test_Sites/New_Test_Sites')
fnam = 'transplanting_date.dat'
with open(fnam,'w') as fp:
    fp.write('# {:.5f} {:.5f} {:4d}\n'.format(t0,t1,nt))
    for i,shaperec in enumerate(r.iterShapeRecords()):
        if i%100 == 0:
            sys.stderr.write('{} {}\n'.format(i,len(r)))
        shp = shaperec.shape
        p = Path(shp.points)
        flags = p.contains_points(np.hstack((xp.flatten()[:,np.newaxis],yp.flatten()[:,np.newaxis]))).reshape(dset[0].shape)
        ndat = flags.sum()
        tmin = np.nan
        vmin = np.nan
        fmin = -1
        tlft = np.nan
        vlft = np.nan
        flft = -1
        trgt = np.nan
        vrgt = np.nan
        frgt = -1
        dmin = np.nan
        dstd = np.nan
        tleg = np.nan
        fleg = -1
        treg = np.nan
        freg = -1
        traw = np.nan
        vraw = np.nan
        fraw = -1
        draw = np.nan
        if ndat > 0:
            dtmp = dset[0][flags]
            cnd = ~np.isnan(dtmp)
            if cnd.sum() > 0:
                yi = dset[:,flags].reshape(dtim.size,-1).mean(axis=1)
                sp = UnivariateCubicSmoothingSpline(ntim,yi,smooth=0.05)
                yy = sp(xx)
                y1 = sp(ntim)
                indx_tmin = np.argmin(yy)
                tmin = xx[indx_tmin]
                vmin = yy[indx_tmin]
                fmin = (1 if indx_tmin == 0 else (2 if indx_tmin == nx-1 else 0))
                dy = np.gradient(yy)
                cnd = (dy >= 0.0) & (xx < tmin)
                if cnd.sum() > 0:
                    indx_tlft = inds[cnd][-1]
                    tlft = xx[indx_tlft]
                    vlft = yy[indx_tlft]
                    flft = 0
                else:
                    indx_tlft = 0
                    tlft = xx[indx_tlft]
                    vlft = yy[indx_tlft]
                    flft = 1
                cnd = (dy <= 0.0) & (xx > tmin)
                if cnd.sum() > 0:
                    indx_trgt = inds[cnd][0]
                    trgt = xx[indx_trgt]
                    vrgt = yy[indx_trgt]
                    frgt = 0
                else:
                    indx_trgt = nx-1
                    trgt = xx[indx_trgt]
                    vrgt = yy[indx_trgt]
                    frgt = 2
                dmin = vmin-splev(xx,splrep(ntim,yi,k=1))[indx_tmin]
                dstd = np.sqrt(np.square(y1-yi).sum()/y1.size)
                cnd0 = (yy >= vmin+dstd)
                cnd = cnd0 & (xx < tmin)
                if cnd.sum() > 0:
                    indx_tleg = inds[cnd][-1]
                    tleg = xx[indx_tleg]
                    vleg = yy[indx_tleg]
                    fleg = 0
                else:
                    indx_tleg = 0
                    tleg = xx[indx_tleg]
                    vleg = yy[indx_tleg]
                    fleg = 1
                cnd = cnd0 & (xx > tmin)
                if cnd.sum() > 0:
                    indx_treg = inds[cnd][0]
                    treg = xx[indx_treg]
                    vreg = yy[indx_treg]
                    freg = 0
                else:
                    indx_treg = nx-1
                    treg = xx[indx_treg]
                    vreg = yy[indx_treg]
                    freg = 2
                indx_traw = np.argmin(yi)
                traw = ntim[indx_traw]
                vraw = yi[indx_traw]
                fraw = (1 if indx_traw == 0 else (2 if indx_traw == nt-1 else 0))
                draw = sp([traw])[0]-vraw
        fp.write('{:6d} {:3d} '.format(i,ndat))
        fp.write('{:13.6e} {:13.6e} {:2d} '.format(tmin,vmin,fmin))
        fp.write('{:13.6e} {:13.6e} {:2d} '.format(tlft,vlft,flft))
        fp.write('{:13.6e} {:13.6e} {:2d} '.format(trgt,vrgt,frgt))
        fp.write('{:13.6e} {:13.6e} '.format(dmin,dstd))
        fp.write('{:13.6e} {:2d} '.format(tleg,fleg))
        fp.write('{:13.6e} {:2d} '.format(treg,freg))
        fp.write('{:13.6e} {:13.6e} {:2d} '.format(traw,vraw,fraw))
        fp.write('{:13.6e}\n'.format(draw))
        if debug and not np.isnan(tmin):
            fig.clear()
            ax1 = plt.subplot(111)
            ax1.plot(ntim,yi,'b-')
            ax1.plot(xx,yy,'r-')
            ax1.plot(tmin,vmin-dmin,'go')
            ax1.plot(traw,vraw+draw,'go')
            ax1.plot(tmin,vmin,'bo')
            ax1.plot(tlft,vlft,'mo')
            ax1.plot(trgt,vrgt,'mo')
            ax1.plot(tleg,vleg,'co')
            ax1.plot(treg,vreg,'co')
            ax1.set_ylim(-30.0,-10.0)
            ax1.set_title('{:d}'.format(i))
            plt.draw()
            plt.savefig(pdf,format='pdf')
        #break
if debug:
    pdf.close()
