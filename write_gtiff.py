#!/usr/bin/env python
import os
import sys
import shutil
from datetime import datetime,timedelta
import numpy as np
try:
    import gdal
except Exception:
    from osgeo import gdal
try:
    import osr
except Exception:
    from osgeo import osr
from matplotlib.dates import date2num

outnam = os.path.join('.','transplanting_date.tif')
d0 = date2num(datetime(2017,3,1))
d1 = d0+180.0

xstp = 10.0
ystp = -10.0
xmin,xmax,ymin,ymax = (743800.0,756800.0,9236000.0,9251800.0)
xg,yg = np.meshgrid(np.arange(xmin,xmax+0.1*xstp,xstp),np.arange(ymax,ymin-0.1*ystp,ystp))
ngrd = xg.size
ny,nx = xg.shape

#sid,xpek,ypek,ns,ymax = np.loadtxt('thinout_peaks.dat',unpack=True)
sid,xpek,ypek = np.loadtxt('thinout_peaks.dat',unpack=True)
sid = (sid+0.1).astype(np.int64)
xpek_sid = [[] for i in range(ngrd)]
ypek_sid = [[] for i in range(ngrd)]
for i,x,y in zip(sid,xpek,ypek):
    if x >= d0 and x <= d1:
        xpek_sid[i].append(x)
        ypek_sid[i].append(y)
for i in range(ngrd):
    ndat = len(xpek_sid[i])
    if ndat > 1:
        indx = np.argmin(xpek_sid[i])
        xpek_sid[i] = [xpek_sid[i][indx]]
        ypek_sid[i] = [ypek_sid[i][indx]]
    elif ndat < 1:
        xpek_sid[i].append(np.nan)
        ypek_sid[i].append(np.nan)
xpek_sid = np.array(xpek_sid)
ypek_sid = np.array(ypek_sid)

drv = gdal.GetDriverByName('GTiff')
ds = drv.Create(outnam,nx,ny,2,gdal.GDT_Float32)
ds.SetGeoTransform((xmin-0.5*xstp,xstp,0.0,ymax-0.5*ystp,0.0,ystp))
srs = osr.SpatialReference()
srs.ImportFromEPSG(32748)
ds.SetProjection(srs.ExportToWkt())
band_data = [xpek_sid,ypek_sid]
band_name = ['xpek','ypek']
for i in range(2):
    band = ds.GetRasterBand(i+1)
    band.WriteArray(band_data[i].reshape(xg.shape))
    band.SetDescription(band_name[i])
band.SetNoDataValue(np.nan) # The TIFFTAG_GDAL_NODATA only support one value per dataset
ds.FlushCache()
ds = None # close dataset
