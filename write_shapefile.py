#!/usr/bin/env python
import os
import sys
import shutil
from datetime import datetime
import numpy as np
import shapefile
from matplotlib.dates import date2num

inpnam = os.path.join('../../SATREPS','New_Test_Sites','New_Test_Sites')
outnam = os.path.join('.','transplanting_date')
d0 = date2num(datetime(2019,3,20))
d1 = date2num(datetime(2019,6,20))

#sid,xpek,ypek,ns,ymax = np.loadtxt('thinout_peaks.dat',unpack=True)
sid,xpek,ypek = np.loadtxt('thinout_peaks.dat',unpack=True)
sid = (sid+0.1).astype(np.int64)

r = shapefile.Reader(inpnam)
w = shapefile.Writer(outnam)
w.shapeType = shapefile.POLYGON
w.fields = r.fields[1:] # skip first deletion field
w.field('trans_date','F',13,6)
w.field('peak_value','F',13,6)
for i,shaperec in enumerate(r.iterShapeRecords()):
    rec = shaperec.record
    shp = shaperec.shape
    cnd = (np.abs(sid-i) < 1.0e-4) & (xpek > d0) & (xpek < d1)
    ndat = cnd.sum()
    if ndat > 1:
        raise ValueError('Error, i={}, ndat={}'.format(i,ndat))
    elif ndat < 1:
        rec.append(0.0)
        rec.append(0.0)
    else:
        rec.append(xpek[cnd][0])
        rec.append(ypek[cnd][0])
    w.shape(shp)
    w.record(*rec)
w.close()
shutil.copy2(inpnam+'.prj',outnam+'.prj')
