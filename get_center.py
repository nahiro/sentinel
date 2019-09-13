#!/usr/bin/env python
import os
import sys
import numpy as np
import shapefile

shpnam = os.path.join('../../SATREPS','New_Test_Sites','New_Test_Sites.shp')

r = shapefile.Reader(shpnam)
for i,shaperec in enumerate(r.iterShapeRecords()):
    rec = shaperec.record
    shp = shaperec.shape
    pp = np.array(shp.points)
    xc = pp[:,0].mean()
    yc = pp[:,1].mean()
    sys.stdout.write('{:6d} {:15.8e} {:15.8e} {:6d} {:15.8e} {:15.8e}\n'.format(i,xc,yc,len(pp),rec[2],rec[3]))
