#!/usr/bin/env python
import sys
import numpy as np

sid,xc,yc,ndat,leng,area = np.loadtxt('get_center.dat',unpack=True)
sid = (sid+0.1).astype(np.int64)
ndat = (ndat+0.1).astype(np.int64)

for i in range(sid.size):
    l2 = np.square(xc-xc[i])+np.square(yc-yc[i])
    indx = np.argsort(l2)
    sys.stdout.write('{:6d}'.format(i))
    for j in indx[1:11]:
        sys.stdout.write(' {:6d} {:4.1f}'.format(sid[j],np.sqrt(l2[j])))
    sys.stdout.write('\n')
        
