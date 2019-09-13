#!/usr/bin/env python
import sys
import numpy as np

sid,xpek,ypek,ns,ymax = np.loadtxt('select_peaks.dat',unpack=True)
sid = (sid+0.1).astype(np.int64)
ns = (ns+0.1).astype(np.int64)

for i in range(sid.max()+1):
    cnd = np.abs(sid-i) < 1.0e-4
    xpek_sid = xpek[cnd]
    ypek_sid = ypek[cnd]
    sind = []
    eind = []
    ind1 = np.arange(xpek_sid.size)
    if ind1.size > 0:
        sind.append(ind1[0])
        dind = np.diff(xpek_sid)
        ind2 = np.where(dind > 45.0)[0]
        if ind2.size > 0:
            for itmp in ind2:
                eind.append(ind1[itmp])
                sind.append(ind1[itmp+1])
        eind.append(ind1[-1])
    for si,ei in zip(sind,eind):
        x = xpek_sid[si:ei+1]
        y = ypek_sid[si:ei+1]
        j = np.argmax(y)
        sys.stdout.write('{:6d} {:15.8e} {:8.3f}\n'.format(i,x[j],y[j]))
    #break
