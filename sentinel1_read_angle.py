#!/usr/bin/env python
import os
import sys
import re
import numpy as np
from subprocess import check_output

srcdir = '/home/naohiro/Work/Sentinel-1/Cihea/sigma0_speckle'
#srcdir = '/home/naohiro/Work/Sentinel-1/Bojongsoang/sigma0_speckle'

angs = []

for f in sorted(os.listdir(srcdir)):
    m = re.search('('+'\d'*8+')_resample.tif',f)
    if not m:
        continue
    fnam = os.path.join(srcdir,f)
    dstr = m.group(1)
    #print(dstr)
    command = 'gdalinfo'
    command += ' '+fnam
    out = check_output(command,shell=True)
    iangle = None
    direction = None
    for line in out.decode().splitlines():
        m = re.search('iangle\s*=\s*(\S+)',line)
        if m:
            iangle = float(m.group(1))
        m = re.search('PASS\s*=\s*(\S+)',line)
        if m:
            direction = m.group(1)
    if iangle is None:
        raise ValueError('Error in finding iangle')
    if direction is None:
        raise ValueError('Error in finding direction')
    if direction.lower() == 'descending':
        iangle *= -1
    elif direction.lower() == 'ascending':
        pass
    else:
        raise ValueError('Error, direction={}'.format(direction))
    if len(angs) < 1:
        angs.append([iangle])
    else:
        angs_diff = []
        for i in range(len(angs)):
            angs_mean = np.nanmean(angs[i])
            angs_diff.append(np.abs(iangle-angs_mean))
        indx = np.argmin(angs_diff)
        if angs_diff[indx] < 0.1:
            angs[indx].append(iangle)
        else:
            angs.append([iangle])
    sys.stdout.write('{} {:10.6f}\n'.format(dstr,iangle))

angs_mean = []
angs_std = []
angs_min = []
angs_max = []
for i in range(len(angs)):
    angs_mean.append(np.nanmean(angs[i]))
    angs_std.append(np.nanstd(angs[i]))
    angs_min.append(np.nanmin(angs[i]))
    angs_max.append(np.nanmax(angs[i]))
indx = np.argsort(np.abs(angs_mean))
for i in range(len(angs)):
    sys.stdout.write('{:2d} {:10.6f} {:10.6f} {:10.6f} {:10.6f}\n'.format(i,angs_mean[i],angs_std[i],angs_min[i]-angs_mean[i],angs_max[i]-angs_mean[i]))
