#!/usr/bin/env python
import numpy as np
from matplotlib.dates import num2date
from csaps import csaps
from scipy.interpolate import splrep,splev

p_smooth = 0.005
vthr1 = 0.06
vthr2 = 0.1

nv_cor = np.load('atcor_nv_corrected.npy')
ntim = np.load('../ntim.npy')
dtim = num2date(ntim)
xorg = ntim.copy()
nobject = nv_cor.shape[1]

cloud_flag = []
for iobj in range(nobject):
    object_id = iobj+1
    yorg = nv_cor[:,iobj]
    ysmo = csaps(xorg,yorg,xorg,smooth=p_smooth)
    inds = []
    for itim in range(xorg.size):
        if yorg[itim] < ysmo[itim]-vthr1:
            inds.append(itim)
    inds = np.array(inds)
    if inds.size > 0:
        xinp = np.delete(xorg,inds)
        yinp = np.delete(yorg,inds)
    else:
        xinp = xorg.copy()
        yinp = yorg.copy()
    yref = csaps(xinp,yinp,xorg,smooth=p_smooth)
    flag = []
    for itim in range(ntim.size):
        if (np.abs(yorg[itim]-yref[itim]) > vthr2):
            flag.append(True)
        else:
            flag.append(False)
    flag = np.array(flag)
    cloud_flag.append(flag)
cloud_flag = np.array(cloud_flag)
np.save('cloud_flag.npy',cloud_flag)
