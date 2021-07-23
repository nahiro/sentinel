#!/usr/bin/env python
import os
import sys
import re
from glob import glob
from datetime import datetime
import numpy as np
from matplotlib.dates import date2num
from csaps import csaps
from scipy.interpolate import splrep,splev
from optparse import OptionParser,IndentedHelpFormatter

# Default values
TMIN = '20190315'
TMAX = '20190615'
SMOOTH = 0.005
VTHR1 = 0.06
VTHR2 = 0.1
SEARCH_KEY = 'correct'
DATDIR = os.curdir
OUTDIR = os.curdir

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('-s','--tmin',default=TMIN,help='Min date of output data in the format YYYYMMDD (%default)')
parser.add_option('-e','--tmax',default=TMAX,help='Max date of output data in the format YYYYMMDD (%default)')
parser.add_option('--data_tmin',default=None,help='Min date of input data in the format YYYYMMDD (%default)')
parser.add_option('--data_tmax',default=None,help='Max date of input data in the format YYYYMMDD (%default)')
parser.add_option('-S','--smooth',default=SMOOTH,type='float',help='Smoothing factor from 0 to 1 (%default)')
parser.add_option('-v','--vthr1',default=VTHR1,type='float',help='Threshold 1 (%default)')
parser.add_option('-V','--vthr2',default=VTHR2,type='float',help='Threshold 2 (%default)')
parser.add_option('--search_key',default=SEARCH_KEY,help='Search key for input data (%default)')
parser.add_option('-D','--datdir',default=DATDIR,help='Input data directory (%default)')
parser.add_option('-O','--outdir',default=OUTDIR,help='Output data directory (%default)')
(opts,args) = parser.parse_args()

dmin = datetime.strptime(opts.tmin,'%Y%m%d')
dmax = datetime.strptime(opts.tmax,'%Y%m%d')
data_dmin = datetime.strptime(opts.data_tmin,'%Y%m%d')
data_dmax = datetime.strptime(opts.data_tmax,'%Y%m%d')

dtim = []
data = []
nobject = None
fs = sorted(glob(os.path.join(opts.datdir,'*'+'[0-9]'*8+'*.npz')))
for fnam in fs:
    f = os.path.basename(fnam)
    if opts.search_key is not None and not re.search(opts.search_key,f):
        continue
    m = re.search('\D('+'\d'*8+')\D',f)
    if not m:
        m = re.search('^('+'\d'*8+')\D',f)
        if not m:
            raise ValueError('Error in finding date >>> '+f)
    dstr = m.group(1)
    d = datetime.strptime(dstr,'%Y%m%d')
    if d < data_dmin or d > data_dmax:
        continue
    #sys.stderr.write(f+' '+dstr+'\n')
    dtmp = np.load(fnam)['data_cor']
    if nobject is None:
        nobject = dtmp.size
    elif dtmp.size != nobject:
        raise ValueError('Error, dtmp.size={}, nobject={}'.format(dtmp.size,nobject))
    dtim.append(d)
    data.append(dtmp)
dtim = np.array(dtim)
data = np.array(data)
ntim = date2num(dtim)
xorg = ntim.copy()

cloud_flag = []
for iobj in range(nobject):
    object_id = iobj+1
    yorg = data[:,iobj]
    ysmo = csaps(xorg,yorg,xorg,smooth=opts.smooth)
    inds = []
    for itim in range(xorg.size):
        if yorg[itim] < ysmo[itim]-opts.vthr1:
            inds.append(itim)
    inds = np.array(inds)
    if inds.size > 0:
        xinp = np.delete(xorg,inds)
        yinp = np.delete(yorg,inds)
    else:
        xinp = xorg.copy()
        yinp = yorg.copy()
    yref = csaps(xinp,yinp,xorg,smooth=opts.smooth)
    flag = []
    for itim in range(ntim.size):
        if (np.abs(yorg[itim]-yref[itim]) > opts.vthr2):
            flag.append(True)
        else:
            flag.append(False)
    flag = np.array(flag)
    cloud_flag.append(flag)
cloud_flag = np.array(cloud_flag).swapaxes(0,1)

for i in range(dtim.size):
    d = dtim[i]
    if d < dmin or d > dmax:
        continue
    np.save(os.path.join(opts.outdir,'{:%Y%m%d}_cflag.npy'.format(d)),cloud_flag[i])
