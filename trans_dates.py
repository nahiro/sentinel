#!/usr/bin/env python
import os
import sys
import glob
from datetime import datetime,timedelta
import numpy as np
from subprocess import call
from vh_2 import *

bindir = '.'
datdir = '../190906/vh'
polarization = 'VH'
period = 60 # days

# Create a filelist
dates = []
files = []
for f in glob.glob(os.path.join(datdir,'*_dB.dim')):
    d = os.path.basename(f)[0:8]
    dates.append(datetime.strptime(d,'%Y%m%d'))
    files.append(f)
dates = np.array(dates)
indx = np.argsort(dates)
datelist = dates[indx]
filelist = np.array(files)[indx]

data_prev = []
for d1 in datelist:
    d0 = d1-timedelta(days=period)
    if d0 < datelist[0]:
        continue
    data = []  # Empty list
    for i in range(len(filelist)):
        date = datelist[i]
        fnam = filelist[i]
        if date < d0 or date > d1:
            continue
        if os.path.exists(fnam):
            data.append(str(fnam))
        #print(d0,d1,date)
    if data == data_prev:
        continue
    gnam = 'collocation_{:%y%m%d}_{:%y%m%d}.tif'.format(d0,d1)
    collocate_all(data,write_product=True,write_fnam=gnam)
    command = os.path.join(bindir,'get_vh_minimum.py')
    command += ' {}'.format(gnam)
    command += ' -e {}'.format(d1.strftime('%Y%m%d'))
    command += ' -s ../../SATREPS/New_Test_Sites/New_Test_Sites.shp'
    command += ' -o collocation_{:%y%m%d}_{:%y%m%d}.dat'.format(d0,d1)
    sys.stderr.write(command+'\n')
    call(command,shell=True)
    data_prev = data
    #break
