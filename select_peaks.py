#!/usr/bin/env python
import sys
import numpy as np
from optparse import OptionParser,IndentedHelpFormatter

# Default values
NTHR = 5
YTHR = 10.0 # dB

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('--ngrd',default=None,type='int',help='Number of image grid for NPZ mode (%default)')
parser.add_option('-n','--nthr',default=NTHR,type='int',help='Minimum number of nearby signal above threshold (%default)')
parser.add_option('-y','--ythr',default=YTHR,type='float',help='Threshold of superposed gaussians in dB for peak selection (%default)')
parser.add_option('-z','--npz',default=False,action='store_true',help='NPZ mode (%default)')
(opts,args) = parser.parse_args()

if opts.npz:
    # read nearby indices
    data = np.load('find_nearest.npz')
    sid_0 = data['sid_0']
    sid_1 = data['sid_1']
    sid_2 = data['sid_2']
    sid_3 = data['sid_3']
    sid_4 = data['sid_4']
    sid_5 = data['sid_5']
    sid_6 = data['sid_6']
    sid_7 = data['sid_7']
    sid_8 = data['sid_8']
    sid_9 = data['sid_9']
    sid_a = data['sid_a']
    sid_b = data['sid_b']
    sid_c = data['sid_c']
    ntmp = sid_0.max()+1
    if opts.ngrd is None:
        opts.ngrd = ntmp
    else:
        if opts.ngrd != ntmp:
            sys.stderr.write('Warning, opts.ngrd={}, ntmp={}\n'.format(opts.ngrd,ntmp))
    # read peaks
    data = np.load('find_peaks.npz')
    sid = data['sid']
    xpek = data['xpek']
    ypek = data['ypek']
    xpek_sid = [[] for i in range(opts.ngrd)]
    ypek_sid = [[] for i in range(opts.ngrd)]
    for i,x,y in zip(sid,xpek,ypek):
        xpek_sid[i].append(x)
        ypek_sid[i].append(y)
    for i in range(opts.ngrd):
        xpek_sid[i] = np.array(xpek_sid[i])
        ypek_sid[i] = np.array(ypek_sid[i])
else:
    # read nearby indices
    sid_0,sid_1,leng_1,sid_2,leng_2,sid_3,leng_3,sid_4,leng_4,sid_5,leng_5,sid_6,leng_6,sid_7,leng_7,sid_8,leng_8,sid_9,leng_9,sid_a,leng_a,sid_b,leng_b,sid_c,leng_c = np.loadtxt('find_nearest.dat',unpack=True)
    sid_0 = (sid_0+0.1).astype(np.int64)
    sid_1 = (sid_1+0.1).astype(np.int64)
    sid_2 = (sid_2+0.1).astype(np.int64)
    sid_3 = (sid_3+0.1).astype(np.int64)
    sid_4 = (sid_4+0.1).astype(np.int64)
    sid_5 = (sid_5+0.1).astype(np.int64)
    sid_6 = (sid_6+0.1).astype(np.int64)
    sid_7 = (sid_7+0.1).astype(np.int64)
    sid_8 = (sid_8+0.1).astype(np.int64)
    sid_9 = (sid_9+0.1).astype(np.int64)
    sid_a = (sid_a+0.1).astype(np.int64)
    sid_b = (sid_b+0.1).astype(np.int64)
    sid_c = (sid_c+0.1).astype(np.int64)
    # read peaks
    sid,xpek,ypek = np.loadtxt('find_peaks.dat',unpack=True)
    sid = (sid+0.1).astype(np.int64)
    xpek_sid = []
    ypek_sid = []
    for i in sid_0:
        cnd = np.abs(sid-i) < 1.0e-4
        xpek_sid.append(xpek[cnd])
        ypek_sid.append(ypek[cnd])

for i in sid_0:
    indx = np.array([sid_1[i],sid_2[i],sid_3[i],sid_4[i],sid_5[i],sid_6[i],sid_7[i],sid_8[i],sid_9[i],sid_a[i],sid_b[i],sid_c[i]])
    for x,y in zip(xpek_sid[i],ypek_sid[i]):
        ys = []
        for j in indx:
            if xpek_sid[j].size < 1:
                continue
            dx = np.abs(xpek_sid[j]-x)
            k = np.argmin(dx)
            if dx[k] < 10.0:
                ys.append(ypek_sid[j][k])
        ys = np.array(ys)
        ns = ys.size
        if ns > 0:
            ymax = ys.max()
        else:
            ymax = 0.0
        if y >= 10.0:
            flag = True
        elif ns >= 5 and ymax >= 10.0:
            flag = True
        else:
            flag = False
        if flag:
            sys.stdout.write('{:6d} {:15.8e} {:8.3f} {:6d} {:8.3f}\n'.format(i,x,y,ns,ymax))
