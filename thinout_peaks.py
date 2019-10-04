#!/usr/bin/env python
import sys
import numpy as np
from optparse import OptionParser,IndentedHelpFormatter

# Default values
XSGM = 5.0 # day
LSGM = 50.0 # m

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('--ngrd',default=None,type='int',help='Number of image grid for NPZ mode (%default)')
parser.add_option('-w','--xsgm',default=XSGM,type='float',help='Standard deviation of gaussian in day (%default)')
parser.add_option('-W','--lsgm',default=LSGM,type='float',help='Standard deviation of gaussian in m (%default)')
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
    leng_1 = data['leng_1']
    leng_2 = data['leng_2']
    leng_3 = data['leng_3']
    leng_4 = data['leng_4']
    leng_5 = data['leng_5']
    leng_6 = data['leng_6']
    leng_7 = data['leng_7']
    leng_8 = data['leng_8']
    leng_9 = data['leng_9']
    leng_a = data['leng_a']
    leng_b = data['leng_b']
    leng_c = data['leng_c']
    # read peaks
    data = np.load('select_peaks.npz')
    sid = data['sid'] # do NOT assume that peak data include all indices.
    xpek = data['xpek']
    ypek = data['ypek']
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
    sid,xpek,ypek,ns,ymax = np.loadtxt('select_peaks.dat',unpack=True)
    sid = (sid+0.1).astype(np.int64)
ntmp = sid_0.max()+1
if opts.ngrd is None:
    opts.ngrd = ntmp
elif opts.ngrd != ntmp:
    sys.stderr.write('Warning, opts.ngrd={}, ntmp={}\n'.format(opts.ngrd,ntmp))
xpek_sid = [[] for i in range(opts.ngrd)]
ypek_sid = [[] for i in range(opts.ngrd)]
for i,x,y in zip(sid,xpek,ypek):
    xpek_sid[i].append(x)
    ypek_sid[i].append(y)
for i in range(opts.ngrd):
    xpek_sid[i] = np.array(xpek_sid[i])
    ypek_sid[i] = np.array(ypek_sid[i])

for i in range(opts.ngrd):
    if i != sid_0[i]:
        raise ValueError('Error, i={}, sid_0={}'.format(i,sid_0[i]))
    sind = []
    eind = []
    ind1 = np.arange(xpek_sid[i].size)
    if ind1.size > 0:
        sind.append(ind1[0])
        dind = np.diff(xpek_sid[i])
        ind2 = np.where(dind > 45.0)[0]
        if ind2.size > 0:
            for itmp in ind2:
                eind.append(ind1[itmp])
                sind.append(ind1[itmp+1])
        eind.append(ind1[-1])
    for si,ei in zip(sind,eind):
        ysum = []
        for xi,yi in zip(xpek_sid[i][si:ei+1],ypek_sid[i][si:ei+1]):
            ys = yi
            for j,leng in zip([sid_1[i],sid_2[i],sid_3[i],sid_4[i],sid_5[i],sid_6[i],sid_7[i],sid_8[i],sid_9[i],sid_a[i],sid_b[i],sid_c[i]],
                              [leng_1[i],leng_2[i],leng_3[i],leng_4[i],leng_5[i],leng_6[i],leng_7[i],leng_8[i],leng_9[i],leng_a[i],leng_b[i],leng_c[i]]):
                for xj,yj in zip(xpek_sid[j],ypek_sid[j]):
                    if np.abs(xj-xi) < opts.xsgm*10.0:
                        ys += yj*np.exp(-0.5*np.square((xj-xi)/opts.xsgm))*np.exp(-0.5*np.square(leng/opts.lsgm))
            ysum.append(ys)
        ysum = np.array(ysum)
        x = xpek_sid[i][si:ei+1]
        y = ypek_sid[i][si:ei+1]
        j = np.argmax(ysum)
        sys.stdout.write('{:8d} {:15.8e} {:8.3f} {:8.3f}\n'.format(i,x[j],y[j],ysum[j]))
    #break
