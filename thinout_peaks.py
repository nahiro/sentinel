#!/usr/bin/env python
import sys
import numpy as np
from optparse import OptionParser,IndentedHelpFormatter

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('--ngrd',default=None,type='int',help='Number of image grid for NPZ mode (%default)')
parser.add_option('-z','--npz',default=False,action='store_true',help='NPZ mode (%default)')
(opts,args) = parser.parse_args()

if opts.npz:
    data = np.load('select_peaks.npz')
    sid = data['sid']
    xpek = data['xpek']
    ypek = data['ypek']
else:
    sid,xpek,ypek,ns,ymax = np.loadtxt('select_peaks.dat',unpack=True)
    sid = (sid+0.1).astype(np.int64)
ntmp = sid.max()+1
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
        x = xpek_sid[i][si:ei+1]
        y = ypek_sid[i][si:ei+1]
        j = np.argmax(y)
        sys.stdout.write('{:8d} {:15.8e} {:8.3f}\n'.format(i,x[j],y[j]))
    #break
