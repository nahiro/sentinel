#!/usr/bin/env python
import os
import sys
import re
from glob import glob
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.dates import date2num,num2date
from matplotlib.backends.backend_pdf import PdfPages
from optparse import OptionParser,IndentedHelpFormatter

# Default values
XSGM = 5.0 # day
STHR = 1.0
NTHR = 1
BTHR = 1.5 # dB
YTHR = 2.0
DTHR = 0.05

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.set_usage('Usage: %prog collocated_geotiff_file [options]')
parser.add_option('-w','--xsgm',default=XSGM,type='float',help='Standard deviation of gaussian in day (%default)')
parser.add_option('-s','--sthr',default=STHR,type='float',help='Threshold of bg-subtracted signal in dB (%default)')
parser.add_option('-n','--nthr',default=NTHR,type='int',help='Minimum number of bg-subtracted signal above threshold (%default)')
parser.add_option('-b','--bthr',default=BTHR,type='float',help='Background level above average in dB (%default)')
parser.add_option('-y','--ythr',default=YTHR,type='float',help='Threshold of superposed gaussians for peak search (%default)')
parser.add_option('-d','--dthr',default=DTHR,type='float',help='Threshold of dy/dt for local maximum search (%default)')
(opts,args) = parser.parse_args()

datdir = '../190906'
fs = sorted(glob(os.path.join(datdir,'collocation_*.dat')))

xmin = 1.0e10
xmax = -1.0e10
jdat = None
tdat = []
ydat = []
sdat = []
for i,f in enumerate(fs):
    m = re.search('collocation_(\d+)_(\d+).dat',os.path.basename(f))
    if not m:
        continue
    d0 = date2num(datetime.strptime(m.group(1),'%y%m%d'))
    d1 = date2num(datetime.strptime(m.group(2),'%y%m%d'))
    try:
        j,ndat,tmin,vmin,fmin,tlft,vlft,flft,trgt,vrgt,frgt,dmin,dstd,tleg,fleg,treg,freg,sstd,scor,traw,vraw,fraw,draw,rstd,rcor,bavg = np.loadtxt(f,unpack=True)
    except Exception:
        continue
    if jdat is None:
        jdat = j
    else:
        if not np.all(j == jdat):
            raise ValueError('Error, different j')
    tdat.append(tmin)
    ydat.append(d1)
    sdat.append(bavg-vmin)
    if d0 < xmin:
        xmin = d0
    if d1 > xmax:
        xmax = d1
jdat = np.array(jdat)
tdat = np.array(tdat)
ydat = np.array(ydat)
sdat = np.array(sdat)

plt.interactive(True)
fig = plt.figure(1,facecolor='w',figsize=(6,3.5))
pdf = PdfPages('find_peaks.pdf')

values = []
labels = []
ticks = []
for y in [2017,2018,2019]:
    for m in range(1,13,3):
        d = datetime(y,m,1)
        values.append(date2num(d))
        labels.append(d.strftime('%Y-%m'))
    for m in range(1,13,1):
        d = datetime(y,m,1)
        ticks.append(date2num(d))

xval = np.arange(xmin,xmax,0.1)
dx_1 = 1.0/np.gradient(xval)
for sid in range(jdat.size):
    if sid%100 == 0:
        sys.stderr.write('{}\n'.format(sid))
    fig.clear()
    plt.subplots_adjust(top=0.85,bottom=0.20,left=0.15,right=0.90)
    ax1 = plt.subplot(111)
    ax1.minorticks_on()
    ax1.grid(True)
    cnd = np.abs(jdat-sid) < 1.0e-4
    tt = tdat[:,cnd].reshape(ydat.shape)
    ss = sdat[:,cnd].reshape(ydat.shape)
    indx = np.argsort(ss)[0:int(ss.size*0.6)]
    bb = ss[indx]
    bavg = np.nanmean(bb)
    bstd = np.nanstd(bb)
    blev = bavg+opts.bthr
    ysig = ss-blev
    cnd2 = (ysig > 0.0)
    yval = (ysig[cnd2].reshape(-1,1)*np.exp(-0.5*np.square((xval.reshape(1,-1)-tt[cnd2].reshape(-1,1))/opts.xsgm))).sum(axis=0)
    y_1d = np.gradient(yval)*dx_1
    y_2d = np.gradient(y_1d)*dx_1
    ind0 = np.arange(yval.size)
    sind = []
    eind = []
    ind1 = np.where(yval > opts.ythr)[0]
    if ind1.size > 0:
        sind.append(ind1[0])
        dind = np.diff(ind1)
        ind2 = np.where(dind > 1)[0]
        if ind2.size > 0:
            for itmp in ind2:
                eind.append(ind1[itmp])
                sind.append(ind1[itmp+1])
        eind.append(ind1[-1])

    for si,ei in zip(sind,eind):
        cnd3 = (tt >= xval[si]) & (tt <= xval[ei]) & (ysig > opts.sthr)
        if cnd3.sum() < opts.nthr:
            continue
        ind3 = ind0[si:ei+1]
        ind4 = np.where((np.abs(y_1d[ind3]) < opts.dthr) & (y_2d[ind3] < 0.0))[0]
        if ind4.size > 0:
            dind = np.diff(ind4)
            ind5 = np.where(dind > 1)[0]
            if ind5.size > 0:
                sj = ind4[0]
                ej = ind4[ind5[0]]
                ind6 = ind3[sj:ej+1][np.argmax(yval[ind3][sj:ej+1])]
            else:
                ind6 = ind0[ind3][np.argmax(yval[ind3])]
        else:
            ind6 = ind0[ind3][np.argmax(yval[ind3])]
        #ax1.plot(xval[ind3],yval[ind3],'c-',zorder=10)
        ax1.hlines(opts.ythr,xval[si],xval[ei],zorder=10)
        ax1.plot(xval[ind6],yval[ind6],'ro',zorder=10)
    ax1.plot(xval,yval,'b-',zorder=1)
    ax1.set_ylabel('Signal value (dB)')
    ax1.set_xlim(xmin,xmax)
    ax1.set_ylim(0.0,50.0)
    ax1.set_title('#{} ({} - {})'.format(sid,num2date(xmin).strftime('%Y-%m-%d'),num2date(xmax).strftime('%Y-%m-%d')))
    ax1.xaxis.set_major_locator(plt.FixedLocator(values))
    ax1.xaxis.set_major_formatter(plt.FixedFormatter(labels))
    ax1.xaxis.set_minor_locator(plt.FixedLocator(ticks))
    ax1.xaxis.set_tick_params(pad=7)
    ax1.yaxis.set_label_coords(-0.10,0.5)
    fig.autofmt_xdate()
    plt.draw()
    plt.savefig(pdf,format='pdf')
    #break
pdf.close()
