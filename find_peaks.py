#!/usr/bin/env python
import os
import sys
import re
from glob import glob
from datetime import datetime,timedelta
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.dates import date2num,num2date
from matplotlib.backends.backend_pdf import PdfPages
from optparse import OptionParser,IndentedHelpFormatter

# Default values
END = datetime.now().strftime('%Y%m%d')
PERIOD = 90 # day
XSGM = 5.0 # day
STHR = 0.5 # dB
NTHR = 1
BTHR = 1.5 # dB
YTHR = 1.0 # dB
DTHR = 0.05 # dB/day
VINT = 100
DATDIR = '.'
OUTNAM = 'find_peaks.dat'
FIGNAM = 'find_peaks.pdf'

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.set_usage('Usage: %prog collocated_geotiff_file [options]')
parser.add_option('-e','--end',default=END,help='End date of the analysis in the format YYYYMMDD (%default)')
parser.add_option('-E','--enum',default=None,type='int',help='End number (%default)')
parser.add_option('-p','--period',default=PERIOD,type='int',help='Observation period in day (%default)')
parser.add_option('-w','--xsgm',default=XSGM,type='float',help='Standard deviation of gaussian in day (%default)')
parser.add_option('-s','--sthr',default=STHR,type='float',help='Threshold of bg-subtracted signal in dB (%default)')
parser.add_option('-n','--nthr',default=NTHR,type='int',help='Minimum number of bg-subtracted signal above threshold (%default)')
parser.add_option('-b','--bthr',default=BTHR,type='float',help='Background level above average in dB (%default)')
parser.add_option('-y','--ythr',default=YTHR,type='float',help='Threshold of superposed gaussians in dB for peak search (%default)')
parser.add_option('--dthr',default=DTHR,type='float',help='Threshold of dy/dt in dB/day for local maximum search (%default)')
parser.add_option('-D','--datdir',default=DATDIR,help='Input data directory (%default)')
parser.add_option('-o','--outnam',default=OUTNAM,help='Output data name (%default)')
parser.add_option('-F','--fignam',default=FIGNAM,help='Output figure name for debug (%default)')
parser.add_option('--ax1_ymax',default=None,type='float',help='Axis1 Y max. (%default)')
parser.add_option('--ax1_ymin',default=None,type='float',help='Axis1 Y min. (%default)')
parser.add_option('--ax1_ystp',default=None,type='float',help='Axis1 Y step. (%default)')
parser.add_option('-S','--subpeak',default=False,action='store_true',help='Output sub peaks (%default)')
parser.add_option('-z','--npz',default=False,action='store_true',help='NPZ mode (%default)')
parser.add_option('--vint',default=VINT,type='int',help='Verbose output interval (%default)')
parser.add_option('-v','--verbose',default=False,action='store_true',help='Verbose mode (%default)')
parser.add_option('-d','--debug',default=False,action='store_true',help='Debug mode (%default)')
(opts,args) = parser.parse_args()
d1 = datetime.strptime(opts.end,'%Y%m%d')
d0 = d1-timedelta(days=opts.period)

if opts.npz:
    fs = sorted(glob(os.path.join(opts.datdir,'collocation_*.npz')))
else:
    fs = sorted(glob(os.path.join(opts.datdir,'collocation_*.dat')))

xmin = 1.0e10
xmax = -1.0e10
jdat = None
tdat = []
ydat = []
sdat = []
for i,f in enumerate(fs):
    if opts.npz:
        m = re.search('collocation_(\d+)_(\d+).npz',os.path.basename(f))
    else:
        m = re.search('collocation_(\d+)_(\d+).dat',os.path.basename(f))
    if not m:
        continue
    d0_tmp = datetime.strptime(m.group(1),'%y%m%d')
    d1_tmp = datetime.strptime(m.group(2),'%y%m%d')
    if d1_tmp < d0 or d1_tmp > d1:
        continue
    try:
        if opts.npz:
            data = np.load(f)
            j = data['j']
            tmin = data['tmin']
            vmin = data['vmin']
            bavg = data['bavg']
        else:
            j,ndat,tmin,vmin,fmin,tlft,vlft,flft,trgt,vrgt,frgt,dmin,dstd,tleg,fleg,treg,freg,sstd,scor,traw,vraw,fraw,draw,rstd,rcor,bavg = np.loadtxt(f,unpack=True)
    except Exception:
        continue
    if jdat is None:
        jdat = j
    else:
        if not np.all(j == jdat):
            raise ValueError('Error, different j')
    tdat.append(tmin)
    ydat.append(d1_tmp)
    sdat.append(bavg-vmin)
    x0_tmp = date2num(d0_tmp)
    x1_tmp = date2num(d1_tmp)
    if x0_tmp < xmin:
        xmin = x0_tmp
    if x1_tmp > xmax:
        xmax = x1_tmp
jdat = np.array(jdat)
tdat = np.array(tdat)
ydat = np.array(ydat)
sdat = np.array(sdat)

if opts.debug:
    plt.interactive(True)
    fig = plt.figure(1,facecolor='w',figsize=(6,3.5))
    pdf = PdfPages(opts.fignam)
    values = []
    labels = []
    ticks = []
    ds = (xmax-xmin)/365
    for y in range(num2date(xmin).year,num2date(xmax).year+1):
        if ds > 2.0:
            for m in range(1,13,3):
                d = datetime(y,m,1)
                values.append(date2num(d))
                labels.append(d.strftime('%Y-%m'))
            for m in range(1,13,1):
                d = datetime(y,m,1)
                ticks.append(date2num(d))
        elif ds > 1.0:
            for m in range(1,13,2):
                d = datetime(y,m,1)
                values.append(date2num(d))
                labels.append(d.strftime('%Y-%m'))
            for m in range(1,13,1):
                d = datetime(y,m,1)
                ticks.append(date2num(d))
        else:
            for m in range(1,13,1):
                d = datetime(y,m,1)
                values.append(date2num(d))
                labels.append(d.strftime('%Y-%m'))

xval = np.arange(xmin,xmax,0.1)
dx_1 = 1.0/np.gradient(xval)
with open(opts.outnam,'w') as fp:
    for sid in range(jdat.size):
        if opts.verbose and sid%opts.vint == 0:
            sys.stderr.write('{}\n'.format(sid))
        if opts.debug:
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
            if ind4.size > 0: # convex
                dind = np.diff(ind4)
                ind5 = np.where(dind > 1)[0]
                if ind5.size > 0: # multiple peak
                    if opts.subpeak:
                        sjs = []
                        ejs = []
                        sjs.append(ind4[0])
                        for jtmp in ind5:
                            ejs.append(ind4[jtmp])
                            sjs.append(ind4[jtmp+1])
                        ejs.append(ind4[-1])
                        xpks = []
                        ypks = []
                        for sj,ej in zip(sjs,ejs):
                            ind6 = ind3[sj:ej+1][np.nanargmax(yval[ind3][sj:ej+1])] # the first peak
                            xpks.append(xval[ind6])
                            ypks.append(yval[ind6])
                    else:
                        sj = ind4[0] # the first peak
                        ej = ind4[ind5[0]] # the first peak
                        ind6 = ind3[sj:ej+1][np.nanargmax(yval[ind3][sj:ej+1])] # the first peak
                        xpks = [xval[ind6]]
                        ypks = [np.nanmax(yval[ind3])] # may not equal to yval[ind6]
                else: # single peak
                    ind6 = ind0[ind3][np.nanargmax(yval[ind3])]
                    xpks = [xval[ind6]]
                    ypks = [yval[ind6]]
            else: # monotonic
                ind6 = ind0[ind3][np.nanargmax(yval[ind3])]
                xpks = [xval[ind6]]
                ypks = [yval[ind6]]
            for xpek,ypek in zip(xpks,ypks):
                fp.write('{:6d} {:15.8e} {:13.6e}\n'.format(sid,xpek,ypek))
            if opts.debug:
                ax1.hlines(opts.ythr,xval[si],xval[ei],zorder=10)
                if opts.subpeak and len(xpks) > 1:
                    for xpek,ypek in zip(xpks,ypks):
                        ax1.plot(xpek,ypek,'ro',zorder=10)
                else:
                    ax1.plot(xval[ind6],yval[ind6],'ro',zorder=10)
        if opts.debug:
            ax1.plot(xval,yval,'b-',zorder=1)
            ax1.set_ylabel('Signal value (dB)')
            ax1.set_xlim(xmin,xmax)
            if opts.ax1_ymin is not None:
                ax1.set_ylim(bottom=opts.ax1_ymin)
            if opts.ax1_ymax is not None:
                ax1.set_ylim(top=opts.ax1_ymax)
            if opts.ax1_ystp is not None:
                ax1.yaxis.set_major_locator(plt.MultipleLocator(opts.ax1_ystp))
            ax1.set_title('#{} ({} - {})'.format(sid,num2date(xmin).strftime('%Y-%m-%d'),num2date(xmax).strftime('%Y-%m-%d')))
            ax1.xaxis.set_major_locator(plt.FixedLocator(values))
            ax1.xaxis.set_major_formatter(plt.FixedFormatter(labels))
            ax1.xaxis.set_minor_locator(plt.FixedLocator(ticks))
            ax1.xaxis.set_tick_params(pad=7)
            ax1.yaxis.set_label_coords(-0.10,0.5)
            fig.autofmt_xdate()
            plt.draw()
            plt.savefig(pdf,format='pdf')
        if opts.enum is not None and sid >= opts.enum:
            break
if opts.debug:
    pdf.close()
