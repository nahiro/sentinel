#!/usr/bin/env python
import os
import sys
import re
import warnings
from glob import glob
from datetime import datetime
import gdal
import osr
import numpy as np
import cartopy.io.shapereader as shpreader
from csaps import csaps
from statsmodels.stats.weightstats import DescrStatsW
from scipy.signal import find_peaks
from optparse import OptionParser,IndentedHelpFormatter

# Default values
TMIN = '20190501'
TMAX = '20190915'
TMRG = 20.0
DATA_TMIN = '20190201'
DATA_TMAX = '20200201'
TMGN = 15.0 # day
FACT1 = 0.035
TSTP = 0.1 # day
TSTR = -20.0 # day
TEND = 20.0 # day
SMOOTH = 0.002
OFFSET_V = 4.8
OFFSET_FACT = 1.25
NDVI_DISTANCE = 10
NDVI_PROMINENCE = 0.002
NDVI2_DISTANCE = 1.0
NDVI2_PROMINENCE = 0.001
DATDIR = os.curdir
ndvi_min = -0.4
ndvi_max = 1.1
NDVI_MEAN = 0.355
ndvi_std = 5.929824e-02
ndvi_inc_mean = 4.564639e-01
ndvi_inc_std = 9.204571e-02
ndvi_dd_mean = 1.595052e-01
ndvi_dd_std = 4.939033e-02

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('-s','--tmin',default=TMIN,help='Min date of transplanting in the format YYYYMMDD (%default)')
parser.add_option('-e','--tmax',default=TMAX,help='Max date of transplanting in the format YYYYMMDD (%default)')
parser.add_option('--tmrg',default=TMRG,type='float',help='Marge width in day (%default)')
parser.add_option('--data_tmin',default=DATA_TMIN,help='Min date of input data in the format YYYYMMDD (%default)')
parser.add_option('--data_tmax',default=DATA_TMAX,help='Max date of input data in the format YYYYMMDD (%default)')
parser.add_option('--tmgn',default=TMGN,type='float',help='Margin of input data in day (%default)')
parser.add_option('--fact1',default=FACT1,type='float',help='Factor for S1 (%default)')
parser.add_option('--ndvi_mean',default=NDVI_MEAN,type='float',help='NDVI mean for S1 (%default)')
parser.add_option('--tstp',default=TSTP,type='float',help='Precision of transplanting date in day (%default)')
parser.add_option('--tstr',default=TSTR,type='float',help='Start day of transplanting period seen from the min. peak (%default)')
parser.add_option('--tend',default=TEND,type='float',help='End day of transplanting period seen from the min. peak (%default)')
parser.add_option('-S','--smooth',default=SMOOTH,type='float',help='Smoothing factor from 0 to 1 (%default)')
parser.add_option('--ndvi_distance',default=NDVI_DISTANCE,type='int',help='Minimum peak distance in day for Sentinel-2 (%default)')
parser.add_option('--ndvi_prominence',default=NDVI_PROMINENCE,type='float',help='Minimum prominence of NDVI for Sentinel-2 (%default)')
parser.add_option('--ndvi2_distance',default=NDVI2_DISTANCE,type='int',help='Minimum peak distance in day for NDVI'' (%default)')
parser.add_option('--ndvi2_prominence',default=NDVI2_PROMINENCE,type='float',help='Minimum prominence of NDVI'' for Sentinel-2 (%default)')
parser.add_option('-D','--datdir',default=DATDIR,help='Input data directory (%default)')
parser.add_option('--cflag_fnam',default=None,help='Cloud flag file name (%default)')
parser.add_option('-o','--output_fnam',default=None,help='Output NPZ name (%default)')
parser.add_option('-F','--fig_fnam',default=None,help='Output figure name for debug (%default)')
parser.add_option('--debug',default=False,action='store_true',help='Debug mode (%default)')
parser.add_option('-b','--batch',default=False,action='store_true',help='Batch mode (%default)')
(opts,args) = parser.parse_args()
warnings.simplefilter('ignore')

if opts.batch:
    import matplotlib
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.dates import date2num,num2date
from matplotlib.backends.backend_pdf import PdfPages

nmin = date2num(datetime.strptime(opts.tmin,'%Y%m%d'))
nmax = date2num(datetime.strptime(opts.tmax,'%Y%m%d'))
if opts.data_tmin is not None:
    dmin = datetime.strptime(opts.data_tmin,'%Y%m%d')
else:
    dmin = num2date(nmin+opts.tstr-opts.tmgn).replace(tzinfo=None)
if opts.data_tmax is not None:
    dmax = datetime.strptime(opts.data_tmax,'%Y%m%d')
else:
    dmax = num2date(nmax+opts.tend+opts.tmgn).replace(tzinfo=None)
data_nmin = date2num(dmin)
data_nmax = date2num(dmax)

values = []
labels = []
ticks = []
for y in range(dmin.year,dmax.year+1):
    for m in range(1,13,2):
        d = datetime(y,m,1)
        values.append(date2num(d))
        labels.append(d.strftime('%Y-%m'))
    for m in range(1,13,1):
        d = datetime(y,m,1)
        ticks.append(date2num(d))

ndvi_flag = np.load('cloud_flag.npy').swapaxes(0,1)
ndvi_data[ndvi_flag] = np.nan
nobject = ndvi_data.shape[1]
data = np.load(opts.cflag_fnam)
cloud_ntim = data['ntim']
cloud_flag = data['cloud_flag']

ndvi_dtim = []
ndvi_data = []
ndvi_flag = []
nobject = None
fs = sorted(glob(os.path.join(opts.datdir,'^atcor_data_ndvi_'+'[0-9]'*8+'\.npz')))
for fnam in fs:
    f = os.path.basename(fnam)
    m = re.search('^atcor_data_ndvi_('+'\d'*8+')\.npz',f)
    if not m:
        raise ValueError('Error in finding date >>> '+f)
    dstr = m.group(1)
    d = datetime.strptime(dstr,'%Y%m%d')
    if d < dmin or d > dmax:
        continue
    sys.stderr.write(f+' '+dstr+'\n')
    dtmp = np.load(fnam)['data_cor']
    if nobject is None:
        nobject = dtmp[0].size
        if nobject != cloud_flag.shape[]
    elif dtmp[0].size != nobject:
        raise ValueError('Error, dtmp[0].size={}, nobject={}'.format(dtmp[0].size,nobject))
    ndvi_dtim.append(d)
    ndvi_data.append(dtmp[i])
ndvi_dtim = np.array(ndvi_dtim)
ndvi_data = np.array(ndvi_data)
ndvi_ntim = date2num(ndvi_dtim)

xx = np.arange(np.floor(data_nmin),np.ceil(data_nmax)+0.1*opts.tstp,opts.tstp)
if not opts.batch:
    plt.interactive(True)
fig = plt.figure(1,facecolor='w')
plt.subplots_adjust(left=0.12,right=0.88,bottom=0.12,top=0.80,hspace=0.5)
if not opts.batch:
    pdf = PdfPages('example_estimation_early.pdf')
for iobj in range(nobject):
#for i in [120]:
    object_id = iobj+1
    zi = ndvi_data[:,iobj] # NDVI
    cnd = ~np.isnan(zi)
    if cnd.sum() < 1:
        continue
    zz = csaps(ndvi_ntim[cnd],zi[cnd],xx,smooth=opts.smooth)
    g1 = csaps(xx,np.gradient(zz)/np.gradient(xx),xx,smooth=opts.smooth)*100.0
    g2 = np.gradient(g1)/np.gradient(xx)*2.0
    max_peaks,properties = find_peaks(g2,distance=opts.ndvi2_distance,prominence=opts.ndvi2_prominence)
    xest = xx[max_peaks]
    yest = g2[max_peaks]
    yrgt = []
    zrgt = []
    #ylft = []
    zlft = []
    for ic in range(xest.size):
        ix = np.argmin(np.abs(xx-xest[ic]))
        cnd = (xx > xest[ic]) & (xx < xest[ic]+45)
        yrgt.append(zz[cnd].max()-zz[ix])
        zrgt.append(zz[ix]-zz[cnd].min())
        cnd = (xx < xest[ic]) & (xx > xest[ic]-45)
        #ylft.append(zz[cnd].max()-zz[ix])
        zlft.append(zz[ix]-zz[cnd].min())
    yrgt = np.array(yrgt)
    zrgt = np.array(zrgt)
    #ylft = np.array(ylft)
    zlft = np.array(zlft)
    zz_inc = []
    for ix in range(xx.size):
        cnd = (xx > xx[ix]) & (xx < xx[ix]+45)
        if cnd.sum() > 0:
            zz_inc.append(zz[cnd].max()-zz[ix])
        else:
            zz_inc.append(0.0)
    zz_inc = np.array(zz_inc)
    ss1 = -((zz-opts.ndvi_mean)/ndvi_std)*0.08
    cnd = (zz < opts.ndvi_mean)
    ss1[cnd] = -((zz[cnd]-opts.ndvi_mean)/ndvi_std)*opts.fact1
    ss2 = ((zz_inc-ndvi_inc_mean)/ndvi_inc_std)*0.08
    ss3 = ((g2-ndvi_dd_mean)/ndvi_dd_std)*0.08
    ss = (ss1+ss2+ss3)*0.4

    yval = []
    yflg = []
    for ic in range(xest.size):
        ix = np.argmin(np.abs(xx-xest[ic]))
        if (zz[ix] < 0.1) or (zz[ix] > 0.5):
            y = -1.0e10
            f = False
        elif (g1[ix] < -1.0):
            y = -1.0e10
            f = False
        elif (yrgt[ic] < 0.1):
            y = -1.0e10
            f = False
        elif (zrgt[ic] > 0.03):
            y = -1.0e10
            f = False
        elif (zlft[ic] > 0.3):
            y = -1.0e10
            f = False
        elif (yest[ic] < 0.05):
            y = -1.0e10
            f = False
        elif (xest[ic] <= nmin) or (xest[ic] >= nmax):
            y = -1.5e10
            f = False
        else:
            y = yest[ic]
            f = True
        yval.append(y)
        yflg.append(f)
    yval = np.array(yval)
    yflg = np.array(yflg)

    for ic in range(0,xest.size-1):
        ix = np.argmin(np.abs(xx-xest[ic]))
        ix2 = np.argmin(np.abs(xx-xest[ic+1]))
        if yflg[ic] and yflg[ic+1] and (xest[ic+1]-xest[ic] < opts.tmrg) and (zz[ix] < zz[ix2]) and (ss[ix] < ss[ix2]):
            ss[ix] += (ss[ix2]-ss[ix])+1.0e-4

    sval = []
    for ic in range(xest.size):
        if yflg[ic]:
            ix = np.argmin(np.abs(xx-xest[ic]))
            s = ss[ix]
        else:
            s = yval[ic]
        sval.append(s)
    sval = np.array(sval)

    # Select three candidates
    cval = sval.copy()
    can_inds = []
    indv = np.argsort(cval)[::-1]
    if cval[indv[0]] > -10.0:
        can_inds.append(indv[0])
        cnd = np.abs(xest-xest[indv[0]]) < 1.0
        cval[cnd] = -2.0e10
        indv = np.argsort(cval)[::-1]
    if cval[indv[0]] > -10.0:
        can_inds.append(indv[0])
        cnd = np.abs(xest-xest[indv[0]]) < 1.0
        cval[cnd] = -3.0e10
        indv = np.argsort(cval)[::-1]
    if cval[indv[0]] > -10.0:
        can_inds.append(indv[0])

    # Sort three candidates
    xans = []
    yans = []
    cval = sval.copy()
    indv = np.argsort(cval)[::-1]
    if cval[indv[0]] > -10.0:
        xans.append(xest[indv[0]])
        yans.append(yest[indv[0]])
        cnd = np.abs(xest-xest[indv[0]]) < 1.0
        cval[cnd] = -2.0e10
        indv = np.argsort(cval)[::-1]
    else:
        xans.append(np.nan)
        yans.append(np.nan)
    if cval[indv[0]] > -10.0:
        xans.append(xest[indv[0]])
        yans.append(yest[indv[0]])
        cnd = np.abs(xest-xest[indv[0]]) < 1.0
        cval[cnd] = -3.0e10
        indv = np.argsort(cval)[::-1]
    else:
        xans.append(np.nan)
        yans.append(np.nan)
    if cval[indv[0]] > -10.0:
        xans.append(xest[indv[0]])
        yans.append(yest[indv[0]])
    else:
        xans.append(np.nan)
        yans.append(np.nan)
    xans = np.array(xans)
    yans = np.array(yans)

    if opts.debug:
        fig.clear()
        ax1 = plt.subplot(111)
        ax1.minorticks_on()
        #ax1.plot(ndvi_ntim,zi,'.',color='#888888')
        l1, = ax1.plot(xx,zz,'-',color='k',label='NDVI')
        l2, = ax1.plot(xx,g2,'-',color='b',label='NDVI$^{\prime\prime}$')
        l3, = ax1.plot(xx,ss1+0.6,'-',color='#666633',label='S$_1$')
        l4, = ax1.plot(xx,ss2+0.6,'-',color='#999900',label='S$_2$')
        l5, = ax1.plot(xx,ss3+0.6,'-',color='#cccc66',label='S$_3$')
        l6, = ax1.plot(xx,ss+0.6,'-',color='#996600',label='S$_{t}$')

        for ic in range(xest.size):
            if yflg[ic]:
                ix = np.argmin(np.abs(xx-xest[ic]))
                y1 = yest[ic]
                y2 = ss[ix]+0.6
                ax1.plot(xest[ic],yest[ic],'bo')
                ax1.vlines(xest[ic],min(y1,y2),max(y1,y2),color='k',ls=':')
            else:
                ax1.plot(xest[ic],yest[ic],'bo')
                ax1.plot(xest[ic],yest[ic],'bx',ms=20)

        cols = ['r','m','c','b','k']
        indv = np.argsort(sval)[::-1]
        for iv in range(xest.size):
            if yval[indv[iv]] > -1000.0:
                ic = indv[iv]
                ix = np.argmin(np.abs(xx-xest[ic]))
                #ax1.plot(xest[ic],yest[ic],'*',color=cols[iv%len(cols)])
                ax1.text(xx[ix],ss[ix]+0.6,'{}'.format(iv+1),ha='center',va='bottom',size=16)
        ax1.plot(xans[2],yans[2],'o',ms=20,mfc='none',color='orange',mew=2,zorder=9)
        ax1.plot(xans[1],yans[1],'o',ms=20,mfc='none',color=cols[1],mew=2,zorder=9)
        ax1.plot(xans[0],yans[0],'o',ms=20,mfc='none',color=cols[0],mew=2,zorder=9)
        l12 = ax1.vlines(xans[2],ndvi_min,yans[2],color='orange',label='T$_{est3}$',zorder=10)
        l11 = ax1.vlines(xans[1],ndvi_min,yans[1],color=cols[1],label='T$_{est2}$',zorder=10)
        l10 = ax1.vlines(xans[0],ndvi_min,yans[0],color=cols[0],label='T$_{est1}$',zorder=10)
        #ax1.plot(xans[0],yans[0],'r*',ms=10)
        #ax1.plot(xans[1],yans[1],'m*',ms=10)

        ax1.set_title('OBJECTID= {}'.format(object_id),y=1.15,x=0.5)
        ax1.set_xlim(xx.min(),xx.max())
        ax1.set_ylim(ndvi_min,ndvi_max)
        ax1.xaxis.set_major_locator(plt.FixedLocator(values))
        ax1.xaxis.set_major_formatter(plt.FixedFormatter(labels))
        ax1.xaxis.set_minor_locator(plt.FixedLocator(ticks))
        #ax1.yaxis.set_major_locator(plt.MultipleLocator(5.0))
        #ax1.set_ylabel('NDVI',color='#888888')
        ax1.set_ylabel(r'NDVI, NDVI$^{\prime\prime}$ $\times$ 200, S$_{x}$ $+$ 0.6')
        ax1.yaxis.set_label_coords(-0.105,0.5)
        for l in ax1.xaxis.get_ticklabels():
            l.set_rotation(30.0)
        lns = [l1,l2,l3,l4,l5,l6,l10,l11,l12]
        lbs = [l.get_label() for l in lns]
        ax1.legend(lns,lbs,prop={'size':12},numpoints=1,loc=8,bbox_to_anchor=(0.5,1.01),ncol=len(lns),frameon=False,handletextpad=0.05,columnspacing=0.40,handlelength=0.8)

        ax1.fill_betweenx([ndvi_min,ndvi_max],xx.min(),nmin,color='k',alpha=0.1)
        ax1.fill_betweenx([ndvi_min,ndvi_max],nmax,xx.max(),color='k',alpha=0.1)

        if not opts.batch:
            plt.savefig(pdf,format='pdf')
            plt.draw()
            plt.pause(0.1)
    #break
if opts.debug:
    if not opts.batch:
        pdf.close()
