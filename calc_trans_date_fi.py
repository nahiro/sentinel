#!/usr/bin/env python
import os
import sys
import re
from glob import glob
from datetime import datetime
import gdal
import osr
import numpy as np
import cartopy.io.shapereader as shpreader
from matplotlib.dates import date2num,num2date
from csaps import csaps
from statsmodels.stats.weightstats import DescrStatsW
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from optparse import OptionParser,IndentedHelpFormatter

# Default values
TMIN = '20190501'
TMAX = '20190915'
TMGN = 30.0 # day
TSTP = 0.1 # day
SMOOTH = 0.01
OFFSET_V = 4.8
OFFSET_FACT = 1.25
DATDIR = '.'
SEARCH_KEY = 'resample'
AREA_FNAM = 'pixel_area_block.dat'

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('-s','--tmin',default=TMIN,help='Min date of transplanting in the format YYYYMMDD (%default)')
parser.add_option('-e','--tmax',default=TMAX,help='Max date of transplanting in the format YYYYMMDD (%default)')
parser.add_option('--data_tmin',default=DATA_TMIN,help='Min date of input data in the format YYYYMMDD (%default)')
parser.add_option('--data_tmax',default=DATA_TMAX,help='Max date of input data in the format YYYYMMDD (%default)')
parser.add_option('--tmgn',default=TMGN,type='float',help='Margin of input data in day (%default)')
parser.add_option('--tstp',default=TSTP,type='float',help='Precision of transplanting date in day (%default)')
parser.add_option('--offset_v',default=OFFSET_V,type='float',help='Max offset-free VH_BSC difference (%default)')
parser.add_option('--offset_fact',default=OFFSET_FACT,type='float',help='Offset factor for VH_BSC difference (%default)')
parser.add_option('-I','--incidence_angle',default=None,help='Incidence angle file, format: date(%Y%m%d) angle(deg) (%default)')
parser.add_option('-l','--incidence_list',default=None,help='Incidence angle list, format: flag(0|1=baseline) pol(VH|VV) angle(deg) filename (%default)')
parser.add_option('-S','--smooth',default=SMOOTH,type='float',help='Smoothing factor from 0 to 1 (%default)')
parser.add_option('--output_epsg',default=None,type='int',help='Output EPSG (guessed from input data)')
parser.add_option('--area_fnam',default=AREA_FNAM,help='Pixel area file name (%default)')
parser.add_option('-D','--datdir',default=DATDIR,help='Input data directory, not used if input_fnam is given (%default)')
parser.add_option('--search_key',default=SEARCH_KEY,help='Search key for input data, not used if input_fnam is given (%default)')
parser.add_option('-i','--inp_fnam',default=None,help='Input GeoTIFF name (%default)')
(opts,args) = parser.parse_args()

nmin = date2num(datetime.strptime(opts.tmin,'%Y%m%d'))
nmax = date2num(datetime.strptime(opts.tmax,'%Y%m%d'))
if opts.data_tmin is not None:
    dmin = datetime.strptime(opts.data_tmin,'%Y%m%d')
else:
    dmin = num2date(nmin-opts.tmgn).replace(tzinfo=None)
if opts.data_tmax is not None:
    dmax = datetime.strptime(opts.data_tmax,'%Y%m%d')
else:
    dmax = num2date(nmax+opts.tmgn).replace(tzinfo=None)

values = []
labels = []
ticks = []
for y in range(2017,2022):
    for m in range(1,13,2):
        d = datetime(y,m,1)
        values.append(date2num(d))
        labels.append(d.strftime('%Y-%m'))
    for m in range(1,13,1):
        d = datetime(y,m,1)
        ticks.append(date2num(d))

object_ids = []
blocks = []
inds = []
areas = []
with open(opts.area_fnam,'r') as fp:
    for line in fp:
        item = line.split()
        if len(item) < 3 or item[0] == '#':
            continue
        object_ids.append(int(item[0]))
        blocks.append(item[1])
        inds.append([])
        areas.append([])
        n = int(item[2])
        if len(item) < 5 and n != 0:
            raise ValueError('Error, len(item)={}, n={}, expected n=0.'.format(len(item),n))
        for nn in range(3,n*2+3,2):
            inds[-1].append(int(item[nn]))
            areas[-1].append(float(item[nn+1]))
        inds[-1] = np.array(inds[-1])
        areas[-1] = np.array(areas[-1])
object_ids = np.array(object_ids)
blocks = np.array(blocks)
inds = np.array(inds)
areas = np.array(areas)
nobject = object_ids.size

if opts.inp_fnam is not None:
    ds = gdal.Open(opts.inp_fnam)
    if opts.output_epsg is None:
        prj = ds.GetProjection()
        srs = osr.SpatialReference(wkt=prj)
        epsg = srs.GetAttrValue('AUTHORITY',1)
        if re.search('\D',epsg):
            raise ValueError('Error in EPSG >>> '+epsg)
        output_epsg = int(epsg)
    else:
        output_epsg = opts.output_epsg
    data = ds.ReadAsArray()
    data_shape = data[0].shape
    data_trans = ds.GetGeoTransform()
    band_list = []
    for i in range(len(data)):
        band_list.append(ds.GetRasterBand(i+1).GetDescription())
    band_list = np.array(band_list)
    ds = None
else:
    output_epsg = None
    data_shape = None
    data_trans = None
    data = []
    band_list = []
    fs = sorted(glob(os.path.join(opts.datdir,'*'+'[0-9]'*8+'*.tif')))
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
        if d < dmin or d > dmax:
            continue
        sys.stderr.write(f+' '+dstr+'\n')
        ds = gdal.Open(fnam)
        if output_epsg is None:
            if opts.output_epsg is None:
                prj = ds.GetProjection()
                srs = osr.SpatialReference(wkt=prj)
                epsg = srs.GetAttrValue('AUTHORITY',1)
                if re.search('\D',epsg):
                    raise ValueError('Error in EPSG >>> '+epsg)
                output_epsg = int(epsg)
            else:
                output_epsg = opts.output_epsg
        dtmp = ds.ReadAsArray()
        if data_shape is None:
            data_shape = dtmp[0].shape
        elif dtmp[0].shape != data_shape:
            raise ValueError('Error, dtmp[0].shape={}, data_shape={}'.format(dtmp[0].shape,data_shape))
        trans = ds.GetGeoTransform()
        if data_trans is None:
            data_trans = trans
        elif trans != data_trans:
            raise ValueError('Error, trans={}, data_trans={}'.format(trans,data_trans))
        for i in range(len(dtmp)):
            data.append(dtmp[i])
            band_list.append(ds.GetRasterBand(i+1).GetDescription())
        ds = None
    data = np.array(data)
    band_list = np.array(band_list)

if opts.incidence_list is not None:
    incidence_flag = []
    incidence_angle = []
    incidence_fnam = []
    incidence_select = []
    with open (opts.incidence_list,'r') as fp:
        for line in fp:
            item = line.split()
            if len(item) < 4:
                continue
            if item[0][0] == '#':
                continue
            pol = item[1].upper()
            if pol != 'VH':
                continue
            incidence_flag.append(int(item[0]))
            incidence_angle.append(float(item[2]))
            incidence_fnam.append(item[3])
            if len(item) > 4:
                if item[4][0].lower() == 'f':
                    incidence_select.append(False)
                else:
                    incidence_select.append(True)
            else:
                incidence_select.append(True)
    incidence_flag = np.array(incidence_flag)
    incidence_angle = np.array(incidence_angle)
    incidence_select = np.array(incidence_select)
    nangle = incidence_angle.size
    cnd = (incidence_flag == 1)
    if cnd.sum() != 1:
        raise ValueError('Error in incidence_flag, cnd.sum()={}'.format(cnd.sum()))
    baseline_indx = np.arange(nangle)[cnd][0]
    signal_avg = []
    for fnam in incidence_fnam:
        avg = np.load(fnam)
        if avg.shape != data_shape:
            raise ValueError('Error, avg.shape={}, data_shape={}, fnam={}'.format(avg.shape,data_shape,fnam))
        signal_avg.append(avg)
    signal_dif = []
    for i in range(nangle):
        if i == baseline_indx:
            signal_dif.append(0.0)
        else:
            signal_dif.append(signal_avg[baseline_indx]-signal_avg[i])
    incidence_indx = {}
    with open(opts.incidence_angle,'r') as fp:
        for line in fp:
            item = line.split()
            if len(item) < 2:
                continue
            if item[0][0] == '#':
                continue
            ang = float(item[1])
            dif = np.abs(incidence_angle-ang)
            indx = np.argmin(dif)
            if dif[indx] > 0.1:
                sys.stderr.write('Warning, ang={}, dif={}\n'.format(ang,dif[indx]))
                continue
            incidence_indx.update({item[0]:indx})

vh_dtim = []
vh_data = []
for i,band in enumerate(band_list):
    if not re.search('VH',band):
        continue
    m = re.search('_('+'\d'*8+')$',band)
    if not m:
        raise ValueError('Error in finding date >>> '+band)
    dstr = m.group(1)
    if opts.incidence_list is not None:
        if not incidence_select[incidence_indx[dstr]]:
            continue
        dtmp = (data[i]+signal_dif[incidence_indx[dstr]]).flatten()
    else:
        dtmp = data[i].flatten()
    data_avg = []
    for i in range(nobject):
        if inds[i].size < 1:
            data_avg.append(np.nan)
            continue
        data_value = dtmp[inds[i]]
        data_weight = areas[i]
        cnd = ~np.isnan(data_value)
        if cnd.sum() <= 1:
            data_avg.append(data_value[cnd].mean())
        else:
            data_weighted_stats = DescrStatsW(data_value[cnd],weights=data_weight[cnd],ddof=0)
            data_avg.append(data_weighted_stats.mean)
    data_avg = np.array(data_avg)
    vh_dtim.append(datetime.strptime(dstr,'%Y%m%d'))
    vh_data.append(data_avg)
vh_dtim = np.array(vh_dtim)
vh_data = np.array(vh_data)
vh_ntim = date2num(vh_dtim)

y_profile = np.load('y_profile.npy')

xx = np.arange(np.floor(vh_ntim.min()),np.ceil(vh_ntim.max())+0.1*opts.tstp,opts.tstp)
plt.interactive(True)
fig = plt.figure(1,facecolor='w')
plt.subplots_adjust(left=0.12,right=0.88,bottom=0.12,top=0.80,hspace=0.5)
pdf = PdfPages('example_estimation.pdf')
for i in range(nobject):
    object_id = object_ids[i]
    if object_id != i+1:
        raise ValueError('Error, object_id={}, i={}'.format(object_id,i))
    if not object_id in object_id_check:
        continue
    if inds[i].size < 1:
        continue
    yi = vh_data[:,i] # VH
    yy = csaps(vh_ntim,yi,xx,smooth=opts.smooth)

    yy_max = yy.max()
    yy_min = yy.min()
    yy_thr = max(yy_min+0.3*(yy_max-yy_min),-21.0)
    #fp_thr = 60.0
    fp_thr = 30.0
    cc = (yy < yy_thr)
    fp_inds = []
    i1 = None
    i2 = None
    flag = False
    for ic in range(cc.size):
        if not flag and cc[ic]:
            i1 = ic
            i2 = None
            flag = True
        elif flag and not cc[ic]:
            i2 = ic
            if i1 is None:
                raise ValueError('Error, i1={}, i2={}'.format(i1,i2))
            if xx[i1:i2].ptp() > fp_thr:
                fp_inds.append([i1,i2])
            i1 = None
            i2 = None
            flag = False
    if i1 is not None:
        i2 = cc.size
        if xx[i1:i2].ptp() > fp_thr:
            fp_inds.append([i1,i2])
    fp_inds = np.array(fp_inds)
    ff1 = np.zeros(xx.size)
    for f in fp_inds:
        j1 = f[0]
        j2 = f[1]
        ff1[j1:j2] = (yy_thr-yy[j1:j2]).sum()
    fishpond_index = (ff1/2000.0).clip(0.0,1.0)

    ss1 = []
    ss2 = []
    ss3 = []
    for i in range(xx.size):
        dmax = yy.max()-yy[i]
        if (dmax > opts.offset_v) and (fishpond_index[i] > 0.01):
            offset = (dmax-opts.offset_v)*opts.offset_fact
        else:
            offset = 0.0
        ytmp = yy-(yy[i]+offset)
        i2 = min(i+y_profile.size,yy.size)
        if i2 > i:
            yfact = y_profile[0:i2-i].copy()
            yfact *= yfact.size/yfact.sum()
            ytmp[i:i2] *= yfact
        cnd = (xx > xx[i]) & (xx < xx[i]+30)
        y2 = ytmp[cnd].mean()
        ss1.append(y2)
        cnd = (xx > xx[i]) & (xx < xx[i]+60)
        y2 = ytmp[cnd].mean()
        ss2.append(y2)
        cnd = (xx > xx[i]) & (xx < xx[i]+90)
        y2 = ytmp[cnd].mean()
        ss3.append(y2)
    ss1 = np.array(ss1)
    ss2 = np.array(ss2)
    ss3 = np.array(ss3)
    f1 = np.gradient(yy)/np.gradient(xx)
    f2 = np.gradient(f1)/np.gradient(xx)
    cc = (f1 >= -0.1) & (f1 <= 0.1) & (f2 > 0.0)
    #cc = (f2 > 0.0)

    xcs = []
    ycs = []
    f1cs = []
    f2cs = []
    fflg = []
    i1 = None
    i2 = None
    flag = False
    for ic in range(cc.size):
        if not flag and cc[ic]:
            i1 = ic
            i2 = None
            flag = True
        elif flag and not cc[ic]:
            i2 = ic
            if i1 is None:
                raise ValueError('Error, i1={}, i2={}'.format(i1,i2))
            xcs.append(xx[i1:i2])
            ycs.append(yy[i1:i2])
            f1cs.append(f1[i1:i2])
            f2cs.append(f2[i1:i2])
            if f1[i1:i2].min()*f1[i1:i2].max() < 0.0:
                fflg.append(True)
            else:
                fflg.append(False)
            i1 = None
            i2 = None
            flag = False
    if i1 is not None:
        i2 = cc.size
        xcs.append(xx[i1:i2])
        ycs.append(yy[i1:i2])
        f1cs.append(f1[i1:i2])
        f2cs.append(f2[i1:i2])
        if f1[i1:i2].min()*f1[i1:i2].max() < 0.0:
            fflg.append(True)
        else:
            fflg.append(False)
    fflg = np.array(fflg)

    fig.clear()
    ax1 = plt.subplot(111)
    ax1.minorticks_on()
    ax2 = ax1.twinx()
    ax3 = ax1.twinx()
    ax3.set_yticks([])
    #ax1.set_zorder(0)
    #ax2.set_zorder(-1)
    #ax1.plot(vh_ntim,yi,'k.')
    l1, = ax1.plot(xx,yy,'k-',label='BSC')
    ss3[-10:] = np.nan
    l2, = ax3.plot(xx,ss3,'y-',lw=1,label='Signal')
    #ax3.plot(xx,ss1,'y-')
    #ax3.plot(xx,ss2,'-',color='orange')
    #ax3.plot(xx,ss3,'-',color='brown')
    
    xest = []
    yest = []
    yflg = []
    for ic in range(fflg.size):
        if f1cs[ic].min()*f1cs[ic].max() < 0:
            indc = np.argmin(np.abs(f1cs[ic]))
        else:
            indc = np.argmax(f2cs[ic])
        xest.append(xcs[ic][indc])
        yest.append(ycs[ic][indc])
        yflg.append(False)
    xest = np.array(xest)
    yest = np.array(yest)
    yflg = np.array(yflg)
    ydif = []
    yinc = []
    ylft = []
    yrgt = []
    for ic in range(fflg.size):
        cnd = (xx > xest[ic]-60) & (xx < xest[ic])
        ylft.append(yy[cnd].max()-yest[ic])
        cnd = (xx > xest[ic]) & (xx < xest[ic]+60)
        yrgt.append(yy[cnd].max()-yest[ic])
        cnd = (xest > xest[ic]) & (fflg)
        if cnd.sum() > 0:
            ix = np.where(cnd)[0][0]
            cnd = (xx > xest[ic]) & (xx < xest[ix])
            yd = yy[cnd].max()-yest[ic]
            ydif.append(yd)
            yinc.append(yest[ix]-yest[ic])
            if (yd > 2.0) & (yest[ix] > yest[ic]+1.0) & (xest[ix] < xest[ic]+60.0):
                yflg[ix] = True
        else:
            cnd = (xx > xest[ic])
            yd = yy[cnd].max()-yest[ic]
            ydif.append(yd)
            yinc.append(yd)
    ydif = np.array(ydif)
    yinc = np.array(yinc)
    ylft = np.array(ylft)
    yrgt = np.array(yrgt)
    sval = []
    for ic in range(fflg.size):
        #if fflg[ic]:
        #    ax1.plot(xest[ic],yest[ic],'ro')
        #else:
        #    ax1.plot(xest[ic],yest[ic],'bo')
        ix = np.argmin(np.abs(xx-xest[ic]))
        if yflg[ic]:
            #ax1.plot(xcs[ic],ycs[ic],'m-')
            #ax1.plot(xest[ic],yest[ic],'m.')
            #ax1.plot(xest[ic],yest[ic],'ro',ms=20,mfc='none')
            if fflg[ic]:
                ax1.plot(xest[ic],yest[ic],'ro')
            else:
                ax1.plot(xest[ic],yest[ic],'bo')
            s = ss3[ix]
            y1 = yest[ic]
            y2 = s-bsc_ofs
            ax1.vlines(xest[ic],min(y1,y2),max(y1,y2),color='k',ls=':')
        elif yrgt[ic] < 1.75:
            #ax1.plot(xcs[ic],ycs[ic],'c-')
            #ax1.plot(xest[ic],yest[ic],'c.')
            if fflg[ic]:
                ax1.plot(xest[ic],yest[ic],'ro')
            else:
                ax1.plot(xest[ic],yest[ic],'bo')
            ax1.plot(xest[ic],yest[ic],'rx',ms=20)
            s = -1.0e10
        elif (ic < fflg.size-1) and (ydif[ic]-yinc[ic] > 4.0) and (xest[ic+1]-xest[ic] < 60.0):
            #ax1.plot(xcs[ic],ycs[ic],'c-')
            #ax1.plot(xest[ic],yest[ic],'c.')
            if fflg[ic]:
                ax1.plot(xest[ic],yest[ic],'ro')
            else:
                ax1.plot(xest[ic],yest[ic],'bo')
            ax1.plot(xest[ic],yest[ic],'rx',ms=20)
            s = -1.0e10
        elif (ydif[ic] > 0.9) or (yinc[ic] > 0.75) or ((ylft[ic] > 3.0) & (yrgt[ic] > 3.0)):
            #ax1.plot(xcs[ic],ycs[ic],'r-')
            #ax1.plot(xest[ic],yest[ic],'r.')
            #ax1.plot(xest[ic],yest[ic],'ro',ms=20,mfc='none')
            if fflg[ic]:
                ax1.plot(xest[ic],yest[ic],'ro')
            else:
                ax1.plot(xest[ic],yest[ic],'bo')
            s = ss3[ix]
            y1 = yest[ic]
            y2 = s-bsc_ofs
            ax1.vlines(xest[ic],min(y1,y2),max(y1,y2),color='k',ls=':')
        else:
            #ax1.plot(xcs[ic],ycs[ic],'r-')
            #ax1.plot(xest[ic],yest[ic],'r.')
            #ax1.plot(xest[ic],yest[ic],'ro',ms=20,mfc='none')
            if fflg[ic]:
                ax1.plot(xest[ic],yest[ic],'ro')
            else:
                ax1.plot(xest[ic],yest[ic],'bo')
            s = ss3[ix]
            y1 = yest[ic]
            y2 = s-bsc_ofs
            ax1.vlines(xest[ic],min(y1,y2),max(y1,y2),color='k',ls=':')
        if (xest[ic] <= nmin) or (xest[ic] >= nmax):
            s = -1.5e10
        sval.append(s)
    sval = np.array(sval)
    indv = np.argsort(sval)[::-1]

    cols = ['r','m','c','b','k']
    for iv in range(fflg.size):
        if sval[indv[iv]] > -1000.0:
            ic = indv[iv]
            ix = np.argmin(np.abs(xx-xest[ic]))
            #ax1.plot(xest[ic],yest[ic],'*',color=cols[iv%len(cols)])
            ax3.text(xx[ix],ss3[ix]+0.1,'{}'.format(iv+1),ha='center',va='bottom',size=16)
            #if iv < 2:
            #    ax1.plot(xest[ic],yest[ic],'o',ms=20,mfc='none',color=cols[iv%len(cols)],mew=2)
            #    ax1.vlines(xest[ic],bsc_min,yest[ic],color=cols[iv%len(cols)])
            #ax3.plot(xx[ix],ss3[ix],'o',color=cols[iv%len(cols)])

    xans = []
    yans = []
    if sval[indv[0]] > -10.0:
        xans.append(xest[indv[0]])
        yans.append(yest[indv[0]])
        cnd = np.abs(xest-xest[indv[0]]) < 1.0
        sval[cnd] = -2.0e10
        indv = np.argsort(sval)[::-1]
    else:
        xans.append(np.nan)
        yans.append(np.nan)
    if sval[indv[0]] > -10.0:
        xans.append(xest[indv[0]])
        yans.append(yest[indv[0]])
        cnd = np.abs(xest-xest[indv[0]]) < 1.0
        sval[cnd] = -3.0e10
        indv = np.argsort(sval)[::-1]
    else:
        xans.append(np.nan)
        yans.append(np.nan)
    if sval[indv[0]] > -10.0:
        xans.append(xest[indv[0]])
        yans.append(yest[indv[0]])
    else:
        xans.append(np.nan)
        yans.append(np.nan)
    xans = np.array(xans)
    yans = np.array(yans)
    ax1.plot(xans[2],yans[2],'o',ms=20,mfc='none',color='orange',mew=2,zorder=9)
    ax1.plot(xans[1],yans[1],'o',ms=20,mfc='none',color=cols[1],mew=2,zorder=9)
    ax1.plot(xans[0],yans[0],'o',ms=20,mfc='none',color=cols[0],mew=2,zorder=9)
    l9 = ax1.vlines(xans[2],bsc_min,yans[2],color='orange',label='T$_{est3}$',zorder=10)
    l8 = ax1.vlines(xans[1],bsc_min,yans[1],color=cols[1],label='T$_{est2}$',zorder=10)
    l7 = ax1.vlines(xans[0],bsc_min,yans[0],color=cols[0],label='T$_{est1}$',zorder=10)
    #ax1.plot(xans[0],yans[0],'r*',ms=10)
    #ax1.plot(xans[1],yans[1],'m*',ms=10)

    #ax2.plot(ndvi_ntim,zi,'.',color='#888888')
    l3, = ax2.plot(xx,zz,'-',color='#888888',label='NDVI')
    l4, = ax2.plot(xx,fishpond_index,'-',color='#cccccc',label='FI',zorder=0)
    indx = list(object_id_check).index(object_id)
    trans_date_uncorrected = trans_date_check[indx]+9.0
    trans_date_corrected = trans_date_check[indx]
    trans_date_survey = trans_date_true_check[indx]
    #ax1.axvline(trans_date_uncorrected,color='b',linestyle=':',label='T$_{uncorrected}$') # back correct for the offset
    #ax1.axvline(trans_date_corrected,color='b',linestyle='-',label='T$_{corrected}$')
    l5 = ax1.axvline(trans_date_survey,color='g',label='T$_{survey}$')
    #if plot_check[indx] != 120:
        #ax1.plot(date2num(datetime.strptime('2019/'+survey_check[indx],'%Y/%m/%d'))-plot_check[indx],-23.75,'rv')
        #ax1.plot(trans_date_survey,-23.65,'gv',ms=10)
    inde = np.argmin(np.abs(trans_date_survey-xest))
    trans_date_estimate = xest[inde]
    trans_signal_estimate = yest[inde]
    trans_date_answer1 = xans[0]
    trans_signal_answer1 = yans[0]
    trans_date_answer2 = xans[1]
    trans_signal_answer2 = yans[1]
    trans_date_answer3 = xans[2]
    trans_signal_answer3 = yans[2]

    ax1.set_title('OBJECTID= {}, $\epsilon$= {:.1f}, Site({})={}, Growth={}'.format(object_id,trans_date_estimate-trans_date_true_check[indx],survey_check[indx],site_check[indx],plot_check[indx]),y=1.15,x=0.5)
    ax1.set_ylim(bsc_min,bsc_max)
    ax2.set_ylim(-0.2,1.1)
    ax3.set_ylim(bsc_min+bsc_ofs,bsc_max+bsc_ofs)
    ax1.set_xlim(xx.min(),xx.max())
    #ax2.set_xlim(xx.min(),xx.max())
    ax1.xaxis.set_major_locator(plt.FixedLocator(values))
    ax1.xaxis.set_major_formatter(plt.FixedFormatter(labels))
    ax1.xaxis.set_minor_locator(plt.FixedLocator(ticks))
    ax1.yaxis.set_major_locator(plt.MultipleLocator(5.0))
    ax1.set_ylabel('BSC, Signal $-$ {} (dB)'.format(int(bsc_ofs+0.1)))
    #ax2.set_ylabel('NDVI',color='#888888')
    ax2.set_ylabel('NDVI, FI')
    ax1.yaxis.set_label_coords(-0.105,0.5)
    ax2.yaxis.set_label_coords(1.10,0.5)
    for l in ax1.xaxis.get_ticklabels():
        l.set_rotation(30.0)
    lns = [l1,l2,l3,l4,l5,l6,l7,l8,l9]
    lbs = [l.get_label() for l in lns]
    ax1.legend(lns,lbs,prop={'size':12},numpoints=1,loc=8,bbox_to_anchor=(0.5,1.01),ncol=9,frameon=False,handletextpad=0.1,columnspacing=0.60,handlelength=1.2)

    plt.savefig(pdf,format='pdf')
    plt.draw()
    #plt.pause(0.1)
    #break
pdf.close()
