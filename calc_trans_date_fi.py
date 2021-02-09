#!/usr/bin/env python
import os
import sys
import shutil
import re
import warnings
from glob import glob
from datetime import datetime
import gdal
import osr
import json
import shapefile
from collections import OrderedDict
import numpy as np
from matplotlib.dates import date2num,num2date
from csaps import csaps
from statsmodels.stats.weightstats import DescrStatsW
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from optparse import OptionParser,IndentedHelpFormatter

# Default values
#TMIN = '20190501'
#TMAX = '20190915'
TMIN = '20190401'
TMAX = '20190615'
TMGN = 30.0 # day
TSTP = 0.1 # day
SMOOTH = 0.01
SIGNAL_V = 4.8 # dB
SIGNAL_FACT = 1.25
OFFSET = 0.0 # day
DATDIR = '.'
SEARCH_KEY = 'resample'
X_PROFILE = 'x_profile.npy'
Y_PROFILE = 'y_profile.npy'
AREA_FNAM = 'pixel_area_block.dat'
HOME = os.environ.get('HOME')
if HOME is None:
    HOME = os.environ.get('HOMEPATH')
SHP_FNAM = os.path.join(HOME,'Work','SATREPS','Shapefile','field_GIS','Bojongsoang','Bojongsoang')
JSON_FNAM = 'output.json'
OUT_FNAM = os.path.join('.','transplanting_date')

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('-s','--tmin',default=TMIN,help='Min date of transplanting in the format YYYYMMDD (%default)')
parser.add_option('-e','--tmax',default=TMAX,help='Max date of transplanting in the format YYYYMMDD (%default)')
parser.add_option('--data_tmin',default=None,help='Min date of input data in the format YYYYMMDD (%default)')
parser.add_option('--data_tmax',default=None,help='Max date of input data in the format YYYYMMDD (%default)')
parser.add_option('--tmgn',default=TMGN,type='float',help='Margin of input data in day (%default)')
parser.add_option('--tstp',default=TSTP,type='float',help='Precision of transplanting date in day (%default)')
parser.add_option('-S','--smooth',default=SMOOTH,type='float',help='Smoothing factor from 0 to 1 (%default)')
parser.add_option('--signal_v',default=SIGNAL_V,type='float',help='Max offset-free VH_BSC difference (%default)')
parser.add_option('--signal_fact',default=SIGNAL_FACT,type='float',help='Offset factor for VH_BSC difference (%default)')
parser.add_option('--offset',default=OFFSET,type='float',help='Offset of transplanting date in day (%default)')
parser.add_option('--output_epsg',default=None,type='int',help='Output EPSG (guessed from input data)')
parser.add_option('-I','--incidence_angle',default=None,help='Incidence angle file, format: date(%Y%m%d) angle(deg) (%default)')
parser.add_option('-l','--incidence_list',default=None,help='Incidence angle list, format: flag(0|1=baseline) pol(VH|VV) angle(deg) filename (%default)')
parser.add_option('-D','--datdir',default=DATDIR,help='Input data directory, not used if input_fnam is given (%default)')
parser.add_option('--search_key',default=SEARCH_KEY,help='Search key for input data, not used if input_fnam is given (%default)')
parser.add_option('--x_profile',default=X_PROFILE,help='Time of BSC profile to calculate post-transplanting signal (%default)')
parser.add_option('--y_profile',default=Y_PROFILE,help='BSC profile to calculate post-transplanting signal (%default)')
parser.add_option('--area_fnam',default=AREA_FNAM,help='Pixel area file name (%default)')
parser.add_option('--shp_fnam',default=SHP_FNAM,help='Input shapefile name (%default)')
parser.add_option('-i','--inp_fnam',default=None,help='Input GeoTIFF name (%default)')
parser.add_option('-j','--json_fnam',default=JSON_FNAM,help='Output JSON name (%default)')
parser.add_option('-o','--out_fnam',default=OUT_FNAM,help='Output shapefile name (%default)')
parser.add_option('--npy_fnam',default=None,help='Output npy file name (%default)')
parser.add_option('-F','--fig_fnam',default=None,help='Output figure name (%default)')
parser.add_option('--sort_post_s',default=False,action='store_true',help='Sort by post-transplanting signal (%default)')
parser.add_option('--early',default=False,action='store_true',help='Early estimation mode (%default)')
parser.add_option('--debug',default=False,action='store_true',help='Debug mode (%default)')
(opts,args) = parser.parse_args()

data_info = OrderedDict()
data_info['tmin'] = opts.tmin
data_info['tmax'] = opts.tmax
data_info['data_tmin'] = opts.data_tmin
data_info['data_tmax'] = opts.data_tmax
data_info['tmgn'] = opts.tmgn
data_info['tstp'] = opts.tstp
data_info['smooth'] = opts.smooth
data_info['signal_v'] = opts.signal_v
data_info['signal_fact'] = opts.signal_fact
data_info['offset'] = opts.offset
data_info['early'] = opts.early
data_info['incidence_angle'] = opts.incidence_angle if opts.incidence_angle is None else os.path.basename(opts.incidence_angle)
data_info['incidence_list'] = opts.incidence_list if opts.incidence_list is None else os.path.basename(opts.incidence_list)
data_info['x_profile'] = os.path.basename(opts.x_profile)
data_info['y_profile'] = os.path.basename(opts.y_profile)
data_info['area_fnam'] = os.path.basename(opts.area_fnam)

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
    iangle = {}
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
        meta = ds.GetMetadata()
        if 'iangle' in meta:
            iangle.update({dstr:float(meta['iangle'])})
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
    if opts.incidence_angle is not None:
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
    else:
        for dstr in iangle:
            ang = iangle[dstr]
            dif = np.abs(incidence_angle-ang)
            indx = np.argmin(dif)
            if dif[indx] > 0.1:
                raise ValueError('Error, dif={} ({})'.format(dif[indx],dstr))
            incidence_indx.update({dstr:indx})

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
        if not dstr in incidence_indx:
            raise ValueError('Error, no incidence angle for '+dstr)
        elif not incidence_select[incidence_indx[dstr]]:
            continue
        dtmp = (data[i]+signal_dif[incidence_indx[dstr]]).flatten()
    else:
        dtmp = data[i].flatten()
    data_avg = []
    for j in range(nobject):
        if inds[j].size < 1:
            data_avg.append(np.nan)
            continue
        data_value = dtmp[inds[j]]
        data_weight = areas[j]
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
data_info['dtim'] = ','.join([d.strftime('%Y%m%d') for d in vh_dtim])
with open(opts.json_fnam,'w') as json_file:
    json.dump(data_info,json_file,indent=4)
    json_file.write('\n') # Add newline cause PyJSON does not

xx = np.arange(np.floor(vh_ntim.min()),np.ceil(vh_ntim.max())+0.1*opts.tstp,opts.tstp)
cnd = (xx >= nmin) & (xx <= nmax)
xc_indx = np.where(cnd)[0]
if xc_indx.size < 3:
    raise ValueError('Error, not enough data, xc_indx.size={}'.format(xc_indx.size))
xc_indx0 = xc_indx[0]
xc_indx1 = xc_indx[1]
xc_indx_1 = xc_indx[-1]
xc_indx_2 = xc_indx[-2]
x_profile = np.load(opts.x_profile)
y_profile = np.load(opts.y_profile)
xstp = np.diff(x_profile).mean()
if np.abs(xstp-opts.tstp) > 1.0e-8:
    raise ValueError('Error, xstp={}, opts.tstp={}'.format(xstp,opts.tstp))

if opts.fig_fnam is not None:
    values = []
    labels = []
    ticks = []
    for y in range(vh_dtim[0].year,vh_dtim[-1].year+1):
        for m in range(1,13,2):
            d = datetime(y,m,1)
            values.append(date2num(d))
            labels.append(d.strftime('%Y-%m'))
        for m in range(1,13,1):
            d = datetime(y,m,1)
            ticks.append(date2num(d))
    plt.interactive(True)
    fig = plt.figure(1,facecolor='w')
    plt.subplots_adjust(left=0.12,right=0.88,bottom=0.12,top=0.80,hspace=0.5)
    pdf = PdfPages(opts.fig_fnam)

nb = 17 # (trans_dN,bsc_minN,fp_offsN,post_sN,fpi_N)x3,fpi_s,fpi_e
output_data = np.full((nb,nobject),np.nan)
if not opts.debug:
    warnings.simplefilter('ignore')
#for ii in range(nobject):
for ii in [1000]:
    if ii%1000 == 0:
        sys.stderr.write('{}/{}\n'.format(ii,nobject))
    object_id = object_ids[ii]
    if object_id != ii+1:
        raise ValueError('Error, object_id={}, ii={}'.format(object_id,ii))
    if inds[ii].size < 1:
        continue
    yi = vh_data[:,ii] # VH
    yy = csaps(vh_ntim,yi,xx,smooth=opts.smooth)

    # Calculate fishpond index
    yy_max = yy.max()
    yy_min = yy.min()
    yy_thr = max(yy_min+0.3*(yy_max-yy_min),-21.0)
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
    ff = np.zeros(xx.size)
    for f in fp_inds:
        j1 = f[0]
        j2 = f[1]
        ff[j1:j2] = (yy_thr-yy[j1:j2]).sum()
    fishpond_index = (ff/2000.0).clip(0.0,1.0)

    # Calculate post-transplanting signal
    ss = []
    oo = []
    for i in range(xx.size):
        dmax = yy.max()-yy[i]
        if (dmax > opts.signal_v) and (fishpond_index[i] > 0.01):
            offset = (dmax-opts.signal_v)*opts.signal_fact
        else:
            offset = 0.0
        ytmp = yy-(yy[i]+offset)
        i2 = min(i+y_profile.size,yy.size)
        if i2 > i:
            yfact = y_profile[0:i2-i].copy()
            yfact *= yfact.size/yfact.sum()
            ytmp[i:i2] *= yfact
        cnd = (xx > xx[i]) & (xx < xx[i]+90)
        if cnd.sum() < 1:
            y2 = 0.0
        else:
            y2 = ytmp[cnd].mean()
        ss.append(y2)
        oo.append(offset)
    ss = np.array(ss)
    oo = np.array(oo)
    vv = yy+oo

    # Search local minima and rising points
    f1 = np.gradient(yy)/np.gradient(xx)
    f2 = np.gradient(f1)/np.gradient(xx)
    cc = (f1 >= -0.1) & (f1 <= 0.1) & (f2 > 0.0)
    xcs = []
    ycs = []
    f1cs = []
    f2cs = []
    fflg = [] # True: local minimum, False: rising point
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

    # List candidates
    xest = []
    yest = []
    xflg = [] # True: latest data point
    yflg = [] # True: valley after a peak
    for ic in range(fflg.size):
        if fflg[ic]:
            indc = np.argmin(np.abs(f1cs[ic]))
        else:
            indc = np.argmax(f2cs[ic])
        xest.append(xcs[ic][indc])
        yest.append(ycs[ic][indc])
        xflg.append(False)
        yflg.append(False)
    xest = np.array(xest)
    yest = np.array(yest)
    xflg = np.array(xflg)
    yflg = np.array(yflg)
    if opts.early and np.abs(xest-xx[xc_indx_1]).min() > 1.0:
        if (yy[xc_indx_1] < yy[xc_indx_2]):
            xest = np.append(xest,xx[xc_indx_1])
            yest = np.append(yest,yy[xc_indx_1])
            xflg = np.append(xflg,True)
            yflg = np.append(yflg,False)
            fflg = np.append(fflg,False)
    test = np.array([xest[ic] if xflg[ic] else xest[ic]+opts.offset for ic in range(fflg.size)])

    # Examine the superiority of the candidates
    ydif = []
    yinc = []
    yrgt = []
    for ic in range(fflg.size):
        cnd = (xx > xest[ic]) & (xx < xest[ic]+45)
        if cnd.sum() < 1:
            yrgt.append(0.0)
        else:
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
            if cnd.sum() < 1:
                yd = 0.0
            else:
                yd = yy[cnd].max()-yest[ic]
            ydif.append(yd)
            yinc.append(yd)
    ydif = np.array(ydif)
    yinc = np.array(yinc)
    yrgt = np.array(yrgt)
    yval = []
    sval = []
    for ic in range(fflg.size):
        ix = np.argmin(np.abs(xx-xest[ic]))
        if xflg[ic]:
            y = -vv[ix]
            s = ss[ix]
        elif yflg[ic]:
            y = -vv[ix]
            s = ss[ix]
        elif yrgt[ic] < oo[ix]:
            y = -1.0e10
            s = -1.0e10
        elif yrgt[ic] < 1.5:
            y = -1.0e10
            s = -1.0e10
        elif (ic < fflg.size-1) and (ydif[ic]-yinc[ic] > 4.0) and (xest[ic+1]-xest[ic] < 60.0):
            y = -1.0e10
            s = -1.0e10
        else:
            y = -vv[ix]
            s = ss[ix]
        if (test[ic] < nmin) or (test[ic] > nmax):
            y = -1.5e10
            s = -1.5e10
        yval.append(y)
        sval.append(s)
    yval = np.array(yval)
    sval = np.array(sval)

    # Select three candidates
    if opts.sort_post_s:
        cval = sval.copy()
    else:
        cval = yval.copy()
    indv = np.argsort(cval)[::-1]
    if cval[indv[0]] > -10.0:
        ix = np.argmin(np.abs(xx-xest[indv[0]]))
        output_data[0,ii] = xx[ix]
        output_data[1,ii] = yy[ix]
        output_data[2,ii] = oo[ix]
        output_data[3,ii] = ss[ix]
        output_data[4,ii] = fishpond_index[ix]
        if not np.all(output_data[[0,1],ii] == np.array([xest[indv[0]],yest[indv[0]]])):
            raise ValueError('Error in result check 1')
        cnd = np.abs(xest-xest[indv[0]]) < 1.0
        cval[cnd] = -2.0e10
        indv = np.argsort(cval)[::-1]
    if cval[indv[0]] > -10.0:
        ix = np.argmin(np.abs(xx-xest[indv[0]]))
        output_data[5,ii] = xx[ix]
        output_data[6,ii] = yy[ix]
        output_data[7,ii] = oo[ix]
        output_data[8,ii] = ss[ix]
        output_data[9,ii] = fishpond_index[ix]
        if not np.all(output_data[[5,6],ii] == np.array([xest[indv[0]],yest[indv[0]]])):
            raise ValueError('Error in result check 2')
        cnd = np.abs(xest-xest[indv[0]]) < 1.0
        cval[cnd] = -3.0e10
        indv = np.argsort(cval)[::-1]
    if cval[indv[0]] > -10.0:
        ix = np.argmin(np.abs(xx-xest[indv[0]]))
        output_data[10,ii] = xx[ix]
        output_data[11,ii] = yy[ix]
        output_data[12,ii] = oo[ix]
        output_data[13,ii] = ss[ix]
        output_data[14,ii] = fishpond_index[ix]
        if not np.all(output_data[[10,11],ii] == np.array([xest[indv[0]],yest[indv[0]]])):
            raise ValueError('Error in result check 3')
    ix = np.argmin(np.abs(xx-nmin))
    output_data[15,ii] = fishpond_index[ix]
    ix = np.argmin(np.abs(xx-nmax))
    output_data[16,ii] = fishpond_index[ix]

    # Plot data
    if opts.fig_fnam is not None:
        ss_offset = 15.0
        bsc_min = -25.0
        bsc_max = -5.0
        fig.clear()
        ax1 = plt.subplot(111)
        ax1.minorticks_on()
        ax2 = ax1.twinx()
        if not np.all(np.abs(oo) < 1.0e-8):
            ax1.plot(xx,vv,'k--')
        l1, = ax1.plot(xx,yy,'k-',label='BSC')
        #ss[-10:] = np.nan
        l2, = ax1.plot(xx,ss-ss_offset,'y-',lw=1,label='Signal')
        l3, = ax2.plot(xx,fishpond_index,'-',color='#cccccc',label='FI',zorder=0)
        for ic in range(fflg.size):
            if fflg[ic]:
                ax1.plot(xest[ic],yest[ic],'ro')
            else:
                ax1.plot(xest[ic],yest[ic],'bo')
            if sval[ic] < -10.0:
                ax1.plot(xest[ic],yest[ic],'rx',ms=20)
        ax1.plot(output_data[10,ii],output_data[11,ii],'o',ms=20,mfc='none',color='orange',mew=2,zorder=9)
        ax1.plot(output_data[5,ii],output_data[6,ii],'o',ms=20,mfc='none',color='m',mew=2,zorder=9)
        ax1.plot(output_data[0,ii],output_data[1,ii],'o',ms=20,mfc='none',color='r',mew=2,zorder=9)
        l4 = ax1.vlines(output_data[10,ii],bsc_min,output_data[11,ii],color='orange',label='T$_{est3}$',zorder=10)
        l5 = ax1.vlines(output_data[5,ii],bsc_min,output_data[6,ii],color='m',label='T$_{est2}$',zorder=10)
        l6 = ax1.vlines(output_data[0,ii],bsc_min,output_data[1,ii],color='r',label='T$_{est1}$',zorder=10)
        ax1.fill_betweenx([bsc_min,bsc_max],xx.min(),nmin,color='k',alpha=0.1)
        ax1.fill_betweenx([bsc_min,bsc_max],nmax,xx.max(),color='k',alpha=0.1)
        ax1.set_title('OBJECTID: {}'.format(object_id),y=1.15,x=0.5)
        ax1.set_ylim(bsc_min,bsc_max)
        ax2.set_ylim(-0.1,1.1)
        ax1.set_xlim(xx.min(),xx.max())
        ax1.xaxis.set_major_locator(plt.FixedLocator(values))
        ax1.xaxis.set_major_formatter(plt.FixedFormatter(labels))
        ax1.xaxis.set_minor_locator(plt.FixedLocator(ticks))
        ax1.yaxis.set_major_locator(plt.MultipleLocator(5.0))
        ax1.set_ylabel('BSC, Signal $-$ {} (dB)'.format(int(ss_offset+0.1)))
        ax2.set_ylabel('FI')
        ax1.yaxis.set_label_coords(-0.105,0.5)
        ax2.yaxis.set_label_coords(1.10,0.5)
        for l in ax1.xaxis.get_ticklabels():
            l.set_rotation(30.0)
        lns = [l1,l2,l3,l4,l5,l6]
        lbs = [l.get_label() for l in lns]
        ax1.legend(lns,lbs,prop={'size':12},numpoints=1,loc=8,bbox_to_anchor=(0.5,1.01),ncol=9,frameon=False,handletextpad=0.1,columnspacing=0.60,handlelength=1.2)
        plt.savefig(pdf,format='pdf')
        plt.draw()
        #plt.pause(0.1)
    #break
if opts.fig_fnam is not None:
    pdf.close()

if opts.npy_fnam is not None:
    np.save(opts.npy_fnam,output_data)

r = shapefile.Reader(opts.shp_fnam)
if len(r) != nobject:
    raise ValueError('Error, len(r)={}, nobject={}'.format(len(r),nobject))
w = shapefile.Writer(opts.out_fnam)
w.shapeType = shapefile.POLYGON
w.fields = r.fields[1:] # skip first deletion field
for i in range(3):
    w.field('trans_d{}'.format(i+1),'F',13,6)
    w.field('bsc_min{}'.format(i+1),'F',13,6)
    w.field('fp_offs{}'.format(i+1),'F',13,6)
    w.field('post_s{}'.format(i+1),'F',13,6)
    w.field('fpi_{}'.format(i+1),'F',13,6)
w.field('fpi_s','F',13,6)
w.field('fpi_e','F',13,6)
for ii,shaperec in enumerate(r.iterShapeRecords()):
    rec = shaperec.record
    shp = shaperec.shape
    rec.extend(list(output_data[:,ii]))
    w.shape(shp)
    w.record(*rec)
w.close()
shutil.copy2(opts.shp_fnam+'.prj',opts.out_fnam+'.prj')
