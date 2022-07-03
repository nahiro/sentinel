#!/usr/bin/env python
import os
import sys
import re
import warnings
try:
    import gdal
except Exception:
    from osgeo import gdal
from datetime import datetime
import numpy as np
#import matplotlib
#matplotlib.use('Agg')
from matplotlib.dates import date2num,num2date
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.offsetbox import AnchoredText
from matplotlib.backends.backend_pdf import PdfPages
from optparse import OptionParser,IndentedHelpFormatter

# Default values
BAND = '4'
BAND_COL = 1
RTHR = 1.0
MTHR = 2.0
BAND4_MAX = 0.35
INDS_FNAM = 'nearest_inds.npy'

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.set_usage('Usage: %prog input_fnam [options]')
parser.add_option('-b','--band',default=BAND,help='Target band (%default)')
parser.add_option('-B','--band_fnam',default=None,help='Band file name (%default)')
parser.add_option('--band_col',default=BAND_COL,help='Band column number (%default)')
parser.add_option('-v','--vthr',default=None,type='float',help='Absolute threshold to remove outliers (%default)')
parser.add_option('-r','--rthr',default=RTHR,type='float',help='Relative threshold for 2-step outlier removal (%default)')
parser.add_option('--mthr',default=MTHR,type='float',help='Multiplying factor of vthr for 2-step outlier removal (%default)')
parser.add_option('--band4_max',default=BAND4_MAX,type='float',help='Band4 threshold (%default)')
parser.add_option('--ax1_xmin',default=None,type='float',help='Axis1 X min (%default)')
parser.add_option('--ax1_xmax',default=None,type='float',help='Axis1 X max (%default)')
parser.add_option('--ax1_ymin',default=None,type='float',help='Axis1 Y min (%default)')
parser.add_option('--ax1_ymax',default=None,type='float',help='Axis1 Y max (%default)')
parser.add_option('--ax1_zmin',default=None,type='float',help='Axis1 Z min (%default)')
parser.add_option('--ax1_zmax',default=None,type='float',help='Axis1 Z max (%default)')
parser.add_option('--mask_fnam',default=None,help='Mask file name (%default)')
parser.add_option('--stat_fnam',default=None,help='Statistic file name (%default)')
parser.add_option('--inds_fnam',default=INDS_FNAM,help='Index file name (%default)')
parser.add_option('-F','--fig_fnam',default=None,help='Output figure name for debug (%default)')
parser.add_option('-o','--output_fnam',default=None,help='Output NPZ name (%default)')
parser.add_option('--ignore_band4',default=False,action='store_true',help='Ignore exceeding the band4 threshold (%default)')
parser.add_option('--outlier_remove2',default=False,action='store_true',help='2-step outlier removal mode (%default)')
parser.add_option('--debug',default=False,action='store_true',help='Debug mode (%default)')
(opts,args) = parser.parse_args()
if len(args) < 1:
    parser.print_help()
    sys.exit(0)
input_fnam = args[0]
m = re.search('^('+'\d'*8+')_',os.path.basename(input_fnam))
if not m:
    raise ValueError('Error in file name >>> '+input_fnam)
dstr = m.group(1)

if opts.band.upper() == 'NDVI':
    band_l = opts.band.lower()
    band_u = opts.band.upper()
    if opts.vthr is None:
        opts.vthr = 0.1
    if opts.ax1_xmin is None:
        opts.ax1_xmin = -0.5
    if opts.ax1_xmax is None:
        opts.ax1_xmax = 1.0
    if opts.ax1_ymin is None:
        opts.ax1_ymin = -0.5
    if opts.ax1_ymax is None:
        opts.ax1_ymax = 1.0
    if opts.ax1_zmin is None:
        opts.ax1_zmin = 0.01
    if opts.ax1_zmax is None:
        opts.ax1_zmax = 0.05
else:
    band_l = 'band'+opts.band
    band_u = 'Band'+opts.band
    if opts.vthr is None:
        opts.vthr = 0.02
    if opts.ax1_xmin is None:
        opts.ax1_xmin = 0.0
    if opts.ax1_xmax is None:
        opts.ax1_xmax = 0.5
    if opts.ax1_ymin is None:
        opts.ax1_ymin = 0.0
    if opts.ax1_ymax is None:
        opts.ax1_ymax = 0.5
    if opts.ax1_zmin is None:
        opts.ax1_zmin = 0.01
    if opts.ax1_zmax is None:
        opts.ax1_zmax = 0.035
if opts.output_fnam is None:
    opts.output_fnam = 'atcor_param_{}_{}.npz'.format(band_l,dstr)
if opts.fig_fnam is None:
    opts.fig_fnam = 'sentinel2_atcor_{}_{}.pdf'.format(band_l,dstr)

stat = np.load(opts.stat_fnam)
data_y_all = stat['mean'].flatten()
data_z_all = stat['std'].flatten()
nearest_inds = np.load(opts.inds_fnam)
nobject = len(nearest_inds)

ds = gdal.Open(input_fnam)
data = ds.ReadAsArray()
data_shape = data[0].shape
band_list = []
if opts.band_fnam is not None:
    with open(opts.band_fnam,'r') as fp:
        for line in fp:
            item = line.split()
            if len(item) <= opts.band_col or item[0][0]=='#':
                continue
            band_list.append(item[opts.band_col])
    if len(data) != len(band_list):
        raise ValueError('Error, len(data)={}, len(band_list)={} >>> {}'.format(len(data),len(band_list),input_fnam))
else:
    for i in range(ds.RasterCount):
        band = ds.GetRasterBand(i+1)
        band_name = band.GetDescription()
        if not band_name:
            raise ValueError('Error, faild to read band name >>> {}'.format(input_fnam))
        band_list.append(band_name)
ds = None
if opts.band.upper() == 'NDVI':
    band_name = 'B4'
    if not band_name in band_list:
        raise ValueError('Error, faild to search index for {}'.format(band_name))
    band4_index = band_list.index(band_name)
    b4_img = data[band4_index].astype(np.float64).flatten()
    band_name = 'B8'
    if not band_name in band_list:
        raise ValueError('Error, faild to search index for {}'.format(band_name))
    band8_index = band_list.index(band_name)
    b8_img = data[band8_index].astype(np.float64).flatten()
    data_img = (b8_img-b4_img)/(b8_img+b4_img)
    b4_img *= 1.0e-4
else:
    band_name = 'B{}'.format(opts.band)
    if not band_name in band_list:
        raise ValueError('Error, faild to search index for {}'.format(band_name))
    band_index = band_list.index(band_name)
    data_img = data[band_index].astype(np.float64).flatten()*1.0e-4
    if opts.band == 4:
        b4_img = data_img.copy()
    else:
        band_name = 'B4'
        if not band_name in band_list:
            raise ValueError('Error, faild to search index for {}'.format(band_name))
        band4_index = band_list.index(band_name)
        b4_img = data[band4_index].astype(np.float64).flatten()*1.0e-4
if opts.mask_fnam is not None:
    ds = gdal.Open(opts.mask_fnam)
    mask = ds.ReadAsArray()
    ds = None
    mask_shape = mask.shape
    if data_shape != mask_shape:
        raise ValueError('Error, data_shape={}, mask_shape={}'.format(data_shape,mask_shape))
    cnd = ~np.isnan(mask.flatten())
    data_x_all = data_img.flatten()[cnd]
    data_b_all = b4_img.flatten()[cnd]
else:
    data_x_all = data_img.flatten()
    data_b_all = b4_img.flatten()
if data_x_all.shape != data_y_all.shape:
    raise ValueError('Error, data_x_all.shape={}, data_y_all.shape={}'.format(data_x_all.shape,data_y_all.shape))

if opts.debug:
    #plt.interactive(True)
    fig = plt.figure(1,facecolor='w',figsize=(6,5))
    plt.subplots_adjust(top=0.85,bottom=0.20,left=0.15,right=0.90)
    pdf = PdfPages(opts.fig_fnam)
    xfit = np.arange(-1.0,1.001,0.01)
else:
    warnings.simplefilter('ignore')

number = []
corcoef = []
b4_mean = []
b4_std = []
factor = []
offset = []
rmse = []
for i in range(nobject):
    if opts.debug and (i%500 == 0):
        sys.stderr.write('{}\n'.format(i))
    object_id = i+1
    indx = nearest_inds[i]
    data_x = data_x_all[indx]
    data_y = data_y_all[indx]
    data_z = data_z_all[indx]
    data_b = data_b_all[indx]
    if opts.outlier_remove2:
        cnd1 = (~np.isnan(data_x)) & (np.abs(data_x-data_y) < (np.abs(data_y)*opts.rthr).clip(min=opts.vthr*opts.mthr))
    else:
        cnd1 = ~np.isnan(data_x)
    if not opts.ignore_band4:
        cnd1 &= ((~np.isnan(data_b)) & (data_b < opts.band4_max))
    xcnd = data_x[cnd1]
    ycnd = data_y[cnd1]
    zcnd = data_z[cnd1]
    bcnd = data_b[cnd1]
    if xcnd.size < 2:
        result = [np.nan,np.nan]
        calc_y = np.array([])
        r_value = np.nan
        b4_value = np.nan
        b4_error = np.nan
        rms_value = np.nan
        flag = False
    else:
        result = np.polyfit(xcnd,ycnd,1)
        calc_y = xcnd*result[0]+result[1]
        cnd2 = np.abs(calc_y-ycnd) < opts.vthr
        xcnd2 = xcnd[cnd2]
        ycnd2 = ycnd[cnd2]
        zcnd2 = zcnd[cnd2]
        bcnd2 = bcnd[cnd2]
        if (xcnd2.size == cnd2.size) or (xcnd2.size < 2):
            r_value = np.corrcoef(xcnd,ycnd)[0,1]
            b4_value = np.nanmean(bcnd)
            b4_error = np.nanstd(bcnd)
            rms_value = np.sqrt(np.square(calc_y-ycnd).sum()/calc_y.size)
            flag = False
        else:
            result = np.polyfit(xcnd2,ycnd2,1)
            calc_y = xcnd2*result[0]+result[1]
            r_value = np.corrcoef(xcnd2,ycnd2)[0,1]
            b4_value = np.nanmean(bcnd2)
            b4_error = np.nanstd(bcnd2)
            rms_value = np.sqrt(np.square(calc_y-ycnd2).sum()/calc_y.size)
            flag = True
    number.append(calc_y.size)
    corcoef.append(r_value)
    b4_mean.append(b4_value)
    b4_std.append(b4_error)
    factor.append(result[0])
    offset.append(result[1])
    rmse.append(rms_value)
    if opts.debug:
        fig.clear()
        ax1 = plt.subplot(111)#,aspect='equal')
        ax1.minorticks_on()
        ax1.grid(True)
        line = 'number: {}\n'.format(calc_y.size)
        line += 'coeff: {:7.4f}\n'.format(r_value)
        line += 'band4: {:6.4f}\n'.format(b4_value)
        line += 'factor: {:5.4f}\n'.format(result[0])
        line += 'offset: {:5.4f}\n'.format(result[1])
        line += 'rmse: {:5.4f}'.format(rms_value)
        at= AnchoredText(line,prop=dict(size=12),loc=2,pad=0.1,borderpad=0.8,frameon=True)
        at.patch.set_boxstyle('round')
        ax1.add_artist(at)
        if xcnd.size != cnd1.size:
            ax1.scatter(data_x,data_y,c='k',marker='.')
        if flag:
            ax1.scatter(xcnd,ycnd,c='#888888',marker='.')
            im = ax1.scatter(xcnd2,ycnd2,c=zcnd2,marker='.',cmap=cm.jet,vmin=opts.ax1_zmin,vmax=opts.ax1_zmax)
        else:
            im = ax1.scatter(xcnd,ycnd,c=zcnd,marker='.',cmap=cm.jet,vmin=opts.ax1_zmin,vmax=opts.ax1_zmax)
        ax1.plot(xfit,np.polyval(result,xfit),'k:')
        ax2 = plt.colorbar(im,ticks=np.arange(0.0,0.101,0.01)).ax
        ax2.minorticks_on()
        ax2.set_ylabel(band_u+' std')
        ax1.set_xlim(opts.ax1_xmin,opts.ax1_xmax)
        ax1.set_ylim(opts.ax1_ymin,opts.ax1_ymax)
        ax1.set_xlabel(band_u)
        ax1.set_ylabel(band_u+' mean')
        ax1.xaxis.set_tick_params(pad=7)
        ax1.xaxis.set_label_coords(0.5,-0.14)
        ax1.yaxis.set_label_coords(-0.15,0.5)
        ax2.yaxis.set_label_coords(5.0,0.5)
        ax1.set_title('{} (OBJECTID={})'.format(dstr,object_id))
        plt.savefig(pdf,format='pdf')
        #plt.draw()
        #plt.pause(0.1)
        #break
if opts.debug:
    pdf.close()
number = np.array(number)
corcoef = np.array(corcoef)
b4_mean = np.array(b4_mean)
b4_std = np.array(b4_std)
factor = np.array(factor)
offset = np.array(offset)
rmse = np.array(rmse)
np.savez(opts.output_fnam,number=number,corcoef=corcoef,b4_mean=b4_mean,b4_std=b4_std,factor=factor,offset=offset,rmse=rmse)
