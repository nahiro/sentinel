#!/usr/bin/env python
import os
import sys
import warnings
import gdal
from datetime import datetime
import numpy as np
#import matplotlib
#matplotlib.use('Agg')
from matplotlib.dates import date2num,num2date
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
from matplotlib.backends.backend_pdf import PdfPages
from optparse import OptionParser,IndentedHelpFormatter

# Default values
BAND = '4'
TMIN = '20160101'
TMAX = '20300101'
RTHR = 1.0
DATDIR = '..'
INDS_FNAM = 'nearest_inds_1000.npy'

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('-b','--band',default=BAND,help='Target band (%default)')
parser.add_option('-s','--tmin',default=TMIN,help='Min date of transplanting in the format YYYYMMDD (%default)')
parser.add_option('-e','--tmax',default=TMAX,help='Max date of transplanting in the format YYYYMMDD (%default)')
parser.add_option('-r','--rthr',default=RTHR,type='float',help='Relative threshold to remove outliers (%default)')
parser.add_option('-v','--vthr',default=None,type='float',help='Absolute threshold to remove outliers (%default)')
parser.add_option('--ax1_xmin',default=None,type='float',help='Axis1 X min (%default)')
parser.add_option('--ax1_xmax',default=None,type='float',help='Axis1 X max (%default)')
parser.add_option('--ax1_ymin',default=None,type='float',help='Axis1 Y min (%default)')
parser.add_option('--ax1_ymax',default=None,type='float',help='Axis1 Y max (%default)')
parser.add_option('--ax1_zmin',default=None,type='float',help='Axis1 Z min (%default)')
parser.add_option('--ax1_zmax',default=None,type='float',help='Axis1 Z max (%default)')
parser.add_option('-D','--datdir',default=DATDIR,help='Input data directory (%default)')
parser.add_option('--inds_fnam',default=INDS_FNAM,help='Index file name (%default)')
parser.add_option('-F','--fig_fnam',default=None,help='Output figure name for debug (%default)')
parser.add_option('--debug',default=False,action='store_true',help='Debug mode (%default)')
(opts,args) = parser.parse_args()

nmin = date2num(datetime.strptime(opts.tmin,'%Y%m%d'))
nmax = date2num(datetime.strptime(opts.tmax,'%Y%m%d'))
ntim = np.load(os.path.join(opts.datdir,'ntim.npy'))

if opts.band.upper() == 'NDVI':
    band_l = opts.band.lower()
    band_u = opts.band.upper()
    band_s = 'nv'
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
    band_s = 'b'+opts.band
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
if opts.fig_fnam is None:
    opts.fig_fnam = 'sentinel2_atcor_'+band_s+'.pdf'

data_x_all = np.load(os.path.join(opts.datdir,band_l+'.npy'))
data_y_all = np.load(os.path.join(opts.datdir,band_s+'_mean.npy'))
data_z_all = np.load(os.path.join(opts.datdir,band_s+'_std.npy'))
nearest_inds = np.load(opts.inds_fnam)
nobject = len(nearest_inds)

if opts.debug:
    fig = plt.figure(1,facecolor='w',figsize=(6,5))
    plt.subplots_adjust(top=0.85,bottom=0.20,left=0.15,right=0.90)
    pdf = PdfPages(opts.fig_fnam)
    xfit = np.arange(-1.0,1.001,0.01)
else:
    warnings.simplefilter('ignore')

for indt in range(ntim.size):
    if ntim[indt] < nmin-1.0e-3 or ntim[indt] > nmax+1.0e-3:
        continue
    dtim = num2date(ntim[indt])
    dstr = dtim.strftime('%Y%m%d')
    sys.stderr.write(dstr+'\n')
    factor = []
    offset = []
    rmse = []
    for i in range(nobject):
        if opts.debug and (i%500 == 0):
            sys.stderr.write('{}\n'.format(i))
        object_id = i+1
        indx = nearest_inds[i]
        data_x = data_x_all[indt,indx]
        data_y = data_y_all[indx]
        data_z = data_z_all[indx]
        cnd1 = ~np.isnan(data_x) & (np.abs(data_x-data_y) < np.abs(data_y)*opts.rthr)
        xcnd = data_x[cnd1]
        ycnd = data_y[cnd1]
        zcnd = data_z[cnd1]
        if xcnd.size < 2:
            result = [np.nan,np.nan]
            calc_y = np.array([])
            rms_value = np.nan
            flag = False
        else:
            result = np.polyfit(xcnd,ycnd,1)
            calc_y = xcnd*result[0]+result[1]
            cnd2 = np.abs(calc_y-ycnd) < opts.vthr
            if cnd2.sum() != cnd2.size:
                flag = True
                xcnd2 = xcnd[cnd2]
                ycnd2 = ycnd[cnd2]
                zcnd2 = zcnd[cnd2]
                result = np.polyfit(xcnd2,ycnd2,1)
                calc_y = xcnd2*result[0]+result[1]
                rms_value = np.sqrt(np.square(calc_y-ycnd2).sum()/calc_y.size)
            else:
                flag = False
                rms_value = np.sqrt(np.square(calc_y-ycnd).sum()/calc_y.size)
        factor.append(result[0])
        offset.append(result[1])
        rmse.append(rms_value)
        if opts.debug:
            fig.clear()
            ax1 = plt.subplot(111,aspect='equal')
            ax1.minorticks_on()
            ax1.grid(True)
            line = 'number: {}\n'.format(calc_y.size)
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
                im = ax1.scatter(data_x,data_y,c=data_z,marker='.',cmap=cm.jet,vmin=opts.ax1_zmin,vmax=opts.ax1_zmax)
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
    factor = np.array(factor)
    offset = np.array(offset)
    rmse = np.array(rmse)
    np.save('atcor_{}_factor_{}.npy'.format(band_s,dstr),factor)
    np.save('atcor_{}_offset_{}.npy'.format(band_s,dstr),offset)
    np.save('atcor_{}_rmse_{}.npy'.format(band_s,dstr),rmse)
