#!/usr/bin/env python
import gdal
import osr
import numpy as np
import warnings
from optparse import OptionParser,IndentedHelpFormatter

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('-I','--src_geotiff',default=None,help='Source GeoTIFF name (%default)')
parser.add_option('-O','--dst_geotiff',default=None,help='Destination GeoTIFF name (%default)')
parser.add_option('-x','--imin',default=None,type='int',help='Start x index (%default)')
parser.add_option('-X','--imax',default=None,type='int',help='Stop x index (%default)')
parser.add_option('-s','--istp',default=None,type='int',help='Step x index (%default)')
parser.add_option('-y','--jmin',default=None,type='int',help='Start y index (%default)')
parser.add_option('-Y','--jmax',default=None,type='int',help='Stop y index (%default)')
parser.add_option('-S','--jstp',default=None,type='int',help='Step y index (%default)')
(opts,args) = parser.parse_args()

ds = gdal.Open(opts.src_geotiff)
src_nx = ds.RasterXSize
src_ny = ds.RasterYSize
src_nb = ds.RasterCount
src_prj = ds.GetProjection()
src_trans = ds.GetGeoTransform()
if src_trans[2] != 0.0 or src_trans[4] != 0.0:
    raise ValueError('Error, src_trans={}'.format(src_trans))
src_meta = ds.GetMetadata()
src_data = ds.ReadAsArray().astype(np.float64)
src_band = []
for iband in range(src_nb):
    band = ds.GetRasterBand(iband+1)
    src_band.append(band.GetDescription())
src_nodata = band.GetNoDataValue()
src_xmin = src_trans[0]
src_xstp = src_trans[1]
src_ymax = src_trans[3]
src_ystp = src_trans[5]
src_xgrd = src_xmin+(np.arange(src_nx)+0.5)*src_xstp
src_ygrd = src_ymax+(np.arange(src_ny)+0.5)*src_ystp
ds = None

if opts.imin is None:
    opts.imin = 0
if opts.imax is None:
    opts.imax = src_nx
if opts.jmin is None:
    opts.jmin = 0
if opts.jmax is None:
    opts.jmax = src_ny

dst_nx = (opts.imax-opts.imin)//opts.istp
dst_ny = (opts.jmax-opts.jmin)//opts.jstp
dst_nb = src_nb
dst_prj = src_prj
dst_trans = [0.0]*len(src_trans)
dst_trans[0] = src_xgrd[opts.imin]-0.5*src_xstp
dst_trans[1] = src_xstp*opts.istp
dst_trans[3] = src_ygrd[opts.jmin]-0.5*src_ystp
dst_trans[5] = src_ystp*opts.jstp
dst_meta = src_meta
dst_data = []
dst_band = []
for iband in range(dst_nb):
    tmp_data = src_data[iband,opts.jmin:opts.jmin+opts.jstp*dst_ny,opts.imin:opts.imin+opts.istp*dst_nx].reshape(dst_ny,opts.jstp,dst_nx,opts.istp)
    if src_nodata is None:
        dst_data.append(tmp_data.mean(axis=-1).mean(axis=1))
    elif np.isnan(src_nodata):
        cnd = np.isnan(tmp_data)
        if cnd.sum() > 0:
            with warnings.catch_warnings():
                warnings.filterwarnings(action='ignore',message='Mean of empty slice')
                dst_data.append(np.nanmean(np.nanmean(tmp_data,axis=-1),axis=1))
        else:
            dst_data.append(tmp_data.mean(axis=-1).mean(axis=1))
    else:
        cnd = (tmp_data == src_nodata)
        if cnd.sum() > 0:
            tmp_data[cnd] = np.nan
            with warnings.catch_warnings():
                warnings.filterwarnings(action='ignore',message='Mean of empty slice')
                avg_data = np.nanmean(np.nanmean(tmp_data,axis=-1),axis=1)
            cnd = np.isnan(avg_data)
            avg_data[cnd] = src_nodata
            dst_data.append(avg_data)
        else:
            dst_data.append(tmp_data.mean(axis=-1).mean(axis=1))
    dst_band.append(src_band[iband])
dst_data = np.array(dst_data).astype(np.float32)
dst_nodata = src_nodata

drv = gdal.GetDriverByName('GTiff')
ds = drv.Create(opts.dst_geotiff,dst_nx,dst_ny,dst_nb,gdal.GDT_Float32)
ds.SetProjection(dst_prj)
ds.SetGeoTransform(dst_trans)
ds.SetMetadata(dst_meta)
for iband in range(dst_nb):
    band = ds.GetRasterBand(iband+1)
    band.WriteArray(dst_data[iband])
    band.SetDescription(dst_band[iband])
if dst_nodata is None:
    band.SetNoDataValue(np.nan) # The TIFFTAG_GDAL_NODATA only support one value per dataset
else:
    band.SetNoDataValue(dst_nodata) # The TIFFTAG_GDAL_NODATA only support one value per dataset
ds.FlushCache()
ds = None # close dataset
