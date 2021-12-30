#!/usr/bin/env python
import gdal
import osr
import numpy as np
from scipy.interpolate import griddata
from optparse import OptionParser,IndentedHelpFormatter

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('-R','--ref_geotiff',default=None,help='Reference GeoTIFF name (%default)')
parser.add_option('-I','--src_geotiff',default=None,help='Source GeoTIFF name (%default)')
parser.add_option('-O','--dst_geotiff',default=None,help='Destination GeoTIFF name (%default)')
(opts,args) = parser.parse_args()

ds = gdal.Open(opts.ref_geotiff)
ref_nx = ds.RasterXSize
ref_ny = ds.RasterYSize
ref_shape = (ref_ny,ref_nx)
ref_prj = ds.GetProjection()
ref_trans = ds.GetGeoTransform()
ref_indy,ref_indx = np.indices(ref_shape)
ref_xp = ref_trans[0]+(ref_indx+0.5)*ref_trans[1]+(ref_indy+0.5)*ref_trans[2]
ref_yp = ref_trans[3]+(ref_indx+0.5)*ref_trans[4]+(ref_indy+0.5)*ref_trans[5]
ref_xp = ref_xp.flatten()
ref_yp = ref_yp.flatten()
ds = None

ds = gdal.Open(opts.src_geotiff)
src_nx = ds.RasterXSize
src_ny = ds.RasterYSize
src_nb = ds.RasterCount
src_shape = (src_ny,src_nx)
src_prj = ds.GetProjection()
src_trans = ds.GetGeoTransform()
src_meta = ds.GetMetadata()
src_data = ds.ReadAsArray().astype(np.float64).reshape(src_nb,src_ny,src_nx)
src_band = []
for iband in range(src_nb):
    band = ds.GetRasterBand(iband+1)
    src_band.append(band.GetDescription())
src_nodata = band.GetNoDataValue()
src_indy,src_indx = np.indices(src_shape)
src_xp = src_trans[0]+(src_indx+0.5)*src_trans[1]+(src_indy+0.5)*src_trans[2]
src_yp = src_trans[3]+(src_indx+0.5)*src_trans[4]+(src_indy+0.5)*src_trans[5]
src_xp = src_xp.flatten()
src_yp = src_yp.flatten()
ds = None

dst_nx = ref_nx
dst_ny = ref_ny
dst_nb = src_nb
dst_prj = ref_prj
dst_trans = ref_trans
dst_meta = src_meta
dst_data = []
dst_band = []
for iband in range(dst_nb):
    dst_data.append(griddata((src_xp,src_yp),src_data[iband].flatten(),(ref_xp,ref_yp),method='nearest').reshape(ref_shape))
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
