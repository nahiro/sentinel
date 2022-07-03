#!/usr/bin/env python
import os
import sys
import re
import numpy as np
try:
    import gdal
except Exception:
    from osgeo import gdal
try:
    import osr
except Exception:
    from osgeo import osr
from optparse import OptionParser,IndentedHelpFormatter

# Default values
DIST_FNAMS = ['pixel_dist.dat']
OUT_FNAM = 'dist_map.tif'

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('-f','--img_fnam',default=None,help='Image file name (%default)')
parser.add_option('-a','--dist_fnams',default=None,action='append',help='Pixel dist file names ({})'.format(DIST_FNAMS))
parser.add_option('-o','--out_fnam',default=OUT_FNAM,help='Output GeoTIFF name (%default)')
(opts,args) = parser.parse_args()
if opts.dist_fnams is None:
    opts.dist_fnams = DIST_FNAMS

ds = gdal.Open(opts.img_fnam)
data = ds.ReadAsArray()
data_shape = (data.shape[-2],data.shape[-1])
data_trans = ds.GetGeoTransform()
prj = ds.GetProjection()
srs = osr.SpatialReference(wkt=prj)
estr = srs.GetAttrValue('AUTHORITY',1)
if re.search('\D',estr):
    raise ValueError('Error in EPSG >>> '+estr)
data_epsg = int(estr)
ds = None

dist_map = np.full(data_shape,2.0e10)

for dist_file in opts.dist_fnams:
    object_ids = []
    inds = []
    dists = []
    with open(dist_file,'r') as fp:
        for line in fp:
            item = line.split()
            if len(item) < 4 or item[0] == '#':
                continue
            object_ids.append(int(item[0]))
            inds.append([])
            dists.append([])
            n = int(item[1])
            for nn in range(2,n*2+2,2):
                inds[-1].append(int(item[nn]))
                dists[-1].append(float(item[nn+1]))
            inds[-1] = np.array(inds[-1])
            dists[-1] = np.array(dists[-1])
    object_ids = np.array(object_ids)
    inds = np.array(inds,dtype='object')
    dists = np.array(dists,dtype='object')

    for i in range(len(object_ids)):
        #object_id = i+1
        #if object_id != object_ids[i]:
        #    raise ValueError('Error, object_id={}, object_ids[{}]={}'.format(object_id,i,object_ids[i]))
        indy,indx = np.unravel_index(inds[i],data_shape)
        dist_map[indy,indx] = np.minimum(dist_map[indy,indx],dists[i])

dist_map[dist_map > 1.0e10] = np.nan

drv = gdal.GetDriverByName('GTiff')
ds = drv.Create(opts.out_fnam,data_shape[1],data_shape[0],1,gdal.GDT_Float32)
ds.SetGeoTransform(data_trans)
srs = osr.SpatialReference()
srs.ImportFromEPSG(data_epsg)
ds.SetProjection(srs.ExportToWkt())
band = ds.GetRasterBand(1)
band.WriteArray(dist_map)
band.SetDescription('dist_map')
band.SetNoDataValue(np.nan) # The TIFFTAG_GDAL_NODATA only support one value per dataset
ds.FlushCache()
ds = None # close dataset
