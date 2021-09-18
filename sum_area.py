#!/usr/bin/env python
import os
import sys
import re
import numpy as np
import gdal
import osr
from optparse import OptionParser,IndentedHelpFormatter

# Default values
DATNAM = ['pixel_area.dat']
OUT_FNAM = 'area_map.tif'

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('-f','--img_fnam',default=None,help='Image file name (%default)')
parser.add_option('-i','--datnam',default=None,action='append',help='Input data name ({})'.format(DATNAM))
parser.add_option('-o','--out_fnam',default=OUT_FNAM,help='Output GeoTIFF name (%default)')
(opts,args) = parser.parse_args()
if opts.datnam is None:
    opts.datnam = DATNAM

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

area_map = np.zeros(data_shape)

for area_file in opts.datnam:
    object_ids = []
    inds = []
    areas = []
    with open(area_file,'r') as fp:
        for line in fp:
            item = line.split()
            if len(item) < 4 or item[0] == '#':
                continue
            object_ids.append(int(item[0]))
            inds.append([])
            areas.append([])
            n = int(item[1])
            for nn in range(2,n*2+2,2):
                inds[-1].append(int(item[nn]))
                areas[-1].append(float(item[nn+1]))
            inds[-1] = np.array(inds[-1])
            areas[-1] = np.array(areas[-1])
    object_ids = np.array(object_ids)
    inds = np.array(inds)
    areas = np.array(areas)

    for i in range(len(object_ids)):
        #object_id = i+1
        #if object_id != object_ids[i]:
        #    raise ValueError('Error, object_id={}, object_ids[{}]={}'.format(object_id,i,object_ids[i]))
        indy,indx = np.unravel_index(inds[i],data_shape)
        area_map[indy,indx] += areas[i]

area_map[area_map < 1.0e-16] = np.nan

drv = gdal.GetDriverByName('GTiff')
ds = drv.Create(opts.out_fnam,data_shape[1],data_shape[0],1,gdal.GDT_Float32)
ds.SetGeoTransform(data_trans)
srs = osr.SpatialReference()
srs.ImportFromEPSG(data_epsg)
ds.SetProjection(srs.ExportToWkt())
band = ds.GetRasterBand(1)
band.WriteArray(area_map)
band.SetDescription('area_map')
band.SetNoDataValue(np.nan) # The TIFFTAG_GDAL_NODATA only support one value per dataset
ds.FlushCache()
ds = None # close dataset
