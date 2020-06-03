#!/usr/bin/env python
import os
import sys
import re
import numpy as np
import tifffile
import xml.etree.ElementTree as ET
import gdal
import osr
from scipy.interpolate import griddata
from optparse import OptionParser,IndentedHelpFormatter

# Defaults
XMIN = 743800.0
XMAX = 756800.0
YMIN = 9236000.0
YMAX = 9251800.0
XSTP = 10.0
YSTP = -10.0
BAND_COL = 1

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.set_usage('Usage: %prog list_of_input_file [options]')
parser.add_option('-b','--output_band',default=None,action='append',help='Output band name (%default)')
parser.add_option('--output_bmin',default=None,type='int',help='Minimum output band index (%default)')
parser.add_option('--output_bmax',default=None,type='int',help='Maximum output band index (%default)')
parser.add_option('-B','--band_fnam',default=None,help='Band file name (%default)')
parser.add_option('-e','--output_epsg',default=None,type='int',help='Output EPSG (guessed from input data)')
parser.add_option('-x','--xmin',default=XMIN,type='float',help='Minimum X in m (%default)')
parser.add_option('-X','--xmax',default=XMAX,type='float',help='Maximum X in m (%default)')
parser.add_option('--xstp',default=XSTP,type='float',help='Step X in m (%default)')
parser.add_option('-y','--ymin',default=YMIN,type='float',help='Minimum Y in m (%default)')
parser.add_option('-Y','--ymax',default=YMAX,type='float',help='Maximum Y in m (%default)')
parser.add_option('--ystp',default=YSTP,type='float',help='Step Y in m (%default)')
parser.add_option('--band_col',default=BAND_COL,help='Band column number (%default)')
parser.add_option('--check_grid',default=False,action='store_true',help='Check grid (%default)')
(opts,args) = parser.parse_args()
if len(args) < 1:
    parser.print_help()
    sys.exit(0)
fnams = args

xg,yg = np.meshgrid(np.arange(opts.xmin,opts.xmax+0.1*opts.xstp,opts.xstp),np.arange(opts.ymax,opts.ymin-0.1*opts.ystp,opts.ystp))
ngrd = xg.size
ny,nx = xg.shape

for input_fnam in fnams:
    f,e = os.path.splitext(os.path.basename(input_fnam))
    output_fnam = f+'_resample'+e
    sys.stderr.write('input: '+input_fnam+', output: '+output_fnam+'\n')
    ds = gdal.Open(input_fnam)
    prj = ds.GetProjection()
    srs = osr.SpatialReference(wkt=prj)
    if opts.output_epsg is None:
        epsg = srs.GetAttrValue('AUTHORITY',1)
        if re.search('\D',epsg):
            raise ValueError('Error in EPSG >>> '+epsg)
        output_epsg = int(epsg)
    else:
        output_epsg = opts.output_epsg
    data = ds.ReadAsArray()
    trans = ds.GetGeoTransform() # maybe obtained from tif_tags['ModelTransformationTag']
    indy,indx = np.indices(data[0].shape)
    xp = trans[0]+(indx+0.5)*trans[1]+(indy+0.5)*trans[2]
    yp = trans[3]+(indx+0.5)*trans[4]+(indy+0.5)*trans[5]
    ndat = len(data)
    if opts.check_grid:
        indx1 = np.argmin(np.abs(xp[0,:]-xg[0,0]))
        indx2 = np.argmin(np.abs(xp[0,:]-xg[0,-1]))+1
        indy1 = np.argmin(np.abs(yp[:,0]-yg[0,0]))
        indy2 = np.argmin(np.abs(yp[:,0]-yg[-1,0]))+1
        if np.all(xg[0,:] == xp[0,indx1:indx2]) and np.all(yg[:,0] == yp[indy1:indy2,0]):
            flag_grid = True
        else:
            flag_grid = False
    else:
        flag_grid = False

    # Get band name
    band_name = []
    if opts.band_fnam is not None:
        with open(opts.band_fnam,'r') as fp:
            for line in fp:
                item = line.split()
                if len(item) <= opts.band_col or item[0][0]=='#':
                    continue
                band_name.append(item[opts.band_col])
    else:
        if ds.GetRasterBand(1).GetDescription() != '':
            for i in range(ndat):
                band = ds.GetRasterBand(i+1)
                band_name.append(band.GetDescription())
        else:
            tif_tags = {}
            with tifffile.TiffFile(input_fnam) as tif:
                for tag in tif.pages[0].tags.values():
                    name,value = tag.name,tag.value
                    tif_tags[name] = value
            if '65000' in tif_tags:
                root = ET.fromstring(tif_tags['65000'])
                for value in root.iter('BAND_NAME'):
                    band_name.append(value.text)
            else:
                for i in range(ndat):
                    band_name.append('band_{}'.format(i))
    nband = len(band_name)
    if nband != ndat:
        raise ValueError('Error, nband={}, ndat={}'.format(nband,ndat))
    ds = None # close dataset

    if opts.output_bmin is not None:
        if opts.output_bmax is not None:
            indxs = list(range(opts.output_bmin,opts.output_bmax+1))
        else:
            indxs = list(range(opts.output_bmin,ndat))
    elif opts.output_bmax is not None:
        indxs = list(range(0,opts.output_bmax+1))
    elif opts.output_band is None:
        indxs = list(range(0,ndat))
    else:
        indxs = []
    if opts.output_band is not None:
        for band in opts.output_band:
            indxs.append(band_name.index(band))
    nset = len(indxs)
    for i in indxs:
        sys.stderr.write('{}\n'.format(band_name[i]))
    if flag_grid:
        dset = data[indxs,indy1:indy2,indx1:indx2]
    else:
        dset = []
        for i in indxs:
            dset.append(griddata((xp.flatten(),yp.flatten()),data[i].flatten(),(xg.flatten(),yg.flatten()),method='nearest').reshape(xg.shape))
        dset = np.array(dset)

    drv = gdal.GetDriverByName('GTiff')
    ds = drv.Create(output_fnam,nx,ny,nset,gdal.GDT_Float32)
    ds.SetGeoTransform((opts.xmin-0.5*opts.xstp,opts.xstp,0.0,opts.ymax-0.5*opts.ystp,0.0,opts.ystp))
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(output_epsg)
    ds.SetProjection(srs.ExportToWkt())
    for i in range(nset):
        band = ds.GetRasterBand(i+1)
        band.WriteArray(dset[i])
        band.SetDescription(band_name[indxs[i]])
    band.SetNoDataValue(np.nan) # The TIFFTAG_GDAL_NODATA only support one value per dataset
    ds.FlushCache()
    ds = None # close dataset
