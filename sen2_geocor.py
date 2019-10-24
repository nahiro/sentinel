#!/usr/bin/env python
import os
import sys
import re
import gdal
import osr
#import tifffile
#import xml.etree.ElementTree as ET
import numpy as np
from scipy.interpolate import interp2d
from scipy.optimize import leastsq
from optparse import OptionParser,IndentedHelpFormatter

# Default values
SUBSET_WIDTH = 100 # pixel
SUBSET_HEIGHT = 100 # pixel
SHIFT_WIDTH = 3 # pixel
SHIFT_HEIGHT = 3 # pixel
MARGIN_WIDTH = 5 # pixel
MARGIN_HEIGHT = 5 # pixel
REF_BAND = 6
TRG_BAND = 7
FEPS = 0.01
RTHR = 0.3

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.set_usage('Usage: %prog reference_geotiff_file target_geotiff_file [options]')
parser.add_option('-b','--ref_band',default=REF_BAND,type='int',help='Reference band# (%default)')
parser.add_option('-B','--trg_band',default=TRG_BAND,type='int',help='Target band# (%default)')
parser.add_option('-W','--subset_width',default=SUBSET_WIDTH,type='int',help='Subset width in target pixel (%default)')
parser.add_option('-H','--subset_height',default=SUBSET_HEIGHT,type='int',help='Subset height in target pixel (%default)')
parser.add_option('--shift_width',default=SHIFT_WIDTH,type='int',help='Max shift width in target pixel (%default)')
parser.add_option('--shift_height',default=SHIFT_HEIGHT,type='int',help='Max shift height in target pixel (%default)')
parser.add_option('--margin_width',default=MARGIN_WIDTH,type='int',help='Margin width in target pixel (%default)')
parser.add_option('--margin_height',default=MARGIN_HEIGHT,type='int',help='Margin height in target pixel (%default)')
parser.add_option('-r','--rthr',default=RTHR,type='float',help='Threshold of correlation coefficient (%default)')
parser.add_option('-E','--feps',default=FEPS,type='float',help='Step length for curve_fit (%default)')
parser.add_option('-d','--debug',default=False,action='store_true',help='Debug mode (%default)')
(opts,args) = parser.parse_args()
if len(args) < 1:
    parser.print_help()
    sys.exit(0)
ref_fnam = args[0]
trg_fnam = args[1]

def residuals(p,refx,refy,refz,trgx,trgy,trgz,pmax):
    if opts.debug:
        sys.stderr.write('{}\n'.format(p))
    pabs = np.abs(p).max()
    if pabs > pmax:
        return np.full(2,2.0+pabs-pmax)
    f = interp2d(trgx+p[0],trgy+p[1],trgz,kind='linear')
    intz = f(refx,refy)[::-1]
    r = np.corrcoef(intz.flatten(),refz.flatten())[0,1]
    return np.full(2,1.0-r)

#tif_tags = {}
#with tifffile.TiffFile(ref_fnam) as tif:
#    for tag in tif.pages[0].tags.values():
#        name,value = tag.name,tag.value
#        tif_tags[name] = value
#ref_bands = []
#for line in tif_tags['ImageDescription'].splitlines():
#    m = re.search('(\d+)\s*;',line);
#    if m:
#        band = m.group(1)
#        ref_bands.append(band)
#ref_band = ref_bands.index('8')
#ref_band = 7

ds = gdal.Open(ref_fnam)
prj = ds.GetProjection()
srs = osr.SpatialReference(wkt=prj)
ref_data = ds.ReadAsArray()[opts.ref_band]
trans = ds.GetGeoTransform() # maybe obtained from tif_tags['ModelTransformationTag']
indy,indx = np.indices(ref_data.shape)
if srs.IsProjected():
    pnam = srs.GetAttrValue('projcs')
    if re.search('UTM',pnam):
        ref_xp = trans[0]+(indx+0.5)*trans[1]+(indy+0.5)*trans[2]
        ref_yp = trans[3]+(indx+0.5)*trans[4]+(indy+0.5)*trans[5]
    else:
        raise ValueError('Unsupported projection >>> '+pnam)
else:
    raise ValueError('Unsupported projection >>> '+prj)
ds = None # close dataset
ref_xp0 = ref_xp[0,:]
ref_yp0 = ref_yp[:,0]
ref_xp_min = ref_xp0.min()
ref_xp_max = ref_xp0.max()
ref_yp_min = ref_yp0.min()
ref_yp_max = ref_yp0.max()
ref_xp_stp = ref_xp0[1]-ref_xp0[0]
ref_yp_stp = ref_yp0[1]-ref_yp0[0]
if ref_xp_stp <= 0.0:
    raise ValueError('Error, ref_xp_stp={}'.format(ref_xp_stp))
if ref_yp_stp >= 0.0:
    raise ValueError('Error, ref_yp_stp={}'.format(ref_yp_stp))

#tif_tags = {}
#with tifffile.TiffFile(trg_fnam) as tif:
#    for tag in tif.pages[0].tags.values():
#        name,value = tag.name,tag.value
#        tif_tags[name] = value
#root = ET.fromstring(tif_tags['65000'])
#trg_bands = []
#for i,value in enumerate(root.iter('BAND_NAME')):
#    band = value.text
#    trg_bands.append(band)
#trg_band = trg_bands.index('B8')

ds = gdal.Open(trg_fnam)
prj = ds.GetProjection()
srs = osr.SpatialReference(wkt=prj)
trg_data = ds.ReadAsArray()[opts.trg_band]
trans = ds.GetGeoTransform() # maybe obtained from tif_tags['ModelTransformationTag']
indy,indx = np.indices(trg_data.shape)
if srs.IsProjected():
    pnam = srs.GetAttrValue('projcs')
    if re.search('UTM',pnam):
        trg_xp = trans[0]+(indx+0.5)*trans[1]+(indy+0.5)*trans[2]
        trg_yp = trans[3]+(indx+0.5)*trans[4]+(indy+0.5)*trans[5]
    else:
        raise ValueError('Unsupported projection >>> '+pnam)
else:
    raise ValueError('Unsupported projection >>> '+prj)
ds = None # close dataset
trg_xp0 = trg_xp[0,:]
trg_yp0 = trg_yp[:,0]
trg_xp_min = trg_xp0.min()
trg_xp_max = trg_xp0.max()
trg_yp_min = trg_yp0.min()
trg_yp_max = trg_yp0.max()
trg_xp_stp = trg_xp0[1]-trg_xp0[0]
trg_yp_stp = trg_yp0[1]-trg_yp0[0]
if trg_xp_stp <= 0.0:
    raise ValueError('Error, trg_xp_stp={}'.format(trg_xp_stp))
if trg_yp_stp >= 0.0:
    raise ValueError('Error, trg_yp_stp={}'.format(trg_yp_stp))

ref_height,ref_width = ref_data.shape
trg_height,trg_width = trg_data.shape

subset_half_width = opts.subset_width//2
subset_half_height = opts.subset_height//2
for trg_indyc in np.arange(0,trg_height,subset_half_height):
#for trg_indyc in [450]:
    trg_indy1 = trg_indyc-subset_half_height-opts.margin_height
    trg_indy2 = trg_indyc+subset_half_height+opts.margin_height+1
    if trg_indy1 < 0:
        continue
    if trg_indy2 > trg_height:
        break
    ref_yp1 = trg_yp0[trg_indyc-subset_half_height] # yp1 > ypc
    ref_yp2 = trg_yp0[trg_indyc+subset_half_height] # yp2 < ypc
    if ref_yp1 > ref_yp_max:
        continue
    if ref_yp2 < ref_yp_min:
        continue
    ref_indy1 = np.where(ref_yp0 <= ref_yp1)[0]
    if ref_indy1.size < 1:
        continue
    ref_indy1 = ref_indy1[0]
    ref_indy2 = np.where(ref_yp0 >= ref_yp2)[0]
    if ref_indy2.size < 1:
        continue
    ref_indy2 = ref_indy2[-1]+1
    for trg_indxc in np.arange(0,trg_width,subset_half_width):
    #for trg_indxc in [1250]:
        trg_indx1 = trg_indxc-subset_half_width-opts.margin_width
        trg_indx2 = trg_indxc+subset_half_width+opts.margin_width+1
        if trg_indx1 < 0:
            continue
        if trg_indx2 > trg_width:
            break
        ref_xp1 = trg_xp0[trg_indxc-subset_half_width]
        ref_xp2 = trg_xp0[trg_indxc+subset_half_width]
        if ref_xp1 < ref_xp_min:
            continue
        if ref_xp2 > ref_xp_max:
            continue
        ref_indx1 = np.where(ref_xp0 >= ref_xp1)[0]
        if ref_indx1.size < 1:
            continue
        ref_indx1 = ref_indx1[0]
        ref_indx2 = np.where(ref_xp0 <= ref_xp2)[0]
        if ref_indx2.size < 1:
            continue
        ref_indx2 = ref_indx2[-1]+1
        # target subset
        trg_subset_xp0 = trg_xp0[trg_indx1:trg_indx2]
        trg_subset_yp0 = trg_yp0[trg_indy1:trg_indy2]
        trg_subset_data = trg_data[trg_indy1:trg_indy2,trg_indx1:trg_indx2]
        # reference subset
        ref_subset_xp0 = ref_xp0[ref_indx1:ref_indx2]
        ref_subset_yp0 = ref_yp0[ref_indy1:ref_indy2]
        ref_subset_data = ref_data[ref_indy1:ref_indy2,ref_indx1:ref_indx2]
        if ref_subset_data.min() < 1:
            continue
        p1 = np.array([0.0,0.0])
        rmax = -1.0e10
        for i in range(-opts.shift_height,opts.shift_height+1):
            for j in range(-opts.shift_width,opts.shift_width+1):
                p2 = np.array([np.abs(trg_yp_stp)*i,np.abs(trg_xp_stp)*j])
                r = 1.0-residuals(p2,ref_subset_xp0,ref_subset_yp0,ref_subset_data,
                                  trg_subset_xp0,trg_subset_yp0,trg_subset_data,1.0e10)[0]
                if r > rmax:
                    rmax = r
                    p1 = p2.copy()
        result = leastsq(residuals,p1,args=(ref_subset_xp0,ref_subset_yp0,ref_subset_data,
                                            trg_subset_xp0,trg_subset_yp0,trg_subset_data,
                                            min(np.abs(trg_xp_stp*opts.shift_width),np.abs(trg_yp_stp*opts.shift_height))),
                                            epsfcn=opts.feps,full_output=True)
        p2 = result[0]
        r = 1.0-result[2]['fvec'][0]
        #f = interp2d(trg_subset_xp0+p2[0],trg_subset_yp0+p2[1],trg_subset_data,kind='linear')
        #trg_interp_data = f(ref_subset_xp0,ref_subset_yp0)[::-1]
        #r = np.corrcoef(trg_interp_data.flatten(),ref_subset_data.flatten())[0,1]
        if r > opts.rthr:
            sys.stdout.write('{:6d} {:6d} {:8.2f} {:8.2f} {:6.2f} {:6.2f} {:8.3f}\n'.format(trg_indxc,trg_indyc,trg_xp0[trg_indxc]+p2[0],trg_yp0[trg_indyc]+p2[1],p2[0],p2[1],r))
        #break
    #break
