import os
import sys
import re
import gdal
import osr
import tifffile
import xml.etree.ElementTree as ET
import numpy as np
from scipy.interpolate import interp2d
from scipy.optimize import leastsq

datdir = '/home/naohiro/Work/Sentinel-2/191022'

ref_fnam = os.path.join(datdir,'wv2_180629_mul.tif')

tif_tags = {}
with tifffile.TiffFile(ref_fnam) as tif:
    for tag in tif.pages[0].tags.values():
        name,value = tag.name,tag.value
        tif_tags[name] = value
ref_bands = []
for line in tif_tags['ImageDescription'].splitlines():
    m = re.search('(\d+)\s*;',line);
    if m:
        band = m.group(1)
        ref_bands.append(band)
ref_band = ref_bands.index('8')

ds = gdal.Open(ref_fnam)
prj = ds.GetProjection()
srs = osr.SpatialReference(wkt=prj)
ref_data = ds.ReadAsArray()
trans = ds.GetGeoTransform() # maybe obtained from tif_tags['ModelTransformationTag']
indy,indx = np.indices(ref_data[ref_band].shape)
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

trg_fnam = os.path.join(datdir,'20190525.tif')

tif_tags = {}
with tifffile.TiffFile(trg_fnam) as tif:
    for tag in tif.pages[0].tags.values():
        name,value = tag.name,tag.value
        tif_tags[name] = value
root = ET.fromstring(tif_tags['65000'])
trg_bands = []
for i,value in enumerate(root.iter('BAND_NAME')):
    band = value.text
    trg_bands.append(band)
trg_band = trg_bands.index('B8')

ds = gdal.Open(trg_fnam)
prj = ds.GetProjection()
srs = osr.SpatialReference(wkt=prj)
trg_data = ds.ReadAsArray()
trans = ds.GetGeoTransform() # maybe obtained from tif_tags['ModelTransformationTag']
indy,indx = np.indices(trg_data[trg_band].shape)
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

ref_height,ref_width = ref_data[ref_band].shape
trg_height,trg_width = trg_data[trg_band].shape

template_width = 100 # target pixel
template_height = 100 # target pixel
template_half_width = template_width//2
template_half_height = template_height//2
shift_width = 3 # target pixel
shift_height = 3 # target pixel
margin_width = 5 # target_pixel
margin_height = 5 # target_pixel
for trg_indyc in np.arange(0,trg_height,template_half_height):
    trg_indy1 = trg_indyc-template_half_height-margin_height
    trg_indy2 = trg_indyc+template_half_height+margin_height+1
    if trg_indy1 < 0:
        continue
    if trg_indy2 > trg_height:
        break
    for trg_indxc in np.arange(0,trg_width,template_half_width):
        trg_indx1 = trg_indxc-template_half_width-margin_width
        trg_indx2 = trg_indxc+template_half_width+margin_width+1
        if trg_indx1 < 0:
            continue
        if trg_indx2 > trg_width:
            break
        #print(trg_indyc,trg_indxc)
        """
        trg_subset_data = trg_data[] # subset of target
        trg_subset_xp = trg_xp[]
        trg_subset_yp = trg_yp[]
        ref_subset_data = ref_data[]
        ref_subset_xp = ref_xp[]
        ref_subset_yp = ref_yp[]
        interp2d(trg_subset_xp,trg_subset_yp,trg_subset_data,ref_subset_xp,ref_subset_yp)
        """
