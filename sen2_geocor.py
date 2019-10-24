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

def residuals(p,refx,refy,refz,trgx,trgy,trgz):
    pmax = np.abs(p).max()
    if pmax > 50.0:
        return np.full(10,1.0+pmax-50.0)
    f = interp2d(trgx+p[0],trgy+p[1],trgz,kind='linear')
    intz = f(refx,refy)[::-1]
    r = np.corrcoef(intz.flatten(),refz.flatten())[0,1]
    return np.full(10,1.0-r)

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
ref_band = 7

ds = gdal.Open(ref_fnam)
prj = ds.GetProjection()
srs = osr.SpatialReference(wkt=prj)
ref_data = ds.ReadAsArray()[ref_band]
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
ref_yp_stp = ref_yp0[1]-ref_yp0[0]
if ref_yp_stp >= 0.0:
    raise ValueError('Error, ref_yp_stp={}'.format(ref_yp_stp))

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
trg_data = ds.ReadAsArray()[trg_band]
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
trg_yp_stp = trg_yp0[1]-trg_yp0[0]
if trg_yp_stp >= 0.0:
    raise ValueError('Error, trg_yp_stp={}'.format(trg_yp_stp))

ref_height,ref_width = ref_data.shape
trg_height,trg_width = trg_data.shape

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
    ref_yp1 = trg_yp0[trg_indyc-template_half_width] # yp1 > ypc
    ref_yp2 = trg_yp0[trg_indyc+template_half_width] # yp2 < ypc
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
    for trg_indxc in np.arange(0,trg_width,template_half_width):
        trg_indx1 = trg_indxc-template_half_width-margin_width
        trg_indx2 = trg_indxc+template_half_width+margin_width+1
        if trg_indx1 < 0:
            continue
        if trg_indx2 > trg_width:
            break
        ref_xp1 = trg_xp0[trg_indxc-template_half_width]
        ref_xp2 = trg_xp0[trg_indxc+template_half_width]
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
        result = leastsq(residuals,p1,args=(ref_subset_xp0,ref_subset_yp0,ref_subset_data,trg_subset_xp0,trg_subset_yp0,trg_subset_data),epsfcn=1.0e-1,full_output=True)
        p2 = result[0]
        r = 1.0-result[2]['fvec'][0]
        sys.stdout.write('{:6d} {:6d} {:8.2f} {:8.2f} {:6.2f} {:6.2f} {:8.3f}\n'.format(trg_indxc,trg_indyc,trg_xp0[trg_indxc]+p2[0],trg_yp0[trg_indyc]+p2[1],p2[0],p2[1],r))
        #break
    #break
