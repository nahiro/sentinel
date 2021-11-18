#!/usr/bin/env python
import sys
import gdal
import numpy as np
import shapefile
from shapely.geometry import Point,Polygon
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
from matplotlib.backends.backend_pdf import PdfPages
from optparse import OptionParser,IndentedHelpFormatter

# Default values
DATNAM = 'pixel_distance.dat'
FIGNAM = 'pixel_distance.pdf'

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('-f','--img_fnam',default=None,help='Image file name (%default)')
parser.add_option('-s','--shp_fnam',default=None,help='Shape file name (%default)')
parser.add_option('-b','--blk_fnam',default=None,help='Block file name (%default)')
parser.add_option('-B','--block',default=None,help='Block name (%default)')
parser.add_option('--max_distance',default=None,type='float',help='Buffer distance (%default)')
parser.add_option('--radius',default=None,type='float',help='Radius to be considered (%default)')
parser.add_option('--xmgn',default=None,type='float',help='X margin (%default)')
parser.add_option('--ymgn',default=None,type='float',help='Y margin (%default)')
parser.add_option('--use_index',default=False,action='store_true',help='Use index instead of OBJECTID (%default)')
parser.add_option('--use_objectid',default=False,action='store_true',help='Use OBJECTID instead of Block (%default)')
parser.add_option('-d','--debug',default=False,action='store_true',help='Debug mode (%default)')
parser.add_option('-c','--check',default=False,action='store_true',help='Check mode (%default)')
parser.add_option('-o','--datnam',default=DATNAM,help='Output data name (%default)')
parser.add_option('-F','--fignam',default=FIGNAM,help='Output figure name for debug (%default)')
(opts,args) = parser.parse_args()

ds = gdal.Open(opts.img_fnam)
data = ds.ReadAsArray()
if ds.RasterCount < 2:
    data_shape = data.shape
else:
    data_shape = data[0].shape
trans = ds.GetGeoTransform() # maybe obtained from tif_tags['ModelTransformationTag']
indy,indx = np.indices(data_shape)
xp = trans[0]+(indx+0.5)*trans[1]+(indy+0.5)*trans[2]
yp = trans[3]+(indx+0.5)*trans[4]+(indy+0.5)*trans[5]
ds = None
xstp = abs(xp[0,1]-xp[0,0])
ystp = abs(yp[1,0]-yp[0,0])
xhlf = 0.5*xstp
yhlf = 0.5*ystp
if opts.max_distance is None:
    opts.max_distance = max(xstp,ystp)
if opts.radius is None:
    opts.radius = opts.max_distance*1.5
if opts.xmgn is None:
    opts.xmgn = xstp*6.0
if opts.ymgn is None:
    opts.ymgn = ystp*6.0

r = shapefile.Reader(opts.shp_fnam)

if opts.blk_fnam is not None:
    block = {}
    with open(opts.blk_fnam,'r') as fp:
        for line in fp:
            item = line.split()
            block.update({int(item[0]):item[1]})
    if len(block) != len(r):
        raise ValueError('Error, len(block)={}, len(r)={}'.format(len(block),len(r)))

if opts.debug or opts.check:
    plt.interactive(True)
    fig = plt.figure(1,facecolor='w',figsize=(6,3.5))
    plt.subplots_adjust(top=0.85,bottom=0.20,left=0.15,right=0.95)
    pdf = PdfPages(opts.fignam)
with open(opts.datnam,'w') as fp:
    for ii,shaperec in enumerate(r.iterShapeRecords()):
        if ii%100 == 0:
            sys.stderr.write('{}\n'.format(ii))
            sys.stderr.flush()
        rec = shaperec.record
        shp = shaperec.shape
        if opts.use_index:
            object_id = ii+1
        else:
            object_id = rec.OBJECTID
        if len(shp.points) < 1:
            sys.stderr.write('Warning, len(shp.points)={}, ii={}\n'.format(len(shp.points),ii))
            continue
        poly_original = Polygon(shp.points)
        path_original = Path(shp.points)
        poly_outer = Polygon(shp.points).buffer(opts.radius)
        path_outer = []
        flags_outer = []
        if poly_outer.area <= 0.0:
            sys.stderr.write('Warning, poly_outer.area={} >>> FID {}, OBJECTID {}\n'.format(poly_outer.area,ii,object_id))
        elif poly_outer.type == 'MultiPolygon':
            sys.stderr.write('Warning, poly_outer.type={} >>> FID {}, OBJECTID {}\n'.format(poly_outer.type,ii,object_id))
            for p in poly_outer:
                p_outer = Path(np.array(p.exterior.coords.xy).swapaxes(0,1))
                path_outer.append(p_outer)
                if len(flags_outer) < 1:
                    flags_outer = p_outer.contains_points(np.hstack((xp.reshape(-1,1),yp.reshape(-1,1))),radius=0.0).reshape(data_shape)
                else:
                    flags_outer |= p_outer.contains_points(np.hstack((xp.reshape(-1,1),yp.reshape(-1,1))),radius=0.0).reshape(data_shape)
        else:
            path_outer = Path(np.array(poly_outer.exterior.coords.xy).swapaxes(0,1))
            flags_outer = path_outer.contains_points(np.hstack((xp.reshape(-1,1),yp.reshape(-1,1))),radius=0.0).reshape(data_shape)
        poly_inner = Polygon(shp.points).buffer(-opts.radius)
        path_inner = []
        flags_inner = []
        if poly_inner.area <= 0.0:
            sys.stderr.write('Warning, poly_inner.area={} >>> FID {}, OBJECTID {}\n'.format(poly_inner.area,ii,object_id))
        elif poly_inner.type == 'MultiPolygon':
            sys.stderr.write('Warning, poly_inner.type={} >>> FID {}, OBJECTID {}\n'.format(poly_inner.type,ii,object_id))
            for p in poly_inner:
                p_inner = Path(np.array(p.exterior.coords.xy).swapaxes(0,1))
                path_inner.append(p_inner)
                if len(flags_inner) < 1:
                    flags_inner = p_inner.contains_points(np.hstack((xp.reshape(-1,1),yp.reshape(-1,1))),radius=0.0).reshape(data_shape)
                else:
                    flags_inner |= p_inner.contains_points(np.hstack((xp.reshape(-1,1),yp.reshape(-1,1))),radius=0.0).reshape(data_shape)
        else:
            path_inner = Path(np.array(poly_inner.exterior.coords.xy).swapaxes(0,1))
            flags_inner = path_inner.contains_points(np.hstack((xp.reshape(-1,1),yp.reshape(-1,1))),radius=0.0).reshape(data_shape)
        if len(flags_inner) < 1:
            flags = flags_outer
        else:
            flags = flags_outer & (~flags_inner)
        if opts.debug or opts.check:
            flags_inside = []
            flags_near = []
            path_pixels = []
        inds = []
        dists = []
        err = False
        for ix,iy in zip(indx[flags],indy[flags]):
            xc = xp[iy,ix]
            yc = yp[iy,ix]
            pc = Point(xc,yc)
            try:
                dist = poly_original.exterior.distance(pc)
            except Exception:
                sys.stderr.write('Warning, error occured at (ix,iy)=({},{}), ii={}\n'.format(ix,iy,ii))
                err = True
                continue
            if dist < opts.max_distance:
                inds.append(np.ravel_multi_index((iy,ix),data_shape))
                dists.append(dist)
            if opts.debug or opts.check:
                flags_inside.append(pc.within(poly_outer) and not pc.within(poly_inner)) # should be True
                flags_near.append(dist < opts.max_distance)
                path_pixel = Path([(xc-xhlf,yc-yhlf),(xc-xhlf,yc+yhlf),(xc+xhlf,yc+yhlf),(xc+xhlf,yc-yhlf),(xc-xhlf,yc-yhlf)])
                path_pixels.append(path_pixel)
        if err:
            continue
        inds = np.array(inds)
        dists = np.array(dists)
        x1,y1,x2,y2 = shp.bbox
        xctr = 0.5*(x1+x2)
        yctr = 0.5*(y1+y2)
        flag_check = Point(xctr,yctr).within(poly_original)
        if not flag_check:
            sys.stderr.write('Warning, center pixel is not included >>> FID: {}, OBJECTID: {}\n'.format(ii,object_id))
        # output results ###
        if opts.use_objectid:
            fp.write('{} {} {}'.format(object_id,object_id,len(inds)))
        elif opts.blk_fnam is not None:
            fp.write('{} {} {}'.format(object_id,block[object_id],len(inds)))
        elif opts.block is not None:
            fp.write('{} {} {}'.format(object_id,opts.block,len(inds)))
        else:
            fp.write('{} {}'.format(object_id,len(inds)))
        isort = np.argsort(dists)
        for ind,dist in zip(inds[isort],dists[isort]):
            fp.write(' {:d} {:.6e}'.format(ind,dist))
        fp.write('\n')
        ####################
        if opts.debug or (opts.check and not flag_check):
            fig.clear()
            ax1 = plt.subplot(111)
            ax1.set_title('OBJECTID: {}'.format(object_id))
            for path_pixel in path_pixels:
                patch = patches.PathPatch(path_pixel,facecolor='none',lw=1)
                ax1.add_patch(patch)
            patch = patches.PathPatch(path_original,facecolor='none',lw=2)
            ax1.add_patch(patch)
            if poly_outer.area <= 0.0:
                pass
            elif poly_outer.type == 'MultiPolygon':
                for p_outer in path_outer:
                    patch = patches.PathPatch(p_outer,facecolor='none',lw=2,ls='--')
                    ax1.add_patch(patch)
            else:
                patch = patches.PathPatch(path_outer,facecolor='none',lw=2,ls='--')
                ax1.add_patch(patch)
            if poly_inner.area <= 0.0:
                pass
            elif poly_inner.type == 'MultiPolygon':
                for p_inner in path_inner:
                    patch = patches.PathPatch(p_inner,facecolor='none',edgecolor='#888888',lw=2)
                    ax1.add_patch(patch)
            else:
                patch = patches.PathPatch(path_inner,facecolor='none',edgecolor='#888888',lw=2)
                ax1.add_patch(patch)
            ax1.plot(xp,yp,'o',color='#888888')
            xmin = 1.0e10
            xmax = -1.0e10
            ymin = 1.0e10
            ymax = -1.0e10
            for j,(x,y) in enumerate(zip(xp[flags],yp[flags])):
                if flags_near[j]:
                    ax1.plot(x,y,'mo')
                elif flags_inside[j]:
                    ax1.plot(x,y,'ro')
                ax1.text(x,y,str(j),ha='center',va='center')
                xmin = min(xmin,x)
                xmax = max(xmax,x)
                ymin = min(ymin,y)
                ymax = max(ymax,y)
            ax1.plot([xctr],[yctr],'rx',ms=10)
            for pp in shp.points:
                ax1.plot(pp[0],pp[1],'bo')
            ax1.set_xlim(xmin-opts.xmgn,xmax+opts.xmgn)
            ax1.set_ylim(ymin-opts.ymgn,ymax+opts.ymgn)
            ax1.ticklabel_format(useOffset=False,style='plain')
            plt.savefig(pdf,format='pdf')
            plt.draw()
            plt.pause(0.1)
        #break
if opts.debug or opts.check:
    pdf.close()
