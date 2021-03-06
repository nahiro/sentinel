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
DATNAM = 'pixel_area.dat'
FIGNAM = 'pixel_area.pdf'

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('-f','--img_fnam',default=None,help='Image file name (%default)')
parser.add_option('-s','--shp_fnam',default=None,help='Shape file name (%default)')
parser.add_option('-b','--blk_fnam',default=None,help='Block file name (%default)')
parser.add_option('-B','--block',default=None,help='Block name (%default)')
parser.add_option('--buffer',default=None,type='float',help='Buffer distance (%default)')
parser.add_option('--use_index',default=False,action='store_true',help='Use index instead of OBJECTID (%default)')
parser.add_option('--use_objectid',default=False,action='store_true',help='Use OBJECTID instead of Block (%default)')
parser.add_option('-d','--debug',default=False,action='store_true',help='Debug mode (%default)')
parser.add_option('-c','--check',default=False,action='store_true',help='Check mode (%default)')
parser.add_option('-o','--datnam',default=DATNAM,help='Output data name (%default)')
parser.add_option('-F','--fignam',default=FIGNAM,help='Output figure name for debug (%default)')
(opts,args) = parser.parse_args()

ds = gdal.Open(opts.img_fnam)
data = ds.ReadAsArray()
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
rstp = max(xstp,ystp)

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
        p = Path(shp.points)
        if opts.buffer is not None:
            p1 = Polygon(shp.points).buffer(opts.buffer)
        else:
            p1 = Polygon(shp.points)
        flags = p.contains_points(np.hstack((xp.reshape(-1,1),yp.reshape(-1,1))),radius=-3.0*rstp).reshape(data_shape)
        if opts.debug or opts.check:
            flags_inside = []
            flags_near = []
            p4_paths = []
        inds = []
        rats = []
        err = False
        for ix,iy in zip(indx[flags],indy[flags]):
            xc = xp[iy,ix]
            yc = yp[iy,ix]
            pc = Point(xc,yc)
            p2 = Polygon([(xc-xhlf,yc-yhlf),(xc-xhlf,yc+yhlf),(xc+xhlf,yc+yhlf),(xc+xhlf,yc-yhlf),(xc-xhlf,yc-yhlf)])
            try:
                p3 = p1.intersection(p2)
            except Exception:
                sys.stderr.write('Warning, error occured at (ix,iy)=({},{}), ii={}\n'.format(ix,iy,ii))
                err = True
                continue
            rat = p3.area/p2.area
            if rat > 1.0e-10:
                inds.append(np.ravel_multi_index((iy,ix),data_shape))
                rats.append(rat)
            if opts.debug or opts.check:
                flags_inside.append(pc.within(p1))
                flags_near.append(rat > 1.0e-10)
                p4 = Path([(xc-xhlf,yc-yhlf),(xc-xhlf,yc+yhlf),(xc+xhlf,yc+yhlf),(xc+xhlf,yc-yhlf),(xc-xhlf,yc-yhlf)])
                p4_paths.append(p4)
        if err:
            continue
        inds = np.array(inds)
        rats = np.array(rats)
        x1,y1,x2,y2 = shp.bbox
        xctr = 0.5*(x1+x2)
        yctr = 0.5*(y1+y2)
        ictr = np.argmin(np.square(xp-xctr)+np.square(yp-yctr))
        if not ictr in inds:
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
        isort = np.argsort(rats)[::-1]
        for ind,rat in zip(inds[isort],rats[isort]):
            fp.write(' {:d} {:.6e}'.format(ind,rat))
        fp.write('\n')
        ####################
        if opts.debug or (opts.check and not ictr in inds):
            fig.clear()
            ax1 = plt.subplot(111)
            ax1.set_title('OBJECTID: {}'.format(object_id))
            for p4 in p4_paths:
                patch = patches.PathPatch(p4,facecolor='none',lw=1)
                ax1.add_patch(patch)

            patch = patches.PathPatch(p,facecolor='none',lw=2)
            ax1.add_patch(patch)
            ax1.plot(xp,yp,'o',color='#888888')
            for j,(x,y) in enumerate(zip(xp[flags],yp[flags])):
                if flags_inside[j]:
                    ax1.plot(x,y,'ro')
                elif flags_near[j]:
                    ax1.plot(x,y,'mo')
                ax1.text(x,y,str(j),ha='center',va='center')
            ax1.plot([xctr],[yctr],'rx',ms=10)
            xmin = 1.0e10
            xmax = -1.0e10
            ymin = 1.0e10
            ymax = -1.0e10
            for pp in shp.points:
                ax1.plot(pp[0],pp[1],'bo')
                xmin = min(xmin,pp[0])
                xmax = max(xmax,pp[0])
                ymin = min(ymin,pp[1])
                ymax = max(ymax,pp[1])
            ax1.set_xlim(xmin-6.0*xstp,xmax+6.0*xstp)
            ax1.set_ylim(ymin-6.0*ystp,ymax+6.0*ystp)
            ax1.ticklabel_format(useOffset=False,style='plain')
            plt.savefig(pdf,format='pdf')
            plt.draw()
            plt.pause(0.1)
if opts.debug or opts.check:
    pdf.close()
