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

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('-f','--img_fnam',default=None,help='Image file name (%default)')
parser.add_option('-s','--shp_fnam',default=None,help='Shape file name (%default)')
parser.add_option('-d','--debug',default=False,action='store_true',help='Debug mode (%default)')
(opts,args) = parser.parse_args()

ds = gdal.Open(opts.img_fnam)
data = ds.ReadAsArray()
data_shape = data[0].shape
trans = ds.GetGeoTransform() # maybe obtained from tif_tags['ModelTransformationTag']
indy,indx = np.indices(data_shape)
xp = trans[0]+(indx+0.5)*trans[1]+(indy+0.5)*trans[2]
yp = trans[3]+(indx+0.5)*trans[4]+(indy+0.5)*trans[5]
ds = None

r = shapefile.Reader(opts.shp_fnam)

if opts.debug:
    plt.interactive(True)
    fig = plt.figure(1,facecolor='w',figsize=(6,3.5))
    plt.subplots_adjust(top=0.85,bottom=0.20,left=0.15,right=0.95)
    pdf = PdfPages('pixel_area.pdf')

for ii,shaperec in enumerate(r.iterShapeRecords()):
    if ii%100 == 0:
        sys.stderr.write('{}\n'.format(ii))
        sys.stderr.flush()
    rec = shaperec.record
    shp = shaperec.shape
    p = Path(shp.points)
    p1 = Polygon(shp.points)
    flags = p.contains_points(np.hstack((xp.reshape(-1,1),yp.reshape(-1,1))),radius=-30.0).reshape(data_shape)
    if opts.debug:
        fig.clear()
        ax1 = plt.subplot(111)
        ax1.set_title('OBJECTID: {}'.format(rec[0]))
        flags_inside = []
        flags_near = []
    inds = []
    rats = []
    for ix,iy in zip(indx[flags],indy[flags]):
        xc = xp[iy,ix]
        yc = yp[iy,ix]
        pc = Point(xc,yc)
        p2 = Polygon([(xc-5.0,yc-5.0),(xc-5.0,yc+5.0),(xc+5.0,yc+5.0),(xc+5.0,yc-5.0),(xc-5.0,yc-5.0)])
        p3 = p1.intersection(p2)
        rat = p3.area/p2.area
        if rat > 1.0e-10:
            inds.append(np.ravel_multi_index((iy,ix),data_shape))
            rats.append(rat)
        if opts.debug:
            flags_inside.append(pc.within(p1))
            flags_near.append(rat > 1.0e-10)
            p4 = Path([(xc-5.0,yc-5.0),(xc-5.0,yc+5.0),(xc+5.0,yc+5.0),(xc+5.0,yc-5.0),(xc-5.0,yc-5.0)])
            patch = patches.PathPatch(p4,facecolor='none',lw=1)
            ax1.add_patch(patch)
    inds = np.array(inds)
    rats = np.array(rats)
    # output results ###
    sys.stdout.write('{} {}'.format(rec[0],len(inds)))
    isort = np.argsort(rats)[::-1]
    for ind,rat in zip(inds[isort],rats[isort]):
        sys.stdout.write(' {:d} {:.6e}'.format(ind,rat))
    sys.stdout.write('\n')
    ####################
    if opts.debug:
        patch = patches.PathPatch(p,facecolor='none',lw=2)
        ax1.add_patch(patch)
        ax1.plot(xp,yp,'o',color='#888888')
        for j,(x,y) in enumerate(zip(xp[flags],yp[flags])):
            if flags_inside[j]:
                ax1.plot(x,y,'ro')
            elif flags_near[j]:
                ax1.plot(x,y,'mo')
            ax1.text(x,y,str(j),ha='center',va='center')
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
        ax1.set_xlim(xmin-60.0,xmax+60.0)
        ax1.set_ylim(ymin-60.0,ymax+60.0)
        ax1.ticklabel_format(useOffset=False,style='plain')
        plt.draw()
        plt.savefig(pdf,format='pdf')
#       plt.pause(0.1)
if opts.debug:
    pdf.close()
