#!/usr/bin/env python
import os
import sys
import re
from datetime import datetime,timedelta
from dateutil.relativedelta import relativedelta
from subprocess import call
from optparse import OptionParser,IndentedHelpFormatter

# Default values
HOME = os.environ.get('HOME')
if HOME is None:
    HOME = os.environ.get('HOMEPATH')
SCRDIR = os.path.join(HOME,'Script')
DATDIR = os.path.join(HOME,'Work','Sentinel-1')
WRKDIR = os.path.join(HOME,'Work','SATREPS','Transplanting_date')
END = datetime.now().strftime('%Y%m%d')
SITES = ['Cihea','Bojongsoang']

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.set_usage('Usage: %prog collocated_geotiff_file [options]')
parser.add_option('--scrdir',default=SCRDIR,help='Script directory (%default)')
parser.add_option('--datdir',default=DATDIR,help='Data directory (%default)')
parser.add_option('--wrkdir',default=WRKDIR,help='Work directory (%default)')
parser.add_option('-s','--str',default=None,help='Start date of download in the format YYYYMMDD (%default)')
parser.add_option('-e','--end',default=END,help='End date of download in the format YYYYMMDD (%default)')
parser.add_option('-S','--sites',default=None,action='append',help='Target sites ({})'.format(SITES))
parser.add_option('--test',default=False,action='store_true',help='Test mode (%default)')
parser.add_option('-d','--debug',default=False,action='store_true',help='Debug mode (%default)')
(opts,args) = parser.parse_args()
if opts.sites is None:
    opts.sites = SITES

dtims = {}
for site in opts.sites:
    dtim_list = []
    for f in sorted(os.listdir(os.path.join(opts.datdir,site,'sigma0_speckle'))):
        m = re.search('('+'\d'*8+')_resample.tif',f)
        if not m:
            continue
        dtim_list.append(datetime.strptime(m.group(1),'%Y%m%d'))
    dtims.update({site:dtim_list})

topdir = os.getcwd()
for site in opts.sites:
    d3 = datetime.strptime(opts.end,'%Y%m%d')
    for d2 in dtims[site][::-1]:
        if d2 < d3:
            break
    if opts.str is None:
        d1 = d2
    else:
        d3 = datetime.strptime(opts.str,'%Y%m%d')
        for d1 in dtims[site][::-1]:
            if d1 > d3:
                break
    for d in dtims[site]:
        if d < d1:
            continue
        if d > d2:
            break
        dstr = d.strftime('%Y%m%d')
        tmin = (d-timedelta(days=90)).strftime('%Y%m%d')
        tmax = dstr
        data_tmin = (d-timedelta(days=150)).strftime('%Y%m%d')
        data_tmax = dstr
        #sys.stderr.write(dstr+'\n')
        wrkdir = os.path.join(opts.wrkdir,site,'preliminary',dstr)
        if not os.path.exists(wrkdir):
            os.makedirs(wrkdir)
        if not os.path.isdir(wrkdir):
            raise IOError('Error, no such directory >>> '+wrkdir)
        os.chdir(wrkdir)
        tif_fnam = os.path.join(wrkdir,'trans_date_{}_{}_preliminary.tif'.format(site,dstr))
        shp_bnam = os.path.join(wrkdir,'trans_date_{}_{}_preliminary'.format(site,dstr))
        shp_fnam = shp_bnam+'.shp'
        trans_pixel_image = os.path.join(wrkdir,'trans_pixel_{}_{}_preliminary.pdf'.format(site,dstr))
        trans_field_image = os.path.join(wrkdir,'trans_field_{}_{}_preliminary.pdf'.format(site,dstr))
        file_list = []
        try:
            command = 'python'
            command += ' '+os.path.join(opts.scrdir,'calc_trans_date.py')
            #command += ' -x 600 -X 900 -y 680 -Y 1030'
            command += ' --tmin '+tmin
            command += ' --tmax '+tmax
            command += ' --data_tmin '+data_tmin
            command += ' --data_tmax '+data_tmax
            command += ' -D '+os.path.join(opts.datdir,site,'sigma0_speckle')
            command += ' --search_key resample'
            command += ' --near_fnam '+os.path.join(opts.wrkdir,site,'find_nearest.npz')
            command += ' --incidence_list '+os.path.join(opts.wrkdir,site,'incidence_list.dat')
            command += ' --out_fnam '+tif_fnam
            command += ' --early'
            command += ' 2>'+os.path.join(wrkdir,'err')
            sys.stderr.write(command+'\n')
            call(command,shell=True)
            if os.path.exists(tif_fnam):
                file_list.append(tif_fnam)
                command = 'python'
                command += ' '+os.path.join(opts.scrdir,'field_mean_shapefile.py')
                command += ' --data_file '+tif_fnam
                command += ' --area_file '+os.path.join(opts.wrkdir,site,'pixel_area_block.dat')
                command += ' --outnam '+shp_bnam
                call(command,shell=True)
                command = 'python'
                command += ' '+os.path.join(opts.scrdir,'draw_trans_pixel_{}.py'.format(site.lower()))
                command += ' --tmin '+tmin
                command += ' --tmax '+tmax
                command += ' --pmin 0'
                command += ' --pmax 300'
                command += ' --title "Search Period: {} - {}"'.format(tmin,tmax)
                command += ' --trans_fnam '+tif_fnam
                command += ' --mask_fnam '+os.path.join(opts.wrkdir,site,'paddy_mask.tif')
                command += ' --output_fnam '+trans_pixel_image
                command += ' --batch'
                call(command,shell=True)
                if os.path.exists(trans_pixel_image):
                    file_list.append(trans_pixel_image)
            if os.path.exists(shp_fnam):
                file_list.append(shp_fnam)
                file_list.append(shp_bnam+'.dbf')
                file_list.append(shp_bnam+'.prj')
                file_list.append(shp_bnam+'.shx')
                command = 'python'
                command += ' '+os.path.join(opts.scrdir,'draw_trans_field_{}.py'.format(site.lower()))
                command += ' --tmin '+tmin
                command += ' --tmax '+tmax
                command += ' --pmin 0'
                command += ' --pmax 300'
                command += ' --title "Search Period: {} - {}"'.format(tmin,tmax)
                command += ' --trans_fnam '+shp_fnam
                command += ' --output_fnam '+trans_field_image
                command += ' --batch'
                call(command,shell=True)
                if os.path.exists(trans_field_image):
                    file_list.append(trans_field_image)
            if len(file_list) > 0:
                command = 'python'
                command += ' '+os.path.join(opts.scrdir,'trans_date_upload.py')
                command += ' --site '+site
                command += ' --date '+dstr
                if opts.test:
                    command += ' --level test'
                else:
                    command += ' --level preliminary'
                command += ' --overwrite'
                command += ' '+' '.join(file_list)
                call(command,shell=True)
        except Exception:
            sys.stderr.write('Warning, error occured.\n')
        d3 = d+relativedelta(months=2)
        d = datetime(d3.year,d3.month,1)-timedelta(days=1)
        os.chdir(topdir)
