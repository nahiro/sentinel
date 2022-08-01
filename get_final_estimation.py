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
    HOME = os.environ.get('USERPROFILE')
SCRDIR = os.path.join(HOME,'Script')
DATDIR = os.path.join(HOME,'Work','Sentinel-1')
WRKDIR = os.path.join(HOME,'Work','SATREPS','Transplanting_date')
END = datetime.now().strftime('%Y%m%d')
SITES = ['Cihea','Bojongsoang']
SUBDIRS = ['Cihea:sigma0','Bojongsoang:sigma0_speckle']
OFFSETS = ['Cihea:-9.0','Bojongsoang:0.0']
VERSIONS = ['Cihea:v1.4','Bojongsoang:v1.0']

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('--scrdir',default=SCRDIR,help='Script directory (%default)')
parser.add_option('--datdir',default=DATDIR,help='Data directory (%default)')
parser.add_option('--wrkdir',default=WRKDIR,help='Work directory (%default)')
parser.add_option('-s','--str',default=None,help='Start date of estimation in the format YYYYMMDD (%default)')
parser.add_option('-e','--end',default=END,help='End date of estimation in the format YYYYMMDD (%default)')
parser.add_option('-S','--sites',default=None,action='append',help='Target sites ({})'.format(SITES))
parser.add_option('--subdirs',default=None,action='append',help='Sub data directory, for example, Cihea:sigma0 ({})'.format(SUBDIRS))
parser.add_option('--offsets',default=None,action='append',help='Offset of transplanting date, for example, Cihea:-9.0 ({})'.format(OFFSETS))
parser.add_option('--versions',default=None,action='append',help='Version of transplanting estimation, for example, Cihea:v1.0 ({})'.format(VERSIONS))
parser.add_option('--test',default=False,action='store_true',help='Test mode (%default)')
parser.add_option('-d','--debug',default=False,action='store_true',help='Debug mode (%default)')
(opts,args) = parser.parse_args()
if opts.str is None:
    opts.str = opts.end
if opts.sites is None:
    opts.sites = SITES
if opts.subdirs is None:
    opts.subdirs = SUBDIRS
if opts.offsets is None:
    opts.offsets = OFFSETS
if opts.versions is None:
    opts.versions = VERSIONS
subdir = {}
for s in opts.subdirs:
    m = re.search('([^:]+):([^:]+)',s)
    if not m:
        raise ValueError('Error in subdir >>> '+s)
    subdir.update({m.group(1).lower():m.group(2)})
offset = {}
for s in opts.offsets:
    m = re.search('([^:]+):([^:]+)',s)
    if not m:
        raise ValueError('Error in offset >>> '+s)
    offset.update({m.group(1).lower():float(m.group(2))})
version = {}
for s in opts.versions:
    m = re.search('([^:]+):([^:]+)',s)
    if not m:
        raise ValueError('Error in version >>> '+s)
    version.update({m.group(1).lower():m.group(2)})

d1 = datetime.strptime(opts.str,'%Y%m%d')
d2 = datetime.strptime(opts.end,'%Y%m%d')
d3 = d1+relativedelta(months=1)
if (d1+timedelta(days=1)).month != d3.month: # not the end of month
    d1 = datetime(d3.year,d3.month,1)-timedelta(days=1)
d3 = d2+relativedelta(months=1)
if (d2+timedelta(days=1)).month != d3.month: # not the end of month
    d2 = datetime(d2.year,d2.month,1)-timedelta(days=1)
if d1 > d2:
    d1 = d2

topdir = os.getcwd()
for site in opts.sites:
    site_low = site.lower()
    d = d1
    while d <= d2:
        dstr = d.strftime('%Y%m%d')
        #sys.stderr.write(dstr+'\n')
        wrkdir = os.path.join(opts.wrkdir,site,'final',version[site_low],dstr)
        if not os.path.exists(wrkdir):
            os.makedirs(wrkdir)
        if not os.path.isdir(wrkdir):
            raise IOError('Error, no such directory >>> '+wrkdir)
        os.chdir(wrkdir)
        file_list = []
        try:
            if site_low == 'cihea':
                t1 = datetime(d.year,d.month,1) # the first day of the month
                tmin = (t1-relativedelta(months=6)).strftime('%Y%m%d')
                tmax = (t1-relativedelta(months=3)-timedelta(days=1)).strftime('%Y%m%d')
                data_tmin = (t1-relativedelta(months=9)).strftime('%Y%m%d')
                data_tmax = dstr
                jsn_fnam = os.path.join(wrkdir,'trans_date_{}_{}_final.json'.format(site_low,dstr))
                tif_fnam = os.path.join(wrkdir,'trans_date_{}_{}_final.tif'.format(site_low,dstr))
                shp_bnam = os.path.join(wrkdir,'trans_date_{}_{}_final'.format(site_low,dstr))
                shp_fnam = shp_bnam+'.shp'
                trans_pixel_image = os.path.join(wrkdir,'trans_pixel_{}_{}_final.pdf'.format(site_low,dstr))
                trans_field_image = os.path.join(wrkdir,'trans_field_{}_{}_final.pdf'.format(site_low,dstr))
                command = 'python'
                command += ' '+os.path.join(opts.scrdir,'calc_trans_date.py')
                #command += ' -x 600 -X 900 -y 680 -Y 1030'
                command += ' --tmin '+tmin
                command += ' --tmax '+tmax
                command += ' --data_tmin '+data_tmin
                command += ' --data_tmax '+data_tmax
                command += ' --offset {:.4f}'.format(offset[site_low])
                command += ' --incidence_list '+os.path.join(opts.wrkdir,site,'incidence_list.dat')
                command += ' --datdir '+os.path.join(opts.datdir,site,subdir[site_low])
                command += ' --search_key resample'
                command += ' --near_fnam '+os.path.join(opts.wrkdir,site,'find_nearest.npz')
                command += ' --mask_fnam '+os.path.join(opts.wrkdir,site,'paddy_mask.tif')
                command += ' --json_fnam '+jsn_fnam
                command += ' --out_fnam '+tif_fnam
                command += ' 2>'+os.path.join(wrkdir,'err')
                sys.stderr.write(command+'\n')
                call(command,shell=True)
                if os.path.exists(tif_fnam):
                    file_list.append(jsn_fnam)
                    file_list.append(tif_fnam)
                    command = 'python'
                    command += ' '+os.path.join(opts.scrdir,'field_mean_shapefile.py')
                    command += ' --data_fnam '+tif_fnam
                    command += ' --area_fnam '+os.path.join(opts.wrkdir,site,'pixel_area_block.dat')
                    command += ' --outnam '+shp_bnam
                    command += ' --use_index'
                    call(command,shell=True)
                    command = 'python'
                    command += ' '+os.path.join(opts.scrdir,'draw_trans_pixel_{}.py'.format(site_low))
                    command += ' --tmin '+tmin
                    command += ' --tmax '+tmax
                    command += ' --smin 0'
                    command += ' --smax 6'
                    command += ' --title "Search Period: {} - {}"'.format(tmin,tmax)
                    command += ' --trans_fnam '+tif_fnam
                    command += ' --mask_fnam '+os.path.join(opts.wrkdir,site,'paddy_mask_studyarea.tif')
                    command += ' --output_fnam '+trans_pixel_image
                    command += ' --add_tmax'
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
                    command += ' '+os.path.join(opts.scrdir,'draw_trans_field_{}.py'.format(site_low))
                    command += ' --tmin '+tmin
                    command += ' --tmax '+tmax
                    command += ' --smin 0'
                    command += ' --smax 6'
                    command += ' --title "Search Period: {} - {}"'.format(tmin,tmax)
                    command += ' --trans_fnam '+shp_fnam
                    command += ' --mask_fnam '+os.path.join(opts.wrkdir,site,'paddy_mask_studyarea.dat')
                    command += ' --output_fnam '+trans_field_image
                    command += ' --use_index'
                    command += ' --add_tmax'
                    command += ' --batch'
                    call(command,shell=True)
                    if os.path.exists(trans_field_image):
                        file_list.append(trans_field_image)
            elif site_low == 'bojongsoang':
                t1 = datetime(d.year,d.month,1) # the first day of the month
                tmin = (t1-relativedelta(months=5)).strftime('%Y%m%d')
                tmax = (t1-relativedelta(months=2)-timedelta(days=1)).strftime('%Y%m%d')
                data_tmin = (t1-relativedelta(months=11)).strftime('%Y%m%d')
                data_tmax = dstr
                jsn_fnam = os.path.join(wrkdir,'trans_date_{}_{}_final.json'.format(site_low,dstr))
                shp_bnam = os.path.join(wrkdir,'trans_date_{}_{}_final'.format(site_low,dstr))
                shp_fnam = shp_bnam+'.shp'
                fpi_s_field_image = os.path.join(wrkdir,'fpi_s_field_{}_{}_final.pdf'.format(site_low,dstr))
                trans_field_image = []
                for ican in range(3):
                    trans_field_image.append(os.path.join(wrkdir,'trans_field_{}_{}_{}_final.pdf'.format(site_low,dstr,ican+1)))
                command = 'python'
                command += ' '+os.path.join(opts.scrdir,'calc_trans_date_fi.py')
                command += ' --tmin '+tmin
                command += ' --tmax '+tmax
                command += ' --data_tmin '+data_tmin
                command += ' --data_tmax '+data_tmax
                command += ' --offset {:.4f}'.format(offset[site_low])
                command += ' --incidence_list '+os.path.join(opts.wrkdir,site,'incidence_list.dat')
                command += ' --datdir '+os.path.join(opts.datdir,site,subdir[site_low])
                command += ' --search_key resample'
                command += ' --x_profile '+os.path.join(opts.wrkdir,site,'x_profile.npy')
                command += ' --y_profile '+os.path.join(opts.wrkdir,site,'y_profile.npy')
                command += ' --area_fnam '+os.path.join(opts.wrkdir,site,'pixel_area_block.dat')
                command += ' --shp_fnam '+os.path.join(opts.wrkdir,site,site)
                command += ' --json_fnam '+jsn_fnam
                command += ' --out_fnam '+shp_bnam
                command += ' --select_post_s'
                command += ' 2>'+os.path.join(wrkdir,'err')
                sys.stderr.write(command+'\n')
                call(command,shell=True)
                if os.path.exists(shp_fnam):
                    file_list.append(jsn_fnam)
                    file_list.append(shp_fnam)
                    file_list.append(shp_bnam+'.dbf')
                    file_list.append(shp_bnam+'.prj')
                    file_list.append(shp_bnam+'.shx')
                    command = 'python'
                    command += ' '+os.path.join(opts.scrdir,'draw_fpi_field_{}.py'.format(site_low))
                    command += ' --title {}'.format(tmin)
                    command += ' --trans_fnam '+shp_fnam
                    command += ' --output_fnam '+fpi_s_field_image
                    command += ' --fpi_s'
                    command += ' --batch'
                    call(command,shell=True)
                    if os.path.exists(fpi_s_field_image):
                        file_list.append(fpi_s_field_image)
                    for ican in range(3):
                        command = 'python'
                        command += ' '+os.path.join(opts.scrdir,'draw_trans_field_{}.py'.format(site_low))
                        command += ' --tmin '+tmin
                        command += ' --tmax '+tmax
                        command += ' --pmin 0'
                        command += ' --pmax 5'
                        command += ' --ncan {}'.format(ican+1)
                        command += ' --title "#{}, Search Period: {} - {}"'.format(ican+1,tmin,tmax)
                        command += ' --trans_fnam '+shp_fnam
                        command += ' --output_fnam '+trans_field_image[ican]
                        command += ' --add_tmax'
                        command += ' --batch'
                        call(command,shell=True)
                        if os.path.exists(trans_field_image[ican]):
                            file_list.append(trans_field_image[ican])
            else:
                raise ValueError('Error, unknown site >>> '+site)
            if len(file_list) > 0:
                command = 'python'
                command += ' '+os.path.join(opts.scrdir,'trans_date_upload_nas.py')
                command += ' --site '+site
                command += ' --date '+dstr
                if opts.test:
                    command += ' --level test'
                else:
                    command += ' --level final'
                command += ' --version '+version[site_low]
                command += ' --overwrite'
                command += ' '+' '.join(file_list)
                call(command,shell=True)
                command = 'python'
                command += ' '+os.path.join(opts.scrdir,'trans_date_copy.py')
                command += ' --site '+site
                command += ' --date '+dstr
                if opts.test:
                    command += ' --level test'
                else:
                    command += ' --level final'
                command += ' --version '+version[site_low]
                command += ' --overwrite'
                command += ' '+' '.join(file_list)
                call(command,shell=True)
        except Exception:
            sys.stderr.write('Warning, error occured.\n')
        d3 = d+relativedelta(months=2)
        d = datetime(d3.year,d3.month,1)-timedelta(days=1)
        os.chdir(topdir)
