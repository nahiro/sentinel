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
OFFSETS = ['Cihea:-9.0','Bojongsoang:0.0']

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('--scrdir',default=SCRDIR,help='Script directory (%default)')
parser.add_option('--datdir',default=DATDIR,help='Data directory (%default)')
parser.add_option('--wrkdir',default=WRKDIR,help='Work directory (%default)')
parser.add_option('-s','--str',default=None,help='Start date of update in the format YYYYMMDD (%default)')
parser.add_option('-e','--end',default=END,help='End date of update in the format YYYYMMDD (%default)')
parser.add_option('-S','--sites',default=None,action='append',help='Target sites ({})'.format(SITES))
parser.add_option('--offsets',default=None,action='append',help='Offset of transplanting date, for example, Cihea:-9.0 ({})'.format(OFFSETS))
parser.add_option('--test',default=False,action='store_true',help='Test mode (%default)')
parser.add_option('-d','--debug',default=False,action='store_true',help='Debug mode (%default)')
(opts,args) = parser.parse_args()
if opts.sites is None:
    opts.sites = SITES
if opts.offsets is None:
    opts.offsets = OFFSETS
offset = {}
for s in opts.offsets:
    m = re.search('([^:]+):([^:]+)',s)
    if not m:
        raise ValueError('Error in offset >>> '+s)
    offset.update({m.group(1):float(m.group(2))})

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
        for d1 in dtims[site]:
            if d1 > d3:
                break
    for d in dtims[site]:
        if d < d1:
            continue
        if d > d2:
            break
        dstr = d.strftime('%Y%m%d')
        #sys.stderr.write(dstr+'\n')
        wrkdir = os.path.join(opts.wrkdir,site,'preliminary',dstr)
        if not os.path.exists(wrkdir):
            continue
        if not os.path.isdir(wrkdir):
            raise IOError('Error, no such directory >>> '+wrkdir)
        os.chdir(wrkdir)
        file_list = []
        try:
            if site.lower() == 'cihea':
                tmin = (d-timedelta(days=90)).strftime('%Y%m%d')
                tmax = dstr
                data_tmin = (d-timedelta(days=150)).strftime('%Y%m%d')
                data_tmax = dstr
                jsn_fnam = os.path.join(wrkdir,'trans_date_{}_{}_preliminary.json'.format(site,dstr))
                tif_fnam = os.path.join(wrkdir,'trans_date_{}_{}_preliminary.tif'.format(site,dstr))
                shp_bnam = os.path.join(wrkdir,'trans_date_{}_{}_preliminary'.format(site,dstr))
                shp_fnam = shp_bnam+'.shp'
                trans_pixel_image = os.path.join(wrkdir,'trans_pixel_{}_{}_preliminary.pdf'.format(site,dstr))
                trans_field_image = os.path.join(wrkdir,'trans_field_{}_{}_preliminary.pdf'.format(site,dstr))
                command = 'python'
                command += ' '+os.path.join(opts.scrdir,'calc_trans_date.py')
                #command += ' -x 600 -X 900 -y 680 -Y 1030'
                command += ' --tmin '+tmin
                command += ' --tmax '+tmax
                command += ' --data_tmin '+data_tmin
                command += ' --data_tmax '+data_tmax
                command += ' --offset {:.4f}'.format(offset[site])
                command += ' --incidence_list '+os.path.join(opts.wrkdir,site,'incidence_list.dat')
                command += ' --datdir '+os.path.join(opts.datdir,site,'sigma0_speckle')
                command += ' --search_key resample'
                command += ' --near_fnam '+os.path.join(opts.wrkdir,site,'find_nearest.npz')
                command += ' --json_fnam '+jsn_fnam
                command += ' --out_fnam '+tif_fnam
                command += ' --early'
                command += ' 2>'+os.path.join(wrkdir,'err')
                #sys.stderr.write(command+'\n')
                #call(command,shell=True)
                if os.path.exists(tif_fnam):
                    #file_list.append(jsn_fnam)
                    #file_list.append(tif_fnam)
                    command = 'python'
                    command += ' '+os.path.join(opts.scrdir,'field_mean_shapefile.py')
                    command += ' --data_fnam '+tif_fnam
                    command += ' --area_fnam '+os.path.join(opts.wrkdir,site,'pixel_area_block.dat')
                    command += ' --outnam '+shp_bnam
                    #call(command,shell=True)
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
                    command += ' --early'
                    command += ' --batch'
                    call(command,shell=True)
                    if os.path.exists(trans_pixel_image):
                        file_list.append(trans_pixel_image)
                if os.path.exists(shp_fnam):
                    #file_list.append(shp_fnam)
                    #file_list.append(shp_bnam+'.dbf')
                    #file_list.append(shp_bnam+'.prj')
                    #file_list.append(shp_bnam+'.shx')
                    command = 'python'
                    command += ' '+os.path.join(opts.scrdir,'draw_trans_field_{}.py'.format(site.lower()))
                    command += ' --tmin '+tmin
                    command += ' --tmax '+tmax
                    command += ' --pmin 0'
                    command += ' --pmax 300'
                    command += ' --title "Search Period: {} - {}"'.format(tmin,tmax)
                    command += ' --trans_fnam '+shp_fnam
                    command += ' --output_fnam '+trans_field_image
                    command += ' --early'
                    command += ' --batch'
                    call(command,shell=True)
                    if os.path.exists(trans_field_image):
                        file_list.append(trans_field_image)
            elif site.lower() == 'bojongsoang':
                tmin = (d-timedelta(days=90)).strftime('%Y%m%d')
                tmax = dstr
                data_tmin = (d-relativedelta(years=1)+timedelta(days=1)).strftime('%Y%m%d')
                data_tmax = dstr
                jsn_fnam = os.path.join(wrkdir,'trans_date_{}_{}_preliminary.json'.format(site,dstr))
                shp_bnam = os.path.join(wrkdir,'trans_date_{}_{}_preliminary'.format(site,dstr))
                shp_fnam = shp_bnam+'.shp'
                fpi_s_field_image = os.path.join(wrkdir,'fpi_s_field_{}_{}_preliminary.pdf'.format(site,dstr))
                fpi_e_field_image = os.path.join(wrkdir,'fpi_e_field_{}_{}_preliminary.pdf'.format(site,dstr))
                trans_field_image = []
                for ican in range(3):
                    trans_field_image.append(os.path.join(wrkdir,'trans_field_{}_{}_{}_preliminary.pdf'.format(site,dstr,ican+1)))
                command = 'python'
                command += ' '+os.path.join(opts.scrdir,'calc_trans_date_fi.py')
                command += ' --tmin '+tmin
                command += ' --tmax '+tmax
                command += ' --data_tmin '+data_tmin
                command += ' --data_tmax '+data_tmax
                command += ' --offset {:.4f}'.format(offset[site])
                command += ' --incidence_list '+os.path.join(opts.wrkdir,site,'incidence_list.dat')
                command += ' --incidence_angle '+os.path.join(opts.wrkdir,site,'incidence_angle.dat')
                command += ' --datdir '+os.path.join(opts.datdir,site,'sigma0_speckle')
                command += ' --search_key resample'
                command += ' --x_profile '+os.path.join(opts.wrkdir,site,'x_profile.npy')
                command += ' --y_profile '+os.path.join(opts.wrkdir,site,'y_profile.npy')
                command += ' --area_fnam '+os.path.join(opts.wrkdir,site,'pixel_area_block.dat')
                command += ' --shp_fnam '+os.path.join(opts.wrkdir,site,site)
                command += ' --json_fnam '+jsn_fnam
                command += ' --out_fnam '+shp_bnam
                command += ' --early'
                command += ' 2>'+os.path.join(wrkdir,'err')
                #sys.stderr.write(command+'\n')
                #call(command,shell=True)
                if os.path.exists(shp_fnam):
                    #file_list.append(shp_fnam)
                    #file_list.append(shp_bnam+'.dbf')
                    #file_list.append(shp_bnam+'.prj')
                    #file_list.append(shp_bnam+'.shx')
                    command = 'python'
                    command += ' '+os.path.join(opts.scrdir,'draw_fpi_field_{}.py'.format(site.lower()))
                    command += ' --title {}'.format(tmin)
                    command += ' --trans_fnam '+shp_fnam
                    command += ' --output_fnam '+fpi_s_field_image
                    command += ' --fpi_s'
                    command += ' --batch'
                    call(command,shell=True)
                    if os.path.exists(fpi_s_field_image):
                        file_list.append(fpi_s_field_image)
                    command = 'python'
                    command += ' '+os.path.join(opts.scrdir,'draw_fpi_field_{}.py'.format(site.lower()))
                    command += ' --title {}'.format(tmax)
                    command += ' --trans_fnam '+shp_fnam
                    command += ' --output_fnam '+fpi_e_field_image
                    command += ' --batch'
                    call(command,shell=True)
                    if os.path.exists(fpi_e_field_image):
                        file_list.append(fpi_e_field_image)
                    for ican in range(3):
                        command = 'python'
                        command += ' '+os.path.join(opts.scrdir,'draw_trans_field_{}.py'.format(site.lower()))
                        command += ' --tmin '+tmin
                        command += ' --tmax '+tmax
                        command += ' --pmin 0'
                        command += ' --pmax 5'
                        command += ' --ncan {}'.format(ican+1)
                        command += ' --title "#{}, Search Period: {} - {}"'.format(ican+1,tmin,tmax)
                        command += ' --trans_fnam '+shp_fnam
                        command += ' --output_fnam '+trans_field_image[ican]
                        command += ' --early'
                        command += ' --batch'
                        call(command,shell=True)
                        if os.path.exists(trans_field_image[ican]):
                            file_list.append(trans_field_image[ican])
            else:
                raise ValueError('Error, unknown site >>> '+site)
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
        os.chdir(topdir)
