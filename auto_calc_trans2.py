#!/usr/bin/env python
import os
import sys
import shutil
import re
from datetime import datetime,timedelta
from subprocess import call
from optparse import OptionParser,IndentedHelpFormatter

# Default values
HOME = os.environ.get('HOME')
if HOME is None:
    HOME = os.environ.get('HOMEPATH')
SCRDIR = os.path.join(HOME,'Script')
DATDIR = os.path.join(HOME,'Work','Sentinel-2','L2A')
WRKDIR = os.path.join(HOME,'Work','SATREPS','Transplanting_date')
SITES = ['Bojongsoang']
DATE_FINAL = 5

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('--scrdir',default=SCRDIR,help='Script directory (%default)')
parser.add_option('--datdir',default=DATDIR,help='Data directory (%default)')
parser.add_option('--wrkdir',default=WRKDIR,help='Work directory (%default)')
parser.add_option('-s','--str',default=None,help='Start date in the format YYYYMMDD (%default)')
parser.add_option('-e','--end',default=None,help='End date in the format YYYYMMDD (%default)')
parser.add_option('-S','--sites',default=None,action='append',help='Target sites ({})'.format(SITES))
parser.add_option('--date_final',default=DATE_FINAL,type='int',help='Date to calculate final estimation (%default)')
parser.add_option('--skip_upload',default=False,action='store_true',help='Skip upload (%default)')
parser.add_option('-d','--debug',default=False,action='store_true',help='Debug mode (%default)')
(opts,args) = parser.parse_args()
if opts.sites is None:
    opts.sites = SITES

for site in opts.sites:
    datdir = os.path.join(opts.datdir,site)
    # Download/Upload data
    log = os.path.join(datdir,site.lower()+'.log')
    command = 'python'
    command += ' '+os.path.join(opts.scrdir,'sentinel2_update.py')
    command += ' --scrdir '+opts.scrdir
    command += ' --datdir '+opts.datdir
    if opts.str is not None:
        command += ' --str '+opts.str
    if opts.end is not None:
        command += ' --end '+opts.end
    if opts.skip_upload:
        command += ' --skip_upload'
    command += ' --sites '+site
    #sys.stderr.write(command+'\n')
    call(command,shell=True)
    fnams = []
    dstrs = []
    if os.path.exists(log):
        with open(log,'r') as fp:
            for line in fp:
                item = line.split()
                if len(item) < 1:
                    continue
                fnam = item[0]
                f = os.path.basename(fnam)
                #S2A_MSIL2A_20210713T025551_N0301_R032_T48MYT_20210713T063300.zip
                #S2B_MSIL2A_20210715T024549_N0301_R132_T48MYT_20210715T055013.zip
                m = re.search('^S2[AB]_MSIL2A_('+'\d'*8+')T\S+\.zip$',f)
                if not m:
                    continue
                fnams.append(fnam)
                dstrs.append(m.group(1))
    if len(dstrs) < 1:
        continue
    # Subset
    subset_dstrs = []
    for fnam,dstr in zip(fnams,dstrs):
        gnam = os.path.join(datdir,'subset',dstr+'.tif')
        if not os.path.exists(gnam):
            unzip_flag = False
            rnam = os.path.splitext(fnam)[0]+'.SAFE'
            if not os.path.exists(rnam):
                try:
                    with zipfile.ZipFile(fnam,'r') as z:
                        z.extractall(dnam)
                    unzip_flag = True
                except Exception:
                    continue
            command = 'python'
            command += ' '+os.path.join(opts.scrdir,'sentinel2_subset.py')
            command += ' '+rnam
            command += ' --datdir '+os.path.join(datdir,'subset')
            command += ' --site '+site
            command += ' --resolution 10.0'
            command += ' --geotiff'
            call(command,shell=True)
            if unzip_flag:
                shutil.rmtree(rnam)
            # Remove cache
            command = 'python'
            command += ' '+os.path.join(opts.scrdir,'remove_snap_cache.py')
            call(command,shell=True)
        if os.path.exists(gnam):
            subset_dstrs.append(dstr)
    # Geocor
    geocor_dstrs = []
    for dstr in subset_dstrs:
        fnam = os.path.join(datdir,'subset',dstr+'.tif')
        gnam = os.path.join(datdir,'geocor',dstr+'_geocor.tif')
        ref_fnam = os.path.join(opts.wrkdir,site,'geocor_reference.tif')
        dat_fnam = os.path.join(datdir,'geocor',dstr+'_geocor.dat')
        sel_fnam = os.path.join(datdir,'geocor',dstr+'_geocor_selected.dat')
        if not os.path.exists(gnam):
            command = 'python'
            command += ' '+os.path.join(opts.scrdir,'find_gcps.py')
            command += ' '+fnam
            command += ' '+ref_fnam
            command += ' --ref_band {}' # Red band
            command += ' --trg_band {}' # Red band
            command += ' --ref_data_min {}' # for WorldView data
            command += ' --exp'
            try:
                out = check_output(command,shell=True).decode()
                with open(dat_fnam,'w') as fp:
                    fp.write(out)
            except Exception:
                continue
            command = 'python'
            command += ' '+os.path.join(opts.scrdir,'select_gcps.py')
            command += ' '+dat_fnam
            command += ' --datdir '+os.path.dirname(dat_fnam)
            command += ' --replace'
            command += ' --exp'
            try:
                call(command,shell=True)
            except Exception:
                continue
            command = 'python'
            command += ' '+os.path.join(opts.scrdir,'auto_geocor.py')
            command += ' '+fnam
            command += ' '+ref_fnam
            command += ' --ref_band {}' # Red band
            command += ' --trg_band {}' # Red band
            command += ' --out_fnam '+gnam
            command += ' --tr 10.0'
            command += ' --use_gcps '+sel_fnam # use
            try:
                call(command,shell=True)
            except Exception:
                continue
        if os.path.exists(sel_fnam):
            os.rename(sel_fnam,dat_fnam)
        if os.path.exists(gnam):
            geocor_dstrs.append(dstr)
    # Resample
    resample_dstrs = []
    for dstr in geocor_dstrs:
        fnam = os.path.join(datdir,'geocor',dstr+'_geocor.tif')
        gnam = os.path.join(datdir,'resample',dstr+'_geocor_resample.tif')
        if not os.path.exists(gnam):
            command = 'python'
            command += ' '+os.path.join(opts.scrdir,'sentinel_resample.py')
            command += ' '+fnam
            command += ' --datdir '+os.path.join(datdir,'resample')
            command += ' --site '+site
            command += ' --read_comments'
            command += ' --band_fnam '+os.path.join(wrkdir,site,'band_names.txt')
            try:
                call(command,shell=True)
            except Exception:
                continue
        if os.path.exists(gnam):
            resample_dstrs.append(dstr)
    # Atcor
    atcor_dstrs = []
    for dstr in resample_dstrs:
        fnam = os.path.join(datdir,'resample',dstr+'_geocor_resample.tif')
        gnam = os.path.join(datdir,'atcor',dstr+'_ndvi_correct.npz')
        fit_fnam = os.path.join(datdir,'atcor',dstr+'_ndvi_fit.npz')
        if not os.path.exists(gnam):
            command = 'python'
            command += ' '+os.path.join(scrdir,'sentinel2_atcor_fit.py')
            command += ' '+fnam
            command += ' --band ndvi'
            command += ' --mask_fnam '+os.path.join(wrkdir,site,'atcor_mask.tif')
            command += ' --stat_fnam '+os.path.join(wrkdir,site,'ndvi_stat.npz')
            command += ' --inds_fnam '+os.path.join(wrkdir,site,'nearest_inds_1000.npy')
            command += ' --output_fnam '+fit_fnam
            call(command,shell=True)
            try:
                call(command,shell=True)
            except Exception:
                continue
            command = 'python'
            command += ' '+os.path.join(scrdir,'sentinel2_atcor_correct.py')
            command += ' '+fnam
            command += ' --band ndvi'
            command += ' --area_fnam '+os.path.join(wrkdir,site,'pixel_area_block.dat')
            command += ' --param_fnam '+fit_fnam
            command += ' --output_fnam '+gnam
            try:
                call(command,shell=True)
            except Exception:
                continue
        if os.path.exists(gnam):
            atcor_dstrs.append(dstr)
    for dstr in dstrs:
        # Cflag
        dtim = datetime.strptime(dstr,'%Y%m%d')
        tmin = (dtim+timedelta(months=-6)).strftime('%Y%m%d')
        tmax = dtim.strftime('%Y%m%d')
        data_tmin = (dtim+timedelta(years=-1)).strftime('%Y%m%d')
        data_tmax = dtim.strftime('%Y%m%d')
        command = 'python'
        command += ' '+os.path.join(scrdir,'sentinel2_cflag.py')
        command += ' --tmin '+tmin
        command += ' --tmax '+tmax
        command += ' --data_tmin '+data_tmin
        command += ' --data_tmax '+data_tmax
        command += ' --datdir '+os.path.join(datdir,'atcor')
        command += ' --outdir '+os.path.join(datdir,'cflag')
        try:
            call(command,shell=True)
        except Exception:
            continue
        # Preliminary estimation
        d1 = (dtim+timedelta(days=-1)).strftime('%Y%m%d')
        d2 = (dtim+timedelta(days=+1)).strftime('%Y%m%d')
        command = 'python'
        command += ' '+os.path.join(opts.scrdir,'get_preliminary_estimation2.py')
        command += ' --scrdir '+opts.scrdir
        command += ' --datdir '+opts.datdir
        command += ' --str '+d1
        command += ' --end '+d2
        command += ' --sites '+site
        try:
            call(command,shell=True)
        except Exception:
            continue
    if os.path.exists(log):
        os.remove(log)

if datetime.now().day == opts.date_final:
    # Final estimation
    for site in opts.sites:
        command = 'python'
        command += ' '+os.path.join(opts.scrdir,'get_final_estimation2.py')
        command += ' --scrdir '+opts.scrdir
        command += ' --datdir '+opts.datdir
        command += ' --sites '+site
        call(command,shell=True)
