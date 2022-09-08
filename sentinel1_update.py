#!/usr/bin/env python
import os
import sys
import re
import time
from datetime import datetime,timedelta
import numpy as np
from subprocess import call
from optparse import OptionParser,IndentedHelpFormatter

# Default values
HOME = os.environ.get('HOME')
if HOME is None:
    HOME = os.environ.get('USERPROFILE')
SCRDIR = os.path.join(HOME,'Script')
DATDIR = os.path.join(HOME,'Work','Sentinel-1_Data')
END = datetime.now().strftime('%Y%m%d')
SITES = ['Cihea','Bojongsoang']

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('--scrdir',default=SCRDIR,help='Script directory (%default)')
parser.add_option('--datdir',default=DATDIR,help='Data directory (%default)')
parser.add_option('-s','--str',default=None,help='Start date of download in the format YYYYMMDD (%default)')
parser.add_option('-e','--end',default=END,help='End date of download in the format YYYYMMDD (%default)')
parser.add_option('-S','--sites',default=None,action='append',help='Target sites ({})'.format(SITES))
parser.add_option('-M','--max_retry',default=None,type='int',help='Maximum number of retries to download data (%default)')
parser.add_option('--skip_download',default=False,action='store_true',help='Skip download (%default)')
parser.add_option('--skip_upload',default=False,action='store_true',help='Skip upload (%default)')
parser.add_option('--skip_copy',default=False,action='store_true',help='Skip copy (%default)')
parser.add_option('-d','--debug',default=False,action='store_true',help='Debug mode (%default)')
(opts,args) = parser.parse_args()
if opts.sites is None:
    opts.sites = SITES

# Determin the start date of download
dmaxs = []
if opts.str is not None:
    for site in opts.sites:
        dmaxs.append(opts.str)
else:
    for site in opts.sites:
        datdir = os.path.join(opts.datdir,site,'GRD')
        if not os.path.exists(datdir):
            os.makedirs(datdir)
        if not os.path.isdir(datdir):
            raise IOError('Error, no such folder >>> {}'.format(datdir))
        dmax = '0'*8
        for d in sorted(os.listdir(datdir)):
            if not re.search('^\d\d\d\d$',d):
                continue
            dnam = os.path.join(datdir,d)
            if not os.path.isdir(dnam):
                continue
            for f in sorted(os.listdir(dnam)):
                #S1A_IW_GRDH_1SDV_20171227T223338_20171227T223405_019894_021DC8_434F.zip
                #S1B_IW_GRDH_1SDV_20200116T223300_20200116T223336_019848_025883_2DEF.zip
                m = re.search('^S1[AB]_IW_GRDH_1SDV_('+'\d'*8+')T\S+\.zip$',f)
                if not m:
                    continue
                dstr = m.group(1)
                if dstr > dmax:
                    dmax = dstr
        if dmax == '0'*8:
            sys.stderr.write('Warning, failed in determining the start date of download >>> {}\n'.format(site))
            sys.stderr.flush()
            dmax = (datetime.strptime(opts.end,'%Y%m%d')-timedelta(days=2)).strftime('%Y%m%d')
        dmaxs.append((datetime.strptime(dmax,'%Y%m%d')+timedelta(days=1)).strftime('%Y%m%d'))
if len(dmaxs) != len(opts.sites):
    raise ValueError('Error, len(dmaxs)={}, len(opts.sites)={}'.format(len(dmaxs),len(opts.sites)))

# Download data
topdir = os.getcwd()
gnams = {}
for site,start in zip(opts.sites,dmaxs):
    site_low = site.lower()
    datdir = os.path.join(opts.datdir,site,'GRD')
    fnam = os.path.join(opts.datdir,site,site_low+'.json')
    if not os.path.exists(fnam):
        raise IOError('No such file >>> '+fnam)
    gnams.update({site:[]})
    command = 'python'
    command += ' '+os.path.join(opts.scrdir,'sentinel_wget.py')
    command += ' -U https://scihub.copernicus.eu/dhus'
    command += ' --netrc "{}"'.format(os.path.join(HOME,'.netrc'))
    command += ' --machine scihub.copernicus.eu'
    command += ' --geometry '+fnam
    command += ' --log '+os.path.join(datdir,site_low+'.log')
    command += ' --producttyp GRD'
    command += ' --start '+start
    command += ' --end '+opts.end
    command += ' --path '+datdir
    command += ' --download'
    command += ' --sort_year'
    if opts.max_retry is not None:
        command += ' --max_retry {}'.format(opts.max_retry)
    if not opts.skip_download:
        sys.stderr.write(command+'\n')
        sys.stderr.flush()
        call(command,shell=True)
    d1 = datetime.strptime(start,'%Y%m%d')
    d2 = datetime.strptime(opts.end,'%Y%m%d')
    for y in range(d1.year,d2.year+1):
        year = '{:04d}'.format(y)
        dnam = os.path.join(datdir,year)
        if not os.path.isdir(dnam):
            continue
        for f in sorted(os.listdir(dnam)):
            #S1A_IW_GRDH_1SDV_20171227T223338_20171227T223405_019894_021DC8_434F.zip
            #S1B_IW_GRDH_1SDV_20200116T223300_20200116T223336_019848_025883_2DEF.zip
            m = re.search('^S1[AB]_IW_GRDH_1SDV_('+'\d'*8+')'+'T\S+\.zip$',f)
            if not m:
                continue
            dstr = m.group(1)
            dtim = datetime.strptime(dstr,'%Y%m%d')
            if (dtim < d1) or (dtim > d2):
                continue
            gnam = os.path.join(dnam,f)
            gnams[site].append(gnam)

# Upload data
if not opts.skip_upload:
    for site in opts.sites:
        for gnam in gnams[site]:
            command = 'python'
            command += ' '+os.path.join(opts.scrdir,'sentinel1_upload_nas.py')
            command += ' --site '+site
            command += ' '+gnam
            call(command,shell=True)
    os.chdir(topdir)

# Copy data
if not opts.skip_copy:
    for site in opts.sites:
        for gnam in gnams[site]:
            command = 'python'
            command += ' '+os.path.join(opts.scrdir,'sentinel1_copy.py')
            command += ' --site '+site
            command += ' '+gnam
            call(command,shell=True)
