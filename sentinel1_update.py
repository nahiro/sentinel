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
SCRDIR = '/home/naohiro/Script'
DATDIR = '/home/naohiro/Work/Sentinel-1'
END = datetime.now().strftime('%Y%m%d')
SITES = ['Cihea','Bojongsoang']

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.set_usage('Usage: %prog collocated_geotiff_file [options]')
parser.add_option('--scrdir',default=SCRDIR,help='Script directory (%default)')
parser.add_option('--datdir',default=DATDIR,help='Data directory (%default)')
parser.add_option('-s','--str',default=None,help='Start date of download in the format YYYYMMDD (%default)')
parser.add_option('-e','--end',default=END,help='End date of download in the format YYYYMMDD (%default)')
parser.add_option('-S','--sites',default=None,action='append',help='Target sites ({})'.format(SITES))
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
        fnam = os.path.join(opts.scrdir,site.lower()+'.json')
        if not os.path.exists(fnam):
            raise IOError('No such file >>> '+fnam)
        datdir = os.path.join(opts.datdir,site)
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
            raise ValueError('Error in determining the start date of download >>> '+site)
        dmaxs.append((datetime.strptime(dmax,'%Y%m%d')+timedelta(days=1)).strftime('%Y%m%d'))
if len(dmaxs) != len(opts.sites):
    raise ValueError('Error, len(dmaxs)={}, len(opts.sites)={}'.format(len(dmaxs),len(opts.sites)))

# Download data
topdir = os.getcwd()
for site,start in zip(opts.sites,dmaxs):
    datdir = os.path.join(opts.datdir,site)
    os.chdir(datdir)
    command = 'python'
    command += ' '+os.path.join(opts.scrdir,'sentinel_download.py')
    command += ' -g '+os.path.join(opts.scrdir,site.lower()+'.json')
    command += ' -t GRD'
    command += ' -s '+start
    command += ' -e '+opts.end
    command += ' -d'
    print(command)
    call(command,shell=True)
    os.chdir(topdir)
    for f in sorted(os.listdir(datdir)):
        #S1A_IW_GRDH_1SDV_20171227T223338_20171227T223405_019894_021DC8_434F.zip
        #S1B_IW_GRDH_1SDV_20200116T223300_20200116T223336_019848_025883_2DEF.zip
        m = re.search('^S1[AB]_IW_GRDH_1SDV_('+'\d'*4+')'+'\d'*4+'T\S+\.zip$',f)
        if not m:
            continue
        year = m.group(1)
        fnam = os.path.join(datdir,f)
        gnam = os.path.join(datdir,year,f)
        os.rename(fnam,gnam)
