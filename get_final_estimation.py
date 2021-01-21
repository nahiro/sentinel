#!/usr/bin/env python
import os
import sys
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
parser.add_option('-d','--debug',default=False,action='store_true',help='Debug mode (%default)')
(opts,args) = parser.parse_args()
if opts.str is None:
    opts.str = opts.end
if opts.sites is None:
    opts.sites = SITES

d1 = datetime.strptime(opts.str,'%Y%m%d')
d2 = datetime.strptime(opts.end,'%Y%m%d')
d3 = d2+relativedelta(months=1)
if (d2+timedelta(days=1)).month != d3.month: # not the end of month
    d2 = datetime(d2.year,d2.month,1)-timedelta(days=1)
if d1 > d2:
    d1 = d2

topdir = os.getcwd()
d = d1
while d <= d2:
    dstr = d.strftime('%Y%m%d')
    t1 = datetime(d.year,d.month,1) # the first day of the month
    tmin = (t1-relativedelta(months=3)).strftime('%Y%m%d')
    tmax = (t1-timedelta(days=1)).strftime('%Y%m%d')
    data_tmin = (t1-relativedelta(months=5)).strftime('%Y%m%d')
    data_tmax = dstr
    #sys.stderr.write(dstr+'\n')
    for site in opts.sites:
        wrkdir = os.path.join(opts.wrkdir,site,'final',dstr)
        if not os.path.exists(wrkdir):
            os.makedirs(wrkdir)
        if not os.path.isdir(wrkdir):
            raise IOError('Error, no such directory >>> '+wrkdir)
        os.chdir(wrkdir)
        command = 'python'
        command += ' '+os.path.join(opts.scrdir,'calc_trans_date.py')
        command += ' --tmin '+tmin
        command += ' --tmax '+tmax
        command += ' --data_tmin '+data_tmin
        command += ' --data_tmax '+data_tmax
        command += ' -D '+os.path.join(opts.datdir,site,'sigma0_speckle')
        command += ' --search_key resample'
        command += ' --near_fnam '+os.path.join(opts.wrkdir,site,'find_nearest.npz')
        command += ' --incidence_list '+os.path.join(opts.wrkdir,site,'incidence_list.dat')
        command += ' 2>'+os.path.join(wrkdir,'err')
        sys.stderr.write(command+'\n')
        #call(command,shell=True)
        os.chdir(topdir)
    d3 = d+relativedelta(months=2)
    d = datetime(d3.year,d3.month,1)-timedelta(days=1)


#./calc_trans_date.py --tmin 20190315 --tmax 20190615 --data_tmin 20190101 --data_tmax 20190831 -D /home/naohiro/Work/Sentinel-1/Cihea/sigma0_speckle --search_key resample --incidence_list incidence_list.dat -x 600 -X 900 -y 680 -Y 1030 2>err
#./field_mean_shapefile.py
