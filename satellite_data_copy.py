import os
import sys
import re
import tempfile
from datetime import datetime,timedelta
from glob import glob
import numpy as np
import pandas as pd
from argparse import ArgumentParser,RawTextHelpFormatter

# Constants
HOME = os.environ.get('USERPROFILE')
if HOME is None:
    HOME = os.environ.get('HOME')
TOPDIR = os.path.join(HOME,'Work')
ITEMS = ['planting','L2A','geocor','indices','parcel','atcor','interp','tentative_interp']

# Default values
PROC_NAME = 'satellite_data_copy.py'
SITES = ['Cihea','Bojongsoang','Testsite']
SERVER = 'satreps-nas'
PORT = 443
PYTHON_PATH = os.path.join(HOME,'miniconda3','python')
RCDIR = HOME
SCRDIR = os.path.join(HOME,'SatelliteTool')
S1_DATA = os.path.join(TOPDIR,'Sentinel-1_Data')
S2_DATA = os.path.join(TOPDIR,'Sentinel-2_Data')
TRANS_PATH = '/SATREPS/ipb/User/1_Spatial-information/Transplanting_date'
S1_PATH = '/SATREPS/ipb/User/1_Spatial-information/Sentinel-1'
S2_PATH = '/SATREPS/ipb/User/1_Spatial-information/Sentinel-2'
END_DATE = datetime.now().strftime('%Y%m%d')

# Read options
parser = ArgumentParser(formatter_class=lambda prog:RawTextHelpFormatter(prog,max_help_position=200,width=200))
parser.add_argument('--proc_name',default=PROC_NAME,help='Process name (%(default)s)')
parser.add_argument('-S','--sites',default=None,action='append',help='Target sites ({})'.format(SITES))
parser.add_argument('--server',default=SERVER,help='Name of the server (%(default)s)')
parser.add_argument('-P','--port',default=PORT,type=int,help='Port# of the server (%(default)s)')
parser.add_argument('--python_path',default=PYTHON_PATH,help='Path to the Python (%(default)s)')
parser.add_argument('--rcdir',default=RCDIR,help='Directory where .netrc exists (%(default)s)')
parser.add_argument('--scrdir',default=SCRDIR,help='Script folder (%(default)s)')
parser.add_argument('--s1_data',default=S1_DATA,help='Sentinel-1 data folder (%(default)s)')
parser.add_argument('--s2_data',default=S2_DATA,help='Sentinel-2 data folder (%(default)s)')
parser.add_argument('--trans_path',default=TRANS_PATH,help='Transplanting data path on NAS (%(default)s)')
parser.add_argument('--s1_path',default=S1_PATH,help='Sentinel-1 data path on NAS (%(default)s)')
parser.add_argument('--s2_path',default=S2_PATH,help='Sentinel-2 data path on NAS (%(default)s)')
parser.add_argument('-s','--start_date',default=None,help='Planting start date in the format YYYYMMDD (%(default)s)')
parser.add_argument('-e','--end_date',default=END_DATE,help='Planting end date in the format YYYYMMDD (%(default)s)')
parser.add_argument('--first_date',default=None,help='Data first date in the format YYYYMMDD (%(default)s)')
parser.add_argument('--last_date',default=None,help='Data last date in the format YYYYMMDD (%(default)s)')
parser.add_argument('--overwrite',default=False,action='store_true',help='Overwrite mode (%(default)s)')
args = parser.parse_args()
if args.sites is None:
    args.sites = SITES

def mktemp(suffix='',prefix=''):
    dnam = tempfile.gettempdir()
    string_seed = string.digits + string.ascii_lowercase + string.ascii_uppercase
    random_string = ''.join(random.choices(string_seed,k=8))
    return os.path.join(dnam,'{}{}_{}{}'.format(prefix,args.proc_name,random_string,suffix))

def run_command(command,message=None,print_command=True,print_time=True):
    if message is not None:
        sys.stderr.write('\n'+message+'\n')
        sys.stderr.flush()
    if print_command:
        sys.stderr.write('\n'+command+'\n')
        sys.stderr.flush()
    if print_time:
        t1 = datetime.now()
        sys.stderr.write('\nStart: {}\n'.format(t1))
        sys.stderr.flush()
    ret = call(command,shell=True)
    if print_time:
        t2 = datetime.now()
        sys.stderr.write('\nEnd: {} ({})\n'.format(t2,t2-t1))
        sys.stderr.flush()
    if ret != 0:
        sys.stderr.write('\nTerminated command in process {}.\n\n'.format(args.proc_name))
        sys.stderr.flush()
        raise ValueError('ERROR')
    return ret

def print_message(message,print_time=True):
    if message is not None:
        sys.stderr.write('\n'+message+'\n')
        sys.stderr.flush()
    if print_time:
        t1 = datetime.now()
        sys.stderr.write('\n{}\n'.format(t1))
        sys.stderr.flush()
    return

# Check files/folders
start_dtim = datetime.strptime(args.start_date,'%Y%m%d')
end_dtim = datetime.strptime(args.end_date,'%Y%m%d')
first_dtim = datetime.strptime(args.first_date,'%Y%m%d')
last_dtim = datetime.strptime(args.last_date,'%Y%m%d')
d1 = start_dtim+timedelta(days=60)
d2 = end_dtim+timedelta(days=240)
if not os.path.exists(args.s2_data):
    os.makedirs(args.s2_data)
if not os.path.isdir(args.s2_data):
    raise ValueError('{}: error, no such folder >>> {}'.format(args.proc_name,args.s2_data))
if not os.path.isdir(args.rcdir):
    raise ValueError('{}: error, no such folder >>> {}'.format(args.proc_name,args.rcdir))

for site in opts.sites:
    site_low = site.lower()

    # Download Sentinel-1 GRD
    if site_low in ['cihea','bojongsoang']:
        data_years = np.arange(first_dtim.year,last_dtim.year+1,1)
        for year in data_years:
            ystr = '{}'.format(year)
            # Make file list
            tmp_fnam = mktemp(suffix='.csv')
            command = args.python_path
            command += ' {}'.format(os.path.join(args.scr_dir,'file_station_query.py'))
            command += ' --server "{}"'.format(args.server)
            command += ' --port "{}"'.format(args.port)
            command += ' --rcdir "{}"'.format(args.rcdir)
            command += ' --srcdir "{}"'.format('{}/{}/{}/{}'.format(args.s1_path,site,'GRD',year))
            command += ' --out_csv "{}"'.format(tmp_fnam)
            command += ' --max_layer 0'
            try:
                run_command(command,message='<<< Make Sentinel-1 GRD list ({}) >>>'.format(ystr))
            except Exception:
                if os.path.exists(tmp_fnam):
                    os.remove(tmp_fnam)
                continue
            df = pd.read_csv(tmp_fnam,comment='#')
            if os.path.exists(tmp_fnam):
                os.remove(tmp_fnam)
            df.columns = df.columns.str.strip()
            inds = []
            for index,row in df.iterrows():
                #fileName,nLayer,fileSize,modifiedDate,folderName,md5Checksum
                src_fnam = row['fileName']
                m = re.search('^S1[AB]_IW_GRDH_1SDV_('+'\d'*8+')T\S+\.zip$',src_fnam)
                if not m:
                    continue
                dstr = m.group(1)
                d = datetime.strptime(dstr,'%Y%m%d')
                if d < first_dtim or d > last_dtim:
                    continue
                inds.append(index)
            if len(inds) < 1:
                print_message('No Sentinel-1 GRD data for download ({})'.format(ystr),print_time=False)
                continue
            tmp_fnam = mktemp(suffix='.csv')
            df.loc[inds].to_csv(tmp_fnam,index=False)
            # Download Data
            command = args.python_path
            command += ' {}'.format(os.path.join(args.scr_dir,'file_station_download_file.py'))
            command += ' --server "{}"'.format(args.server)
            command += ' --port "{}"'.format(args.port)
            command += ' --rcdir "{}"'.format(args.rcdir)
            command += ' --inp_list "{}"'.format(tmp_fnam)
            command += ' --dstdir "{}"'.format(os.path.join(args.s1_data,site,'GRD',ystr))
            command += ' --verbose'
            if args.overwrite:
                command += ' --overwrite'
            run_command(command,message='<<< Download Sentinel-1 GRD ({}) >>>'.format(ystr))
            if os.path.exists(tmp_fnam):
                os.remove(tmp_fnam)

    # Download Planting data
    if site_low in ['cihea','bojongsoang']:
        if site_low in ['bojongsoang']:
            enam_list = ['.json','.dbf','.prj','.shp','.shx']
        else:
            enam_list = ['.tif','.json']
        data_years = np.arange(d1.year,d2.year+1,1)
        for year in data_years:
            ystr = '{}'.format(year)
            # Make file list
            tmp_fnam = mktemp(suffix='.csv')
            command = args.python_path
            command += ' {}'.format(os.path.join(args.scr_dir,'file_station_query.py'))
            command += ' --server "{}"'.format(args.server)
            command += ' --port "{}"'.format(args.port)
            command += ' --rcdir "{}"'.format(args.rcdir)
            command += ' --srcdir "{}"'.format('{}/{}/{}/{}'.format(args.trans_path,site,'planting',year))
            command += ' --out_csv "{}"'.format(tmp_fnam)
            command += ' --max_layer 1'
            try:
                run_command(command,message='<<< Make Planting data list ({}) >>>'.format(ystr))
            except Exception:
                if os.path.exists(tmp_fnam):
                    os.remove(tmp_fnam)
                continue
            df = pd.read_csv(tmp_fnam,comment='#')
            if os.path.exists(tmp_fnam):
                os.remove(tmp_fnam)
            df.columns = df.columns.str.strip()
            inds = []
            for index,row in df.iterrows():
                #fileName,nLayer,fileSize,modifiedDate,folderName,md5Checksum
                items = row['fileName'].strip().split('/')
                if len(items) != 2:
                    continue
                src_dnam = items[0]
                src_fnam = items[1]
                m = re.search('_('+'\d'*8+')_final(\.\S+)$',src_fnam)
                if not m:
                    continue
                dstr = m.group(1)
                enam = m.group(2)
                if not enam in enam_list:
                    continue
                d = datetime.strptime(dstr,'%Y%m%d')
                if d < d1 or d > d2:
                    continue
                src_pnam = row['folderName'].strip()
                df.loc[index,'fileName'] = src_fnam
                df.loc[index,'folderName'] = '{}/{}'.format(src_pnam,src_dnam)
                df.loc[index,'nLayer'] = 0
                inds.append(index)
            if len(inds) < 1:
                print_message('No planting data for download ({})'.format(ystr),print_time=False)
                continue
            tmp_fnam = mktemp(suffix='.csv')
            df.loc[inds].to_csv(tmp_fnam,index=False)
            # Download Data
            command = args.python_path
            command += ' {}'.format(os.path.join(args.scr_dir,'file_station_download_file.py'))
            command += ' --server "{}"'.format(args.server)
            command += ' --port "{}"'.format(args.port)
            command += ' --rcdir "{}"'.format(args.rcdir)
            command += ' --inp_list "{}"'.format(tmp_fnam)
            command += ' --dstdir "{}"'.format(os.path.join(args.s1_data,site,'planting',ystr))
            command += ' --verbose'
            if args.overwrite:
                command += ' --overwrite'
            run_command(command,message='<<< Download planting data ({}) >>>'.format(ystr))
            if os.path.exists(tmp_fnam):
                os.remove(tmp_fnam)

    # Download Sentinel-2 L2A
    if site_low in ['bojongsoang']:
        data_years = np.arange(first_dtim.year,last_dtim.year+1,1)
        for year in data_years:
            ystr = '{}'.format(year)
            # Make file list
            tmp_fnam = mktemp(suffix='.csv')
            command = args.python_path
            command += ' {}'.format(os.path.join(args.scr_dir,'file_station_query.py'))
            command += ' --server "{}"'.format(args.server)
            command += ' --port "{}"'.format(args.port)
            command += ' --rcdir "{}"'.format(args.rcdir)
            command += ' --srcdir "{}"'.format('{}/{}/{}/{}'.format(args.s2_path,site,'L2A',year))
            command += ' --out_csv "{}"'.format(tmp_fnam)
            command += ' --max_layer 0'
            try:
                run_command(command,message='<<< Make Sentinel-2 L2A list ({}) >>>'.format(ystr))
            except Exception:
                if os.path.exists(tmp_fnam):
                    os.remove(tmp_fnam)
                continue
            df = pd.read_csv(tmp_fnam,comment='#')
            if os.path.exists(tmp_fnam):
                os.remove(tmp_fnam)
            df.columns = df.columns.str.strip()
            inds = []
            for index,row in df.iterrows():
                #fileName,nLayer,fileSize,modifiedDate,folderName,md5Checksum
                src_fnam = row['fileName']
                m = re.search('^S2[AB]_MSIL2A_('+'\d'*8+')T\S+\.zip$',src_fnam)
                if not m:
                    continue
                dstr = m.group(1)
                d = datetime.strptime(dstr,'%Y%m%d')
                if d < first_dtim or d > last_dtim:
                    continue
                inds.append(index)
            if len(inds) < 1:
                print_message('No Sentinel-2 L2A data for download ({})'.format(ystr),print_time=False)
                continue
            tmp_fnam = mktemp(suffix='.csv')
            df.loc[inds].to_csv(tmp_fnam,index=False)
            # Download Data
            command = args.python_path
            command += ' {}'.format(os.path.join(args.scr_dir,'file_station_download_file.py'))
            command += ' --server "{}"'.format(args.server)
            command += ' --port "{}"'.format(args.port)
            command += ' --rcdir "{}"'.format(args.rcdir)
            command += ' --inp_list "{}"'.format(tmp_fnam)
            command += ' --dstdir "{}"'.format(os.path.join(args.s2_data,site,'L2A',ystr))
            command += ' --verbose'
            if args.overwrite:
                command += ' --overwrite'
            run_command(command,message='<<< Download Sentinel-2 L2A ({}) >>>'.format(ystr))
            if os.path.exists(tmp_fnam):
                os.remove(tmp_fnam)

    # Download Sentinel-2 geocor/indices/parcel/atcor/interp/tentative_interp
    for i,(targ,enam) in enumerate(zip(['geocor','indices','parcel','atcor','interp','tentative_interp'],
                                       ['geocor\S*\.\S*','indices\S*\.\S*','parcel\S*\.\S*','\S*\.\S*','interp\S*\.\S*','interp\S*\.\S*'])):
        if site_low in ['testsite'] and targ in ['geocor','indices']:
            continue
        data_years = np.arange(first_dtim.year,last_dtim.year+1,1)
        for year in data_years:
            ystr = '{}'.format(year)
            # Make file list
            tmp_fnam = mktemp(suffix='.csv')
            command = args.python_path
            command += ' {}'.format(os.path.join(args.scr_dir,'file_station_query.py'))
            command += ' --server "{}"'.format(args.server)
            command += ' --port "{}"'.format(args.port)
            command += ' --rcdir "{}"'.format(args.rcdir)
            command += ' --srcdir "{}"'.format('{}/{}/{}/{}'.format(args.s2_path,site,targ,year))
            command += ' --out_csv "{}"'.format(tmp_fnam)
            command += ' --max_layer 0'
            try:
                run_command(command,message='<<< Make Sentinel-2 {} list ({}) >>>'.format(targ,ystr))
            except Exception:
                if os.path.exists(tmp_fnam):
                    os.remove(tmp_fnam)
                continue
            df = pd.read_csv(tmp_fnam,comment='#')
            if os.path.exists(tmp_fnam):
                os.remove(tmp_fnam)
            df.columns = df.columns.str.strip()
            inds = []
            for index,row in df.iterrows():
                #fileName,nLayer,fileSize,modifiedDate,folderName,md5Checksum
                src_fnam = row['fileName']
                m = re.search('^('+'\d'*8+')_'+enam,src_fnam)
                if not m:
                    continue
                dstr = m.group(1)
                d = datetime.strptime(dstr,'%Y%m%d')
                if d < first_dtim or d > last_dtim:
                    continue
                inds.append(index)
            if len(inds) < 1:
                print_message('No Sentinel-2 {} data for download ({})'.format(targ,ystr),print_time=False)
                continue
            tmp_fnam = mktemp(suffix='.csv')
            df.loc[inds].to_csv(tmp_fnam,index=False)
            # Download Data
            command = args.python_path
            command += ' {}'.format(os.path.join(args.scr_dir,'file_station_download_file.py'))
            command += ' --server "{}"'.format(args.server)
            command += ' --port "{}"'.format(args.port)
            command += ' --rcdir "{}"'.format(args.rcdir)
            command += ' --inp_list "{}"'.format(tmp_fnam)
            command += ' --dstdir "{}"'.format(os.path.join(args.s2_data,site,targ,ystr))
            command += ' --verbose'
            if args.overwrite:
                command += ' --overwrite'
            run_command(command,message='<<< Download Sentinel-2 {} ({}) >>>'.format(targ,ystr))
            if os.path.exists(tmp_fnam):
                os.remove(tmp_fnam)
