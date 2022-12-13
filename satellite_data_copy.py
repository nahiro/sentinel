import os
import sys
import re
from datetime import datetime,timedelta
from glob import glob
import numpy as np
import pandas as pd

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
start_dtim = datetime.strptime(args.start_date,args.date_fmt)
end_dtim = datetime.strptime(args.end_date,args.date_fmt)
first_dtim = datetime.strptime(args.first_date,args.date_fmt)
last_dtim = datetime.strptime(args.last_date,args.date_fmt)
d1 = start_dtim+timedelta(days=60)
d2 = end_dtim+timedelta(days=240)
if not os.path.exists(args.s2_data):
    os.makedirs(args.s2_data)
if not os.path.isdir(args.s2_data):
    raise ValueError('{}: error, no such folder >>> {}'.format(args.proc_name,args.s2_data))
netrc_dir = os.path.dirname(args.netrc_fnam)
if not os.path.isdir(netrc_dir):
    raise ValueError('{}: error, no such folder >>> {}'.format(args.proc_name,netrc_dir))

# Download Planting data
itarg = self.list_labels['dflag'].index('planting')
iflag = self.list_labels['oflag'].index('planting')
if args.dflag[itarg]:
    if 'bojongsoang' in args.trans_path.lower():
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
        command += ' --rcdir "{}"'.format(netrc_dir)
        command += ' --srcdir "{}"'.format('{}/{}'.format(args.trans_path,year))
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
        command += ' --rcdir "{}"'.format(netrc_dir)
        command += ' --inp_list "{}"'.format(tmp_fnam)
        command += ' --dstdir "{}"'.format(os.path.join(args.s1_data,'planting',ystr))
        command += ' --verbose'
        if args.oflag[iflag]:
            command += ' --overwrite'
        run_command(command,message='<<< Download planting data ({}) >>>'.format(ystr))
        if os.path.exists(tmp_fnam):
            os.remove(tmp_fnam)

# Download Sentinel-2 L2A
itarg = self.list_labels['dflag'].index('L2A')
iflag = self.list_labels['oflag'].index('L2A')
if args.dflag[itarg]:
    keys = []
    for key in [s.strip() for s in args.search_key.split(',')]:
        if key:
            keys.append(key)
    if len(keys) < 1:
        keys = None
    data_years = np.arange(first_dtim.year,last_dtim.year+1,1)
    for year in data_years:
        ystr = '{}'.format(year)
        # Make file list
        tmp_fnam = mktemp(suffix='.csv')
        command = args.python_path
        command += ' {}'.format(os.path.join(args.scr_dir,'file_station_query.py'))
        command += ' --server "{}"'.format(args.server)
        command += ' --port "{}"'.format(args.port)
        command += ' --rcdir "{}"'.format(netrc_dir)
        command += ' --srcdir "{}"'.format('{}/{}'.format(args.l2a_path,year))
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
            if keys is not None:
                flag = False
                for key in keys:
                    if not key in src_fnam:
                        flag = True
                        break
                if flag:
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
        command += ' --rcdir "{}"'.format(netrc_dir)
        command += ' --inp_list "{}"'.format(tmp_fnam)
        command += ' --dstdir "{}"'.format(os.path.join(args.s2_data,'L2A',ystr))
        command += ' --verbose'
        if args.oflag[iflag]:
            command += ' --overwrite'
        run_command(command,message='<<< Download Sentinel-2 L2A ({}) >>>'.format(ystr))
        if os.path.exists(tmp_fnam):
            os.remove(tmp_fnam)

# Download Sentinel-2 geocor/indices/parcel/atcor/interp/tentative_interp
for i,(targ,pnam,enam) in enumerate(zip(['geocor','indices','parcel','atcor','interp','tentative_interp'],
                                        ['geocor','indices','parcel','atcor','interp','tentative'],
                                        ['geocor\S*\.\S*','indices\S*\.\S*','parcel\S*\.\S*','\S*\.\S*','interp\S*\.\S*','interp\S*\.\S*'])):
    itarg = self.list_labels['dflag'].index(targ)
    iflag = self.list_labels['oflag'].index(targ)
    if args.dflag[itarg]:
        data_years = np.arange(first_dtim.year,last_dtim.year+1,1)
        for year in data_years:
            ystr = '{}'.format(year)
            # Make file list
            tmp_fnam = mktemp(suffix='.csv')
            command = args.python_path
            command += ' {}'.format(os.path.join(args.scr_dir,'file_station_query.py'))
            command += ' --server "{}"'.format(args.server)
            command += ' --port "{}"'.format(args.port)
            command += ' --rcdir "{}"'.format(netrc_dir)
            command += ' --srcdir "{}"'.format('{}/{}'.format(getattr(args,'{}_path'.format(pnam)),year))
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
            command += ' --rcdir "{}"'.format(netrc_dir)
            command += ' --inp_list "{}"'.format(tmp_fnam)
            command += ' --dstdir "{}"'.format(os.path.join(args.s2_data,targ,ystr))
            command += ' --verbose'
            if args.oflag[iflag]:
                command += ' --overwrite'
            run_command(command,message='<<< Download Sentinel-2 {} ({}) >>>'.format(targ,ystr))
            if os.path.exists(tmp_fnam):
                os.remove(tmp_fnam)
