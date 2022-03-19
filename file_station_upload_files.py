#!/usr/bin/env python
import os
import sys
import re
import shutil
import hashlib
from base64 import b64encode
import requests
from requests.packages.urllib3.exceptions import InsecureRequestWarning
import logging
from http.client import HTTPConnection
import time
from datetime import datetime,timedelta
import pytz
from glob import glob
from pathlib import Path
from optparse import OptionParser,IndentedHelpFormatter

# Constants
B = 1
KB = 1024
MB = KB*1024
GB = MB*1024
TB = GB*1024

# Default values
HOME = os.environ.get('HOME')
if HOME is None:
    HOME = os.environ.get('HOMEPATH')
RCDIR = HOME
PORT = 443
MAX_ITEM = 10000
CHUNK_SIZE = GB

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.set_usage('Usage: %prog list_of_input_file [options]')
parser.add_option('-D','--dstdir',default=None,help='Destination directory (%default)')
parser.add_option('--rcdir',default=RCDIR,help='Directory where .netrc exists (%default)')
parser.add_option('-N','--server',default=None,help='Name of the server (%default)')
parser.add_option('-P','--port',default=PORT,type='int',help='Port# of the server (%default)')
parser.add_option('-M','--max_item',default=MAX_ITEM,type='int',help='Max# of items for listing (%default)')
parser.add_option('-C','--chunk_size',default=CHUNK_SIZE,type='int',help='Chunk size in byte (%default)')
parser.add_option('-l','--logging',default=False,action='store_true',help='Logging mode (%default)')
parser.add_option('-v','--verbose',default=False,action='store_true',help='Verbose mode (%default)')
parser.add_option('-d','--debug',default=False,action='store_true',help='Debug mode (%default)')
parser.add_option('--overwrite',default=False,action='store_true',help='Overwrite mode (%default)')
(opts,args) = parser.parse_args()
if len(args) < 1:
    parser.print_help()
    sys.exit(0)
fnams = args
if opts.server is None or opts.dstdir is None:
    raise ValueError('Error, server={}, dstdir={}'.format(opts.server,opts.dstdir))

def get_size(fnam):
    s = 'None'
    if os.path.exists(fnam):
        size = os.path.getsize(fnam)
        if size > TB:
            s = '{:.2f} TB'.format(size/TB)
        elif size > GB:
            s = '{:.2f} GB'.format(size/GB)
        elif size > MB:
            s = '{:.2f} MB'.format(size/MB)
        elif size > KB:
            s = '{:.2f} KB'.format(size/KB)
        else:
            s = '{} B'.format(size)
    return s

def get_time(s):
    try:
        t = datetime.utcfromtimestamp(s).replace(tzinfo=pytz.UTC)
    except Exception:
        raise ValueError('Error in time >>> '+s)
    return t

def list_file(path=None):
    ds = []
    fs = {}
    try:
        resp = session.get(common_url+'&func=get_list&path={}&limit={}'.format(path,opts.max_item))
        params = resp.json()
        if 'datas' in params:
            items = params['datas']
            for item in items:
                if item['isfolder'] == 1:
                    ds.append(item['filename'])
                else:
                    fs.update({item['filename']:{'size':int(item['filesize']),'mtime':get_time(item['epochmt'])}})
        else:
            status = params['status']
            raise ValueError('Error, status={}'.format(status))
    except Exception as e:
        sys.stderr.write(str(e)+'\n')
        sys.stderr.write('Error in listing file >>> {}\n'.format(path))
        sys.stderr.flush()
        return None,None
    return ds,fs

def query_file(path):
    parent = os.path.dirname(path)
    target = os.path.basename(path)
    try:
        resp = session.get(common_url+'&func=checksum&file_total=1&path={}&file_name={}'.format(parent,target))
        params = resp.json()
        if 'datas' in params:
            item = params['datas'][0]
        else:
            status = params['status']
            raise ValueError('Error, status={}'.format(status))
    except Exception as e:
        sys.stderr.write(str(e)+'\n')
        sys.stderr.write('Error in querying file >>> {}\n'.format(path))
        sys.stderr.flush()
        return None
    return item['filename'],int(item['filesize']),get_time(item['mt']),item['checksum']

def delete_file(path):
    parent = os.path.dirname(path)
    target = os.path.basename(path)
    try:
        resp = session.get(common_url+'&func=delete&file_total=1&path={}&file_name={}'.format(parent,target))
        status = resp.json()['status']
        if status != 1:
            raise ValueError('Error, status={}'.format(status))
    except Exception as e:
        sys.stderr.write(str(e)+'\n')
        sys.stderr.write('Error in deleting folder >>> {}\n'.format(path))
        sys.stderr.flush()
        return -1
    return 0

folders = []

def make_folder(path):
    global folders
    if path in folders:
        return 0
    parent = os.path.dirname(path)
    target = os.path.basename(path)
    flag = False
    try:
        resp = session.get(common_url+'&func=stat&file_total=1&path={}&file_name={}'.format(parent,target))
        params = resp.json()
        if 'datas' in params:
            flag = True
            item = params['datas'][0]
            if item['isfolder'] == 1:
                folders.append(path)
                return 0
            elif item['owner'] != '':
                raise ValueError('Error, file exists >>> {}'.format(path))
    except Exception as e:
        if flag:
            sys.stderr.write(str(e)+'\n')
            sys.stderr.write('Error in making folder >>> {}\n'.format(path))
            sys.stderr.flush()
            return -1
    try:
        resp = session.get(common_url+'&func=createdir&dest_path={}&dest_folder={}'.format(parent,target))
        status = resp.json()['status']
        if status != 1:
            raise ValueError('Error, status={}'.format(status))
    except Exception as e:
        sys.stderr.write(str(e)+'\n')
        sys.stderr.write('Error in making folder >>> {}\n'.format(path))
        sys.stderr.flush()
        return -1
    folders.append(path)
    return 0

def make_folders(path):
    normalized_path = os.path.normpath(path)
    path_components = normalized_path.split(os.sep)
    if path[0] == os.sep:
        dnam = os.sep
    else:
        dnam = ''
    for p in path_components:
        if not p:
            continue
        dnam = os.path.join(dnam,p)
        make_folder(dnam)

def read_in_chunks(file_object,chunk_size=GB):
    while True:
        data = file_object.read(chunk_size)
        if not data:
            break
        yield data

def upload_file(fnam,gnam,chunk_size=GB):
    parent = os.path.dirname(gnam)
    target = os.path.basename(gnam)
    if parent == '':
        raise ValueError('Error in parent >>> {}'.format(gnam))
    elif not parent in folders:
        raise IOError('Error, no such folder >>> '+parent)
    ds,fs = list_file(parent)
    if target in fs:
        if opts.overwrite:
            if opts.verbose:
                sys.stderr.write('File exists, delete >>> '+target+'\n')
                sys.stderr.flush()
            if delete_file(gnam) < 0:
                raise IOError('Error in deleting file >>> '+target)
        else:
            if opts.verbose:
                sys.stderr.write('File exists, skip >>> '+target+'\n')
                sys.stderr.flush()
            return query_file(gnam)
    # Upload file
    byte_size = os.path.getsize(fnam)
    ftim = int(os.path.getmtime(fnam)+0.5)
    if opts.verbose:
        tstr = datetime.now()
        sys.stderr.write('{:%Y-%m-%dT%H:%M:%S} Uploading file ({}) >>> {}\n'.format(tstr,get_size(fnam),fnam))
        sys.stderr.flush()
    flag = False
    try:
        resp = session.get(common_url+'&func=start_chunked_upload&upload_root_dir={}'.format(parent))
        status = resp.json()['status']
        if status != 0:
            raise ValueError('Error, status={}'.format(status))
        upload_id = resp.json()['upload_id']
        with open(fnam,'rb') as fp:
            t1 = datetime.now()
            offset = 0
            for chunk in read_in_chunks(fp,chunk_size):
                data_size = len(chunk)
                url = common_url+'&func=chunked_upload&upload_id={}&upload_root_dir={}&dest_path={}&upload_name={}&filesize={}&offset={}&overwrite={}&settime=1&mtime={}'.format(upload_id,parent,parent,target,byte_size,offset,(1 if opts.overwrite else 0),ftim)
                offset += data_size
                if offset < byte_size:
                    url += '&multipart=1'
                resp = session.post(url,files=(('fileName',(None,target)),('file',('blob',chunk,'application/octet-stream'))))
                status = resp.json()['status']
                if status != 1:
                    raise ValueError('Error, status={}'.format(status))
                if offset < byte_size:
                    t2 = datetime.now()
                    dt = (t2-t1).total_seconds()
                    rate = data_size/dt
                    t3 = t2+timedelta(seconds=(byte_size-offset)/rate)
                    sys.stderr.write('{:%Y-%m-%dT%H:%M:%S} Upload {:6.2f} % @ {:8.3f} Mbps, Expected completion at {:%Y-%m-%dT%H:%M:%S}\n'.format(t2,100.0*offset/byte_size,rate*8.0e-6,t3))
                    sys.stderr.flush()
                    t1 = t2
        flag = True
    except Exception as e:
        sys.stderr.write(str(e)+'\n')
        sys.stderr.write('Error in uploading file >>> {}\n'.format(fnam))
        sys.stderr.flush()
    if opts.verbose:
        tend = datetime.now()
        dt = (tend-tstr).total_seconds()
        rate = byte_size/dt
        if flag:
            sys.stderr.write('{:%Y-%m-%dT%H:%M:%S} Upload completed in {:.2f} seconds ({:8.3f} Mbps) >>> {}\n'.format(tend,dt,rate*8.0e-6,fnam))
        else:
            sys.stderr.write('{:%Y-%m-%dT%H:%M:%S} Upload stopped in {:.2f} seconds >>> {}\n'.format(tend,dt,fnam))
        sys.stderr.flush()
    # Check uploaded file
    result = query_file(gnam)
    if result is None:
        sys.stderr.write('Warning, faild in uploading file >>> {}\n'.format(gnam))
        sys.stderr.flush()
        return '',-1,'',''
    else:
        return result

def upload_and_check_file(fnam,gnam,chunk_size=GB):
    title = os.path.basename(gnam)
    size = os.path.getsize(fnam)
    title_dst,size_dst,mdate_dst,md5_dst = upload_file(fnam,gnam,chunk_size)
    if (title_dst.lower() == title.lower()) and (size_dst == size):
        with open(fnam,'rb') as fp:
            md5 = hashlib.md5(fp.read()).hexdigest()
        if (md5_dst.lower() == md5.lower()):
            if opts.verbose:
                sys.stderr.write('Successfully uploaded >>> {}\n'.format(fnam))
                sys.stderr.flush()
            return 0
        else:
            if opts.verbose:
                sys.stderr.write('Warning, title={} ({}), size={} ({}), md5={} ({})\n'.format(title_dst,title,size_dst,size,md5_dst,md5))
                sys.stderr.flush()
            return -1
    else:
        if opts.verbose:
            sys.stderr.write('Warning, title={} ({}), size={} ({})\n'.format(title_dst,title,size_dst,size))
            sys.stderr.flush()
        return -1

fnam = os.path.join(opts.rcdir,'.netrc')
if not os.path.exists(fnam):
    raise IOError('Error, no such file >>> '+fnam)
server = None
username = None
password = None
flag = False
with open(fnam,'r') as fp:
    for line in fp:
        m = re.search('machine\s+(\S+)',line)
        if m:
            if re.search(opts.server,m.group(1)):
                flag = True
                server = m.group(1)
            else:
                flag = False
            continue
        m = re.search('login\s+(\S+)',line)
        if m:
            if flag:
                username = m.group(1)
            continue
        m = re.search('password\s+(\S+)',line)
        if m:
            if flag:
                password = m.group(1)
            continue
if server is None or username is None or password is None:
    raise ValueError('Error, server={}, username={}, password={}'.format(server,username,password))
encode_string = b64encode(password.encode()).decode()

# Create a requests session
requests.packages.urllib3.disable_warnings(InsecureRequestWarning)
session = requests.sessions.Session()
session.verify = False
# Get session ID
url = 'https://{}:{}/cgi-bin/authLogin.cgi?user={}&pwd={}'.format(server,opts.port,username,encode_string)
try:
    resp = session.get(url)
    m = re.search('<authSid><!\[CDATA\[(\S+)\]\]><\/authSid>',resp.text)
    sid = m.group(1)
except Exception as e:
    sys.stderr.write('Login failed >>> '+str(e)+'\n')
    sys.stderr.flush()
    sys.exit()
common_url = 'https://{}:{}/cgi-bin/filemanager/utilRequest.cgi?sid={}'.format(server,opts.port,sid)
if opts.logging:
    log = logging.getLogger('urllib3')
    log.setLevel(logging.DEBUG)
    stream = logging.StreamHandler()
    stream.setLevel(logging.DEBUG)
    log.addHandler(stream)
    HTTPConnection.debuglevel = 1
# Upload file
make_folders(opts.dstdir)
for fnam in fnams:
    f = os.path.basename(fnam)
    gnam = os.path.join(opts.dstdir,f)
    upload_and_check_file(fnam,gnam,opts.chunk_size)
