#!/usr/bin/env python
import os
import sys
import re
import shutil
import hashlib
import requests
from requests.packages.urllib3.exceptions import InsecureRequestWarning
import atexit
import time
from datetime import datetime,timedelta
from glob import glob
from pathlib import Path
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.keys import Keys
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
DRVDIR = os.path.join(HOME,'local','bin')

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('-S','--srcdir',default=None,help='Source directory (%default)')
parser.add_option('-s','--subdir',default=None,action='append',help='Sub directory (%default)')
parser.add_option('-D','--dstdir',default=None,help='Destination directory (%default)')
parser.add_option('-L','--locdir',default=None,help='Local destination directory (%default)')
parser.add_option('--drvdir',default=DRVDIR,help='Directory where chromedriver exists (%default)')
parser.add_option('--rcdir',default=RCDIR,help='Directory where .netrc exists (%default)')
parser.add_option('-K','--keep_folder',default=None,action='append',help='Directory to keep (%default)')
parser.add_option('-I','--ignore_file',default=None,action='append',help='File to ignore (%default)')
parser.add_option('-p','--port',default=None,type='int',help='Port# of Chrome to use (%default)')
parser.add_option('--skip_login',default=False,action='store_true',help='Skip login procedure (%default)')
parser.add_option('-H','--headless',default=False,action='store_true',help='Headless mode (%default)')
parser.add_option('-v','--verbose',default=False,action='store_true',help='Verbose mode (%default)')
parser.add_option('-d','--debug',default=False,action='store_true',help='Debug mode (%default)')
parser.add_option('--overwrite',default=False,action='store_true',help='Overwrite mode (%default)')
(opts,args) = parser.parse_args()
if opts.srcdir is None or opts.subdir is None or opts.dstdir is None or opts.locdir is None:
    raise ValueError('Error, srcdir={}, subdir={}, dstdir={}, locdir={}'.format(opts.srcdir,opts.subdir,opts.dstdir,opts.locdir))
keep_folder = []
keep_folder_lower = []
if opts.keep_folder is not None:
    for f in opts.keep_folder:
        for p in glob(os.path.normpath(os.path.join(opts.srcdir,f))):
            if os.path.isdir(p):
                keep_folder.append(p)
                keep_folder_lower.append(p.lower())
ignore_file = []
ignore_file_lower = []
if opts.ignore_file is not None:
    for f in opts.ignore_file:
        for p in glob(os.path.normpath(os.path.join(opts.srcdir,f))):
            if not os.path.isdir(p):
                ignore_file.append(p)
                ignore_file_lower.append(p.lower())
if opts.verbose:
    if len(keep_folder) > 0:
        sys.stderr.write('keep_folder:\n')
        for f in keep_folder:
            sys.stderr.write(f+'\n')
        sys.stderr.flush()
    if len(ignore_file) > 0:
        sys.stderr.write('ignore_file:\n')
        for f in ignore_file:
            sys.stderr.write(f+'\n')
        sys.stderr.flush()

opts.srcdir = os.path.abspath(opts.srcdir)
topdir = os.getcwd()
#os.chdir(opts.drvdir)
def clean_up():
    sys.stderr.write('clean_up called.\n\n')
    sys.stderr.flush()
    os.chdir(topdir)
atexit.register(clean_up)

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
    # '2021-11-26T09:24:23.01891815Z'
    m = re.search('(\d\d\d\d-\d\d-\d\dT\d\d:\d\d):(\d\d\.\d+)Z',s)
    if m:
        t = datetime.strptime(m.group(1)+'Z','%Y-%m-%dT%H:%M%z')+timedelta(seconds=float(m.group(2)))
    else:
        try:
            t = datetime.strptime(s,'%Y-%m-%dT%H:%M:%S%z')
        except Exception:
            raise ValueError('Error in time >>> '+s)
    return t

def get_url(path,root='files'):
    if path is None or path == '.' or len(path) < 1: # current folder
        m = re.search('https://'+server+'/files(\S+)',driver.current_url)
        if m:
            path = m.group(1)
        else:
            path = '/'
    if path[0] != '/': # relative path
        m = re.search('https://'+server+'/files(\S+)',driver.current_url)
        if m:
            parent = m.group(1)
        else:
            parent = '/'
        if parent[-1] != '/':
            parent += '/'
        path = parent+path
    if root == 'files':
        url = 'https://{}/files{}'.format(server,path)
    elif root == 'resources':
        url = 'https://{}/api/resources{}'.format(server,path)
    else:
        raise ValueError('Error, in root >>> '+root)
    return url

def change_folder(path):
    url = get_url(path)
    if driver.current_url != url:
        driver.get(url)
    time.sleep(1)
    if re.search('This location can\'t be reached',driver.page_source):
        return -1
    return 0

def list_file(path=None):
    ds = []
    fs = {}
    url = get_url(path,root='resources')
    try:
        resp = session.get(url,verify=False)
        params = resp.json()
        items = params['items']
        for item in items:
            if item['isDir']:
                ds.append(item['name'])
            else:
                fs.update({item['name']:{'size':item['size'],'mtime':get_time(item['modified'])}})
    except Exception:
        return None,None
    return ds,fs

def query_file(path):
    url = get_url(path,root='resources')+'?checksum=md5'
    try:
        resp = session.get(url,verify=False)
        item = resp.json()
    except Exception:
        return None
    return item['name'],item['size'],get_time(item['modified']),item['checksums']['md5']

def delete_file(path):
    url = get_url(path,root='resources')
    try:
        resp = session.delete(url,verify=False)
    except Exception:
        return -1
    return 0

folders = []

def make_folder(path):
    global folders
    if path in folders:
        return 0
    parent = os.path.dirname(path)
    target = os.path.basename(path)
    ds,fs = list_file(parent)
    if ds is None or fs is None:
        raise IOError('Error, no such folder >>> '+parent)
    if target in ds:
        folders.append(path)
        return 0
    url = get_url(path,root='resources')
    if url[-1] != '/':
        url += '/'
    try:
        resp = session.post(url,verify=False)
    except Exception:
        return -1
    folders.append(path)
    return 0

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
    if opts.verbose:
        tstr = datetime.now()
        sys.stderr.write('{:%Y-%m-%dT%H:%M:%S} Uploading file ({}) >>> {}\n'.format(tstr,get_size(fnam),fnam))
        sys.stderr.flush()
    url = get_url(gnam,root='resources')
    try:
        with open(fnam,'rb') as fp:
            session.post(url,data=read_in_chunks(fp,chunk_size),verify=False)
    except Exception as e:
        sys.stderr.write(str(e)+'\n')
        sys.stderr.flush()
    if opts.verbose:
        tend = datetime.now()
        dt = (tend-tstr).total_seconds()
        bit_size = os.path.getsize(fnam)*8
        rate = bit_size/dt
        sys.stderr.write('{:%Y-%m-%dT%H:%M:%S} Upload completed in {:.2f} seconds ({:8.3f} Mbps) >>> {}\n'.format(tend,dt,rate*1.0e-6,fnam))
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
            if re.search('satreps',m.group(1)):
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

if opts.port is not None or opts.headless:
    options = Options()
    if opts.port:
        options.add_experimental_option('debuggerAddress','localhost:{}'.format(opts.port))
    if opts.headless:
        options.add_argument('--headless')
        options.add_argument('--disable-gpu')
        options.add_argument('--window-size=1920,1080')
    driver = webdriver.Chrome(os.path.join(opts.drvdir,'chromedriver'),options=options)
else:
    driver = webdriver.Chrome(os.path.join(opts.drvdir,'chromedriver'))

# Login to the NAS server
if not opts.skip_login:
    driver.get('https://{}/login'.format(server))
    time.sleep(1)
    inputs = driver.find_elements_by_class_name('input')
    if len(inputs) != 2:
        raise ValueError('Error, len(inputs)={}'.format(len(inputs)))
    inputs[0].send_keys((Keys.CONTROL+'a'))
    inputs[0].send_keys(Keys.DELETE)
    inputs[0].send_keys(username)
    time.sleep(1)
    inputs[1].send_keys((Keys.CONTROL+'a'))
    inputs[1].send_keys(Keys.DELETE)
    inputs[1].send_keys(password)
    time.sleep(1)
    buttons = driver.find_elements_by_class_name('button')
    if len(buttons) != 1:
        raise ValueError('Error, len(buttons)={}'.format(len(buttons)))
    buttons[0].click()
    time.sleep(1)
# Create a requests session
requests.packages.urllib3.disable_warnings(InsecureRequestWarning)
session = requests.sessions.Session()
for cookie in driver.get_cookies():
    c = {cookie['name']:cookie['value']}
    session.cookies.update(c)
# Upload file
for subdir in opts.subdir:
    make_folder(os.path.join(opts.dstdir,subdir))
    for root,ds,fs in os.walk(os.path.join(opts.srcdir,subdir)):
        curdir = os.path.relpath(root,opts.srcdir)
        if opts.verbose:
            sys.stderr.write('#####################\n')
            sys.stderr.write(curdir+'\n')
            if len(ds) > 0:
                sys.stderr.write('Folders: {}\n'.format(ds))
            if len(fs) > 0:
                sys.stderr.write('Files  : {}\n'.format(fs))
            sys.stderr.flush()
        if curdir == os.curdir:
            srcdir = opts.srcdir
            dstdir = opts.dstdir
            locdir = opts.locdir
        else:
            srcdir = os.path.join(opts.srcdir,curdir)
            dstdir = os.path.join(opts.dstdir,curdir)
            locdir = os.path.join(opts.locdir,curdir)
        #print(srcdir,'-----',dstdir)
        if not dstdir in folders:
            if make_folder(dstdir) != 0:
                raise IOError('Error, faild in making folder >>> '+dstdir)
        if not os.path.exists(locdir):
            os.makedirs(locdir)
        if not os.path.isdir(locdir):
            raise IOError('Error, no such folder >>> '+locdir)
        for f in sorted(fs):
            fnam = os.path.join(srcdir,f)
            gnam = os.path.join(dstdir,f)
            lnam = os.path.join(locdir,f)
            if fnam.lower() in ignore_file_lower:
                continue
            if upload_and_check_file(fnam,gnam) == 0:
                shutil.move(fnam,lnam)
                if opts.debug and not os.path.exists(fnam) and os.path.exists(lnam):
                    sys.stderr.write('Moved {} to {}\n'.format(fnam,lnam))
                    sys.stderr.flush()
        if len(os.listdir(srcdir)) == 0:
            if not srcdir.lower() in keep_folder_lower:
                if os.path.isdir(locdir):
                    os.rmdir(srcdir)
                    if opts.debug and not os.path.exists(srcdir):
                        sys.stderr.write('Removed {}\n'.format(srcdir))
                        sys.stderr.flush()
                else:
                    sys.stderr.write('Warning, no such directory >>> {}\n'.format(locdir))
                    sys.stderr.flush()
# Remove empty directories
for subdir in opts.subdir:
    for root,ds,fs in os.walk(os.path.join(opts.srcdir,subdir),topdown=False):
        curdir = os.path.relpath(root,opts.srcdir)
        if curdir == os.curdir:
            srcdir = opts.srcdir
            locdir = opts.locdir
        else:
            srcdir = os.path.join(opts.srcdir,curdir)
            locdir = os.path.join(opts.locdir,curdir)
        if len(os.listdir(srcdir)) == 0:
            if not srcdir.lower() in keep_folder_lower:
                if os.path.isdir(locdir):
                    os.rmdir(srcdir)
                    if opts.debug and not os.path.exists(srcdir):
                        sys.stderr.write('Removed {}\n'.format(srcdir))
                        sys.stderr.flush()
                else:
                    sys.stderr.write('Warning, no such directory >>> {}\n'.format(locdir))
                    sys.stderr.flush()
