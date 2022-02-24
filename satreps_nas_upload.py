#!/usr/bin/env python
import os
import sys
import re
import shutil
import hashlib
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

def get_size(s):
    m = re.search('(\S+)\s+(\S+)',s)
    if m:
        size_unit = m.group(2).upper()
        if size_unit == 'B':
            size = int(m.group(1))
            size_error = 0
        elif size_unit == 'KB':
            size = float(m.group(1))*KB
            size_error = KB//10
        elif size_unit == 'MB':
            size = float(m.group(1))*MB
            size_error = MB//10
        elif size_unit == 'GB':
            size = float(m.group(1))*GB
            size_error = GB//10
        elif size_unit == 'TB':
            size = float(m.group(1))*TB
            size_error = TB//10
        else:
            raise ValueError('Error in size unit >>> '+s)
    else:
        raise ValueError('Error in size >>> '+s)
    return size,size_error

def get_time(s):
    # '2021-11-26T09:24:23.01891815Z'
    m = re.search('(\d\d\d\d-\d\d-\d\dT\d\d:\d\d):(\d\d\.\d+)Z',s)
    if m:
        t = datetime.strptime(m.group(1)+'Z','%Y-%m-%dT%H:%M%z')+timedelta(seconds=float(m.group(2)))
    else:
        raise ValueError('Error in time >>> '+s)
    return t

def list_file(path=None):
    ds = {}
    fs = {}
    if path is None or path == '.':
        pass
    else:
        if len(path) < 1:
            url = 'https://{}/files/'.format(server)
        elif path[0] == '/':
            url = 'https://{}/files{}'.format(server,path)
        else:
            url = 'https://{}/files/{}'.format(server,path)
        if driver.current_url != url:
            driver.get(url)
        time.sleep(1)
        if re.search('This location can\'t be reached',driver.page_source):
            return None,None
    items = driver.find_elements_by_class_name('item')
    for item in items:
        lines = item.text.splitlines()
        nline = len(lines)
        if nline == 4:
            if lines[0] == 'folder':
                ds.update({lines[1]:item})
            elif lines[0] == 'insert_drive_file':
                size,size_error = get_size(lines[2])
                t = get_time(item.find_element_by_tag_name('time').get_attribute('datetime'))
                fs.update({lines[1]:{'element':item,'size':size,'size_error':size_error,'mtime':t}})
            else:
                raise ValueError('Error, lines[0]={}'.format(lines[0]))
        elif nline == 3 and lines[1] != 'Size':
            size,size_error = get_size(lines[1])
            t = get_time(item.find_element_by_tag_name('time').get_attribute('datetime'))
            fs.update({lines[0]:{'element':item,'size':size,'size_error':size_error,'mtime':t}})
    return ds,fs

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
    else:
        action_new_folder = None
        actions = driver.find_elements_by_class_name('action')
        for action in actions:
            if action.get_attribute('title') == 'New folder':
                action_new_folder = action
                break
        if action_new_folder is None:
            raise ValueError('Error in finding New folder.')
        action_new_folder.click()
        time.sleep(1)
        inputs = driver.find_elements_by_class_name('input')
        if len(inputs) != 1:
            raise ValueError('Error, len(inputs)={}'.format(len(inputs)))
        inputs[0].send_keys(target)
        buttons = driver.find_elements_by_class_name('button')
        if len(buttons) != 2:
            raise ValueError('Error, len(buttons)={}'.format(len(buttons)))
        if buttons[1].text != 'CREATE':
            raise ValueError('Error, buttons[1].text={}'.format(buttons[1].text))
        buttons[1].click()
        time.sleep(1)
        folders.append(path)
        return 0
    return -1

def upload_file(fnam,gnam):
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
            fs[target]['element'].click()
            time.sleep(1)
            action_delete = None
            actions = driver.find_elements_by_class_name('action')
            for action in actions:
                if action.get_attribute('title') == 'Delete':
                    action_delete = action
                    break
            if action_delete is None:
                raise ValueError('Error in finding Delete.')
            action_delete.click()
            time.sleep(1)
            buttons = driver.find_elements_by_class_name('button')
            if len(buttons) != 2:
                raise ValueError('Error, len(buttons)={}'.format(len(buttons)))
            if buttons[1].text != 'DELETE':
                raise ValueError('Error, buttons[1].text={}'.format(buttons[1].text))
            buttons[1].click()
            time.sleep(1)
        else:
            if opts.verbose:
                sys.stderr.write('File exists, skip >>> '+target+'\n')
                sys.stderr.flush()
            return fs[target]
    # Upload file
    if opts.verbose:
        tstr = datetime.now()
        sys.stderr.write('{:%Y-%m-%dT%H:%M:%S} Uploading file >>> {}\n'.format(tstr,fnam))
        sys.stderr.flush()
    sender = driver.find_element_by_id('upload-input')
    sender.send_keys(fnam)
    time.sleep(1)
    bit_size = os.path.getsize(fnam)*8
    progress = driver.find_element_by_id('progress')
    bar = progress.find_element_by_xpath('div')
    previous_size = -1
    t1 = tstr
    while True:
        total_size = float(progress.value_of_css_property('width').replace('px',''))
        transfered_size = float(bar.value_of_css_property('width').replace('px',''))
        if transfered_size == 0.0:
            break
        elif transfered_size != previous_size:
            t2 = datetime.now()
            dt = (t2-t1).total_seconds()
            if previous_size < 0:
                ds = transfered_size/total_size
            else:
                ds = (transfered_size-previous_size)/total_size
            if ds > 0.99e-2:
                ds *= bit_size
                remaining_size = bit_size*(total_size-transfered_size)/total_size
                rate = ds/dt
                t3 = t2+timedelta(seconds=remaining_size/rate)
                sys.stderr.write('{:%Y-%m-%dT%H:%M:%S} Uploaded {:6.2f} % @ {:8.3f} Mbps, Expected completion at {:%Y-%m-%dT%H:%M:%S}\n'.format(t2,100.0*transfered_size/total_size,rate*1.0e-6,t3))
                sys.stderr.flush()
                previous_size = transfered_size
                t1 = t2
        time.sleep(1)
    if opts.verbose:
        tend = datetime.now()
        sys.stderr.write('{:%Y-%m-%dT%H:%M:%S} Upload completed in {:.2f} seconds >>> {}\n'.format(tend,(tend-tstr).total_seconds(),fnam))
        sys.stderr.flush()
    # Check uploaded file
    ds,fs = list_file(parent)
    if target in fs:
        return fs[target]
    else:
        sys.stderr.write('Warning, faild in uploading file >>> {}\n'.format(gnam))
        sys.stderr.flush()
        return None

def upload_and_check_file(fnam,gnam):
    title = os.path.basename(gnam)
    size = os.path.getsize(fnam)
    f = upload_file(fnam,gnam)
    if f is not None and (abs(size-f['size']) <= f['size_error']):
        itry = 0
        while f['element'].get_attribute('aria-selected') != 'true':
            itry += 1
            try:
                f['element'].click()
            except Exception:
                pass
            if itry%10 == 0:
                driver.refresh()
            if opts.verbose:
                # for debug
                sys.stderr.write('current_url={}\n'.format(driver.current_url))
                sys.stderr.write('is_displayed={}\n'.format(f['element'].is_displayed()))
                sys.stderr.write('is_enabled={}\n'.format(f['element'].is_enabled()))
                sys.stderr.write('is_selected={}\n'.format(f['element'].is_selected()))
                # for debug
                if itry > 1:
                    sys.stderr.write('Wating for button to be ready >>> {}\n'.format(gnam))
                    sys.stderr.flush()
            time.sleep(1)
        action_info = None
        dropdown = driver.find_element_by_id('dropdown')
        actions = dropdown.find_elements_by_class_name('action')
        for action in actions:
            if action.get_attribute('title') == 'Info':
                action_info = action
                break
        if action_info is None:
            raise ValueError('Error in finding Info.')
        action_info.click()
        time.sleep(1)
        card_contents = driver.find_elements_by_class_name('card-content')
        if len(card_contents) != 1:
            raise ValueError('Error, len(card_contents)={}'.format(len(card_contents)))
        p_md5 = None
        ps = card_contents[0].find_elements_by_tag_name('p')
        for p in ps:
            if re.search('MD5',p.text):
                p_md5 = p
                break
        if p_md5 is None:
            raise ValueError('Error in finding MD5.')
        aa = p_md5.find_elements_by_tag_name('a')
        if len(aa) != 1:
            raise ValueError('Error, len(aa)={}'.format(len(aa)))
        itry = 0
        while aa[0].text == 'Show':
            itry += 1
            aa[0].click()
            if opts.verbose and itry > 1:
                sys.stderr.write('Waiting for MD5 to be ready >>> {}\n'.format(gnam))
                sys.stderr.flush()
            time.sleep(1)
        md5_dst = aa[0].text
        buttons = driver.find_elements_by_class_name('button')
        if len(buttons) != 1:
            raise ValueError('Error, len(buttons)={}'.format(len(buttons)))
        if buttons[0].text != 'OK':
            raise ValueError('Error, buttons[0].text={}'.format(buttons[0].text))
        buttons[0].click()
        time.sleep(1)
        with open(fnam,'rb') as fp:
            md5 = hashlib.md5(fp.read()).hexdigest()
        if (md5_dst.lower() == md5.lower()):
            if opts.verbose:
                sys.stderr.write('Successfully uploaded >>> {}\n'.format(fnam))
                sys.stderr.flush()
            return 0
        else:
            if opts.verbose:
                sys.stderr.write('Warning, size={}, size_error={} ({}), md5={} ({}) >>> {}\n'.format(f['size'],f['size_error'],size,md5_dst,md5,gnam))
                sys.stderr.flush()
            return -1
    elif f is None:
        if opts.verbose:
            sys.stderr.write('Warning, failed in checking file >>> {}\n'.format(gnam))
            sys.stderr.flush()
        return -1
    else:
        if opts.verbose:
            sys.stderr.write('Warning, size={}, size_error={} ({}) >>> {}\n'.format(f['size'],f['size_error'],size,gnam))
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
                os.rmdir(srcdir)
                if opts.debug and not os.path.exists(srcdir):
                    sys.stderr.write('Removed {}\n'.format(srcdir))
                    sys.stderr.flush()
