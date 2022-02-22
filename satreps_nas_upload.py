#!/usr/bin/env python
import os
import sys
import re
import shutil
import hashlib
import atexit
import time
from glob import glob
from pathlib import Path
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.keys import Keys
from optparse import OptionParser,IndentedHelpFormatter

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
os.chdir(opts.drvdir)
def clean_up():
    sys.stderr.write('clean_up called.\n')
    sys.stderr.flush()
    os.chdir(topdir)
atexit.register(clean_up)

def list_file(path):
    ds = []
    fs = []
    if len(path) < 1:
        url = 'https://{}/files/'.format(server)
    elif path[0] == '/':
        url = 'https://{}/files{}'.format(server,path)
    else:
        url = 'https://{}/files/{}'.format(server,path)
    driver.get(url)
    time.sleep(1)
    if re.search('This location can\'t be reached',driver.page_source):
        return ds,fs
    items = driver.find_elements_by_class_name('item')
    for item in items:
        lines = item.text.splitlines()
        nline = len(lines)
        if nline == 3:
            if lines[0] == 'folder':
                ds.append(lines[1])
            elif lines[0] == 'insert_drive_file':
                fs.append(lines[1])
            else:
                raise ValueError('Error, lines[0]={}'.format(lines[0]))
        elif nline == 2:
            fs.append(lines[0])
    return ds,fs

folders = []

def make_folder(path):
    global folders
    if path in folders:
        return 0
    parent = os.path.dirname(path)
    target = os.path.basename(path)
    ds,fs = list_file(parent)
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
        folders.append(path)
        return 0
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

if opts.port is not None:
    options = Options()
    options.add_experimental_option('debuggerAddress','localhost:{}'.format(opts.port))
    driver = webdriver.Chrome(os.path.join(opts.drvdir,'chromedriver'),options=options)
else:
    driver = webdriver.Chrome(os.path.join(opts.drvdir,'chromedriver'))

# Login to the NAS server
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
url = 'https://{}/files/{}'.format(server,opts.dstdir)
driver.get(url)
time.sleep(1)
if re.search('This location can\'t be reached.',driver.page_source):
    raise IOError('No such page >>> '+url)
"""
sender = driver.find_element_by_id('upload-input')
sender.send_keys(os.path.join(opts.srcdir,'test.txt'))
time.sleep(1)
progress = driver.find_element_by_id('progress')
bar = progress.find_element_by_xpath('div')
while True:
    total_size = float(progress.value_of_css_property('width').replace('px',''))
    transfered_size = float(bar.value_of_css_property('width').replace('px',''))
    if transfered_size == 0.0:
        break
    print(transfered_size,'/',total_size)
items = driver.find_elements_by_class_name('item')
for item in items:
    lines = item.text.splitlines()
    if (len(lines) > 2) and lines[1] == 'test.txt':
        item.click()
dropdown = driver.find_element_by_id('dropdown')
actions = dropdown.find_elements_by_class_name('action')
actions[8].click()
card_content = driver.find_element_by_class_name('card-content')
ps = card_content.find_elements_by_tag_name('p')
a = ps[4].find_element_by_tag_name('a')
a.click()
while True:
    dst_md5 = a.text
    if dst_md5 == 'Show':
        time.sleep(1)
        continue
    break
dst_size = ps[1].text
"""
