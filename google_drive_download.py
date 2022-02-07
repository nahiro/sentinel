#!/usr/bin/env python
import os
import sys
from pydrive2.auth import GoogleAuth
from pydrive2.drive import GoogleDrive
from subprocess import call
from optparse import OptionParser,IndentedHelpFormatter

# Default values
HOME = os.environ.get('HOME')
if HOME is None:
    HOME = os.environ.get('HOMEPATH')
DRVDIR = os.path.join(HOME,'Work','SATREPS','IPB_Satreps')

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('-S','--srcdir',default=None,help='GoogleDrive source directory (%default)')
parser.add_option('-D','--dstdir',default=None,help='Local destination directory (%default)')
parser.add_option('--drvdir',default=DRVDIR,help='GoogleDrive directory (%default)')
parser.add_option('-v','--verbose',default=False,action='store_true',help='Verbose mode (%default)')
parser.add_option('--overwrite',default=False,action='store_true',help='Overwrite mode (%default)')
(opts,args) = parser.parse_args()

opts.dstdir = os.path.abspath(opts.dstdir)
topdir = os.getcwd()
os.chdir(opts.drvdir)

folders = {}

def query_folder(path):
    global folders
    parent = os.path.dirname(path)
    target = os.path.basename(path)
    if (parent == '') or (parent == '/'):
        l = drive.ListFile({'q': '"root" in parents and trashed = false and mimeType = "application/vnd.google-apps.folder" and title = "{}"'.format(target)}).GetList()
        n_list = len(l)
        if n_list != 1:
            os.chdir(topdir)
            raise ValueError('Error in finding folder, n_list={} >>> {}'.format(n_list,path))
        folders.update({path:l[0]})
        return 0
    elif not parent in folders:
        os.chdir(topdir)
        raise IOError('Error, no such folder >>> '+parent)
    l = drive.ListFile({'q': '"{}" in parents and trashed = false and mimeType = "application/vnd.google-apps.folder" and title = "{}"'.format(folders[parent]['id'],target)}).GetList()
    n_list = len(l)
    if n_list > 1:
        os.chdir(topdir)
        raise ValueError('Error in finding folder, n_list={} >>> {}'.format(n_list,path))
    elif n_list == 1:
        if not path in folders:
            folders.update({path:l[0]})
        return 0
    else:
        os.chdir(topdir)
        raise ValueError('Error, no such folder >>> {}'.format(target))
    return -1

gauth = GoogleAuth()
gauth.LocalWebserverAuth()
drive = GoogleDrive(gauth)

l = opts.srcdir.split(os.sep)
for i in range(len(l)):
    if not l[i]:
        continue
    d = os.sep.join(l[:i+1])
    query_folder(d)
    #print(d)

files = {}
qs = [folders[opts.srcdir]['id']]
ds = [opts.dstdir]
cnt = 0
while len(qs) != 0:
    srcdir = qs.pop(0)
    dstdir = ds.pop(0)
    fs = drive.ListFile({'q':'"{}" in parents and trashed = false'.format(srcdir)}).GetList()
    for f in fs:
        files[cnt] = {}
        files[cnt]['id'] = f['id']
        files[cnt]['title'] = f['title']
        files[cnt]['dir'] = os.path.join(dstdir,f['title'])
        if f['mimeType'] == 'application/vnd.google-apps.folder':
            files[cnt]['type'] = 'folder'
            qs.append(f['id'])
            ds.append(files[cnt]['dir'])
        else:
            files[cnt]['type'] = 'file'
        cnt += 1

# This is the base wget command that we will use. This might change in the future due to changes in Google drive
wget_text = '"wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&amp;confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate \'https://docs.google.com/uc?export=download&amp;id=FILE_ID\' -O- | sed -rn \'s/.*confirm=([0-9A-Za-z_]+).*/\\1\\n/p\')&id=FILE_ID" -O FILE_NAME && rm -rf /tmp/cookies.txt"'.replace('&amp;','&')
for cnt in files.keys():
    if files[cnt]['type'] == 'folder':
        os.makedirs(files[cnt]['dir'])
    else:
        command = wget_text[1:-1].replace('FILE_ID',files[cnt]['id']).replace('FILE_NAME',files[cnt]['dir'])
        call(command,shell=True)

os.chdir(topdir)
