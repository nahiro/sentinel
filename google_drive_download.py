#!/usr/bin/env python
import os
import sys
import hashlib
from pydrive2.auth import GoogleAuth
from pydrive2.drive import GoogleDrive
from subprocess import call
from optparse import OptionParser,IndentedHelpFormatter

# Default values
HOME = os.environ.get('HOME')
if HOME is None:
    HOME = os.environ.get('HOMEPATH')
DRVDIR = os.path.join(HOME,'Work','SATREPS','IPB_Satreps')
MAX_RETRY = 10

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('-S','--srcdir',default=None,help='GoogleDrive source directory (%default)')
parser.add_option('-D','--dstdir',default=None,help='Local destination directory (%default)')
parser.add_option('--drvdir',default=DRVDIR,help='GoogleDrive directory (%default)')
parser.add_option('-M','--max_retry',default=MAX_RETRY,type='int',help='Maximum number of retries to download data (%default)')
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

qs = [folders[opts.srcdir]['id']]
ds = [os.path.join(opts.dstdir,os.path.basename(opts.srcdir))]
while len(qs) != 0:
    srcdir = qs.pop(0)
    dstdir = ds.pop(0)
    fs = drive.ListFile({'q':'"{}" in parents and trashed = false'.format(srcdir)}).GetList()
    for f in fs:
        dst_nam = os.path.join(dstdir,f['title'])
        if opts.verbose:
            sys.stderr.write(dst_nam+'\n')
            sys.stderr.flush()
        if f['mimeType'] == 'application/vnd.google-apps.folder':
            dnam = dst_nam
            if not os.path.exists(dnam):
                os.makedirs(dnam)
            qs.append(f['id'])
            ds.append(dnam)
        else:
            fnam = dst_nam
            dnam = os.path.dirname(fnam)
            if not os.path.exists(dnam):
                os.makedirs(dnam)
            flag = False
            src_md5 = f['md5Checksum'].upper()
            for ntry in range(opts.max_retry): # loop to download 1 file
                f.GetContentFile(fnam)
                if os.path.exists(fnam):
                    with open(fnam,'rb') as fp:
                        dst_md5 = hashlib.md5(fp.read()).hexdigest().upper()
                    if dst_md5 != src_md5:
                        sys.stderr.write('Warning, dst_md5={}, src_md5={} >>> {}\n'.format(dst_md5,src_md5,fnam))
                        sys.stderr.flush()
                    else:
                        flag = True
                        break
            if not flag:
                sys.stderr.write('Warning, faild in downloading >>> {}\n'.format(fnam))
                sys.stderr.flush()
os.chdir(topdir)
