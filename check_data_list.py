#!/usr/bin/env python
import os
import sys
import re
import numpy as np
from optparse import OptionParser,IndentedHelpFormatter

# Constants
KB = 1024
MB = KB*1024
GB = MB*1024

# Default values
DATDIR = os.curdir

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.set_usage('Usage: %prog input_fnam [options]')
parser.add_option('-D','--datdir',default=DATDIR,help='Input data directory (%default)')
(opts,args) = parser.parse_args()
if len(args) < 1:
    parser.print_help()
    sys.exit(0)
input_fnam = args[0]

fnams = []
sizes = []
fsizs = []
flags = []
size_labels = []
size_values = []
size_errors = []
with open(input_fnam,'r') as fp:
    for line in fp:
        m = re.search('\s*\d+\s+\S+\s+(S\S+)\s+(\d+)[\s!]*\(\s*([\d\.]+)\s+([KMG]B)\)',line)
        if not m:
            continue
        label = m.group(3)+' '+m.group(4)
        value = float(m.group(3))
        unit = m.group(4).upper()
        if unit == 'KB':
            value *= KB
            error = KB//10
        elif unit == 'MB':
            value *= MB
            error = MB//10
        elif unit == 'GB':
            value *= GB
            error = GB//10
        else:
            raise ValueError('Error in size unit >>> '+line)
        fnams.append(m.group(1))
        sizes.append(int(m.group(2)))
        fsizs.append(-1)
        flags.append(-1)
        size_labels.append(label)
        size_values.append(value)
        size_errors.append(error)
extra_list = []
extra_size = []
for f in sorted(os.listdir(opts.datdir)):
    m = re.search('^S',f)
    if not m:
        continue
    fnam = os.path.join(opts.datdir,f)
    bnam = os.path.splitext(f)[0]
    fsiz = os.path.getsize(fnam)
    if bnam in fnams:
        indx = fnams.index(bnam)
        fsizs[indx] = fsiz
        if (fsiz == sizes[indx]) or (np.abs(fsiz-size_values[indx]) <= size_errors[indx]):
            flags[indx] = 0
        else:
            flags[indx] = 1
    else:
        extra_list.append(f)
        extra_size.append(fsiz)
fnams = np.array(fnams)
sizes = np.array(sizes)
fsizs = np.array(fsizs)
flags = np.array(flags)

cnd = (flags == 0)
sys.stderr.write('OK files: {}\n'.format(cnd.sum()))
cnd = (flags == 1)
sys.stderr.write('NG files: {}\n'.format(cnd.sum()))
for f,s1,s2 in zip(fnams[cnd],sizes[cnd],fsizs[cnd]):
    sys.stderr.write('{} {} {}\n'.format(f,s1,s2))
cnd = (flags == -1)
sys.stderr.write('Missing files: {}\n'.format(cnd.sum()))
for f,s1 in zip(fnams[cnd],sizes[cnd]):
    sys.stderr.write('{} {}\n'.format(f,s1))
sys.stderr.write('Extra files: {}\n'.format(len(extra_list)))
for f,s2 in zip(extra_list,extra_size):
    sys.stderr.write('{} {}\n'.format(f,s2))
