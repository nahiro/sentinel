#!/usr/bin/env python
import sys
import re
from optparse import OptionParser,IndentedHelpFormatter

# Default values
BLK_FNAM = 'field_block.dat'
AREA_FNAM = 'pixel_area.dat'

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('-b','--blk_fnam',default=BLK_FNAM,help='Block file name (%default)')
parser.add_option('-a','--area_fnam',default=AREA_FNAM,help='Pixel area file name (%default)')
(opts,args) = parser.parse_args()

blocks = []
with open(opts.blk_fnam,'r') as fp:
    iobj = 0
    for line in fp:
        item = line.split()
        if len(item) !=  3:
            raise ValueError('Error, len(item)={}'.format(len(item)))
        object_id = int(item[0])
        block = item[1]
        if object_id != iobj+1:
            raise ValueError('Error, object_id={}, iobj={}'.format(object_id,iobj))
        blocks.append(block)
        iobj += 1

with open(opts.area_fnam,'r') as fp:
    iobj = 0
    for line in fp:
        m = re.search('^\s*(\d+)\s+(\S+.*$)',line)
        if not m:
            raise ValueError('ERROR')
        object_id = int(m.group(1))
        if object_id != iobj+1:
            raise ValueError('Error, object_id={}, iobj={}'.format(object_id,iobj))
        data = m.group(2)
        sys.stdout.write('{} {} {}\n'.format(object_id,blocks[iobj],data))
        #break
        iobj += 1

if iobj != len(blocks):
    raise ValueError('Error, iobj={}, len(blocks)={}'.format(iobj,len(blocks)))
