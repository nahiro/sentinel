#!/usr/bin/env python
import sys
import re

blocks = []
with open('field_block.dat','r') as fp:
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

with open('pixel_area.dat','r') as fp:
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
