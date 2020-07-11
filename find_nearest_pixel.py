#!/usr/bin/env python
import sys
import numpy as np

xstp = 10.0
ystp = -10.0
xmin,xmax,ymin,ymax = (743800.0,756800.0,9236000.0,9251800.0)
xg,yg = np.meshgrid(np.arange(xmin,xmax+0.1*xstp,xstp),np.arange(ymax,ymin-0.1*ystp,ystp))
ngrd = xg.size
nx = xg.shape[1]
ny = xg.shape[0]
sid = np.arange(ngrd).reshape(xg.shape)

sid_0 = []
sid_1 = []
sid_2 = []
sid_3 = []
sid_4 = []
sid_5 = []
sid_6 = []
sid_7 = []
sid_8 = []
sid_9 = []
sid_a = []
sid_b = []
sid_c = []
leng_1 = []
leng_2 = []
leng_3 = []
leng_4 = []
leng_5 = []
leng_6 = []
leng_7 = []
leng_8 = []
leng_9 = []
leng_a = []
leng_b = []
leng_c = []
for n in range(ngrd):
    i,j = np.unravel_index(n,xg.shape)
    x1 = max(0,j-2)
    x2 = min(nx,x1+5)
    x1 = x2-5
    y1 = max(0,i-2)
    y2 = min(ny,y1+5)
    y1 = y2-5
    l2 = (np.square(xg[y1:y2,x1:x2]-xg[i,j])+np.square(yg[y1:y2,x1:x2]-yg[i,j])).flatten()
    inds = sid[y1:y2,x1:x2].flatten()
    indx = np.argsort(l2)
    sid_0.append(n)
    sid_1.append(inds[indx[1]])
    sid_2.append(inds[indx[2]])
    sid_3.append(inds[indx[3]])
    sid_4.append(inds[indx[4]])
    sid_5.append(inds[indx[5]])
    sid_6.append(inds[indx[6]])
    sid_7.append(inds[indx[7]])
    sid_8.append(inds[indx[8]])
    sid_9.append(inds[indx[9]])
    sid_a.append(inds[indx[10]])
    sid_b.append(inds[indx[11]])
    sid_c.append(inds[indx[12]])
    leng_1.append(np.sqrt(l2[indx[1]]))
    leng_2.append(np.sqrt(l2[indx[2]]))
    leng_3.append(np.sqrt(l2[indx[3]]))
    leng_4.append(np.sqrt(l2[indx[4]]))
    leng_5.append(np.sqrt(l2[indx[5]]))
    leng_6.append(np.sqrt(l2[indx[6]]))
    leng_7.append(np.sqrt(l2[indx[7]]))
    leng_8.append(np.sqrt(l2[indx[8]]))
    leng_9.append(np.sqrt(l2[indx[9]]))
    leng_a.append(np.sqrt(l2[indx[10]]))
    leng_b.append(np.sqrt(l2[indx[11]]))
    leng_c.append(np.sqrt(l2[indx[12]]))
    #break
sid_0 = np.array(sid_0)
sid_1 = np.array(sid_1)
sid_2 = np.array(sid_2)
sid_3 = np.array(sid_3)
sid_4 = np.array(sid_4)
sid_5 = np.array(sid_5)
sid_6 = np.array(sid_6)
sid_7 = np.array(sid_7)
sid_8 = np.array(sid_8)
sid_9 = np.array(sid_9)
sid_a = np.array(sid_a)
sid_b = np.array(sid_b)
sid_c = np.array(sid_c)
leng_1 = np.array(leng_1)
leng_2 = np.array(leng_2)
leng_3 = np.array(leng_3)
leng_4 = np.array(leng_4)
leng_5 = np.array(leng_5)
leng_6 = np.array(leng_6)
leng_7 = np.array(leng_7)
leng_8 = np.array(leng_8)
leng_9 = np.array(leng_9)
leng_a = np.array(leng_a)
leng_b = np.array(leng_b)
leng_c = np.array(leng_c)
np.savez('find_nearest.npz',sid_0=sid_0,
sid_1=sid_1,leng_1=leng_1,
sid_2=sid_2,leng_2=leng_2,
sid_3=sid_3,leng_3=leng_3,
sid_4=sid_4,leng_4=leng_4,
sid_5=sid_5,leng_5=leng_5,
sid_6=sid_6,leng_6=leng_6,
sid_7=sid_7,leng_7=leng_7,
sid_8=sid_8,leng_8=leng_8,
sid_9=sid_9,leng_9=leng_9,
sid_a=sid_a,leng_a=leng_a,
sid_b=sid_b,leng_b=leng_b,
sid_c=sid_c,leng_c=leng_c)
