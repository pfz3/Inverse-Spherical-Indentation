# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 19:03:09 2016

@author: Patxi
"""
import os
import re
from scipy.optimize import least_squares
import numpy as np

def func_unload(x,ht,P,E,):
    
    # ht and P must be equal length (and shape) np arrays
    return ht - x[0] - ((1-0.3**2)*3*P/4/E/(x[1]*1e3)**0.5)**(2.0/3.0)

raw = []


dir = os.getcwd()

for ii in range(1,51):
    print('Post Processing in Python... ind %i' %ii)
    with open(os.path.join(dir,'%i_indentation.txt' %ii) ,'r') as f:
        raw=[[float(re.sub('\n','',j)) for j in line.split(',')] for line in f]
        
    raw = np.array(raw)
    t = raw[:,0]
    P = raw[:,1]
    D = raw[:,2]    

    unload_inds = [i for i in range(1,raw.shape[0]-1) if D[i]>D[i-1] and D[i]>D[i+1]]


    # get it into Pa
    E = float([[float(i) for i in re.sub('\n','',line).split(',')] for line in open('n50d3.txt','r')][ii-1][0])*1e9
    print('E read from input file... %1.2f GPa' %(E*1e-9))

    au = [0]
    x0=[0.0,6350] #plastic displacement and initial contact radius
    xv = []
    R = []
    count = 0
    for start in unload_inds:
        try:
            end = next(i for i in range(start+1,raw.shape[0]-1) if D[i+1]>D[i])+1
            count+=1
        except:
            end = raw.shape[0]-1
            break
        bounds = [(0,0),(20000,np.inf)]
        x = least_squares(func_unload,x0,bounds=bounds,args=(D[start:end]*1e9,P[start:end]*1e3,E*1e-15))
        x0=x['x']    
        xv.append(x0)
        au.append(((1-0.3**2)*3*P[start]*x['x'][1]*1e-6/4/(E))**(1.0/3.0))
        R.append(x['cost'])
    au = np.array(au)
    Du= D[[0]+unload_inds[:count]]
    tu= t[[0]+unload_inds[:count]]
    Pu= P[[0]+unload_inds[:count]]

    text_file = open(os.path.join(dir,'%i_output.txt' %ii) ,'w')
    for x,y,z,w in zip(tu,Pu,Du,au):
        text_file.write('%s,%s,%s,%s\n' %(x,y,z,w))
    text_file.close()
