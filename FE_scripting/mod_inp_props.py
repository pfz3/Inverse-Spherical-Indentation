# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 19:22:44 2016

@author: Patxi

take input file and modify the power law hardening paramters
input file is assumed to be contained in the same directory
abaqus however does not have a strain-power-law hardening module
therefore we are simpling providing tabulated values of the hardenling
law

"""
import sys
import numpy as np

raw = []
file = open(r'LH_u.inp','r+')
raw = [line for line in file]

# find ind where *plastic is located
# and find the following line where there is **
# delete all these entries
p_start = next(i+1 for i in range(len(raw)) if '*Plastic' in raw[i])
p_end   = next(i for i in range(p_start,len(raw)) if '**' in raw[i])

del raw[p_start:p_end]

E  = float(sys.argv[1]) # input in GPa
sy = float(sys.argv[2]) # input in MPa
K  = float(sys.argv[3]) # input in MPa
n  = float(sys.argv[4]) # dimensionless

# have to be careful here... discritization is important
# otherwise material won't flow according to power law 
# but rather multi-linear behavior
e_pl        = np.logspace(-300,-1,30000)  
e_pl        = np.hstack(([0],e_pl))
flow_stress = sy + K * e_pl ** n 

el = next(i for i in range(len(raw)) if '*Elastic' in raw[i])
raw[el+1] = '%1.8f,0.3\n' %(E*1e3)

pl = next(i for i in range(el,len(raw)) if '*Plastic' in raw[i])

for i,(p,s) in enumerate(zip(e_pl,flow_stress)):
    raw.insert(pl+i+1,'%1.8f,%1.8e\n' %(s,p))

file = open(r'LH_u.inp','r+')
file.seek(0)
file.write(''.join(raw))
file.truncate()
file.close()

