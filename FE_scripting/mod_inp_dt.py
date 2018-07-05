# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 19:22:44 2016

@author: Patxi

sys.argv[1] = minimum dt time step to be used in analysis
"""
import sys

raw = []
file = open(r'rambergosgood.inp','r+')
raw = [line for line in file]


text_prev = ''
for i in range(len(raw)):
    line = raw[i]
    if '*Static' in text_prev:
        raw[i] = '0.01,100.,1e-05,%1.6e\n' %(float(sys.argv[1]))
    text_prev = line

file = open(r'rambergosgood.inp','r+')
file.seek(0)
file.write(''.join(raw))
file.truncate()
file.close()

