# python script that populates all the folders of the scratch folder w corresponding input file

import os
import sys
import numpy as np

dir = "/nv/hp22/pfz3/scratch/LH_u_ofhc_cons_3"
file = 'n30d3.txt'
os.chdir(dir)

#bounds = [[80,140],[50,800],[0,70000]]

T = [[float(x) for x in line.split(',')] for line in open(file,'r')]
T = np.array(T)
T[:,2]*=1000

#for i,b in enumerate(bounds):
#	T[:,i] = b[0] + (b[1]-b[0])*T[:,i]

subdirs = [name for name in os.listdir(".") if os.path.isdir(name)]

for i,row in enumerate(T):
        k = i + 1
        if ('%i' %k) not in subdirs:
                os.system('mkdir %i' %k)

os.chdir(r'/nv/hp22/pfz3/scripts/LH_unloading')

print('begining sequential design input file population....')

for i,row in enumerate(T):
         k = 1+i
         print("%s/%i/ E=%1.2f, sig=%1.2f, K=%1.2f, n=%1.2f" %(dir,k,row[0],row[1],row[2],1))
         os.system('python mod_inp_props.py %1.2f %1.2f %1.2f %1.2f' %(row[0],row[1],row[2],1))
         os.system('python gen_batch.py %i %s' %(k,dir))
         os.system('cp batch.pbs "%s/%i"' %(dir,k))
         os.system('cp odbreader.py "%s/%i"' %(dir,k))
         os.system('cp postprocess.py "%s/%i"' %(dir,k))
         os.system('cp LH_u.inp "%s/%i"' %(dir,k))


##############################################################
#  generate bash file to submit all relevant jobs
##############################################################

file = open(os.path.join(dir,'runall.sh'),'w')
file.write('echo "-----begin submitting jobs ------"\n')
for i,row in enumerate(T):
	k=i+1
        file.write('echo "job: %i"\n' %k)
        file.write('cd %i\n' %k)
        file.write('qsub batch.pbs\n')
        file.write('cd ..\n')

file.write('echo "-----end submitting jobs ------\n"')
file.close()
         
