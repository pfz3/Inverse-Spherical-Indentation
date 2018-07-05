from odbAccess import *
import os
import sys
import operator

dir = os.getcwd()
print('Reading Odb in Python...%s' %dir)
odb = openOdb(os.path.join(dir,'sim_%s.odb' %sys.argv[1]))

# get the node label that corresponds to the rigib body reference point
ref = odb.rootAssembly.nodeSets['SET-3'].nodes[0][0].label

TPD_file = open(os.path.join(dir,'indentation.txt')  , "w")

#loop through each frame and get the load and displacement data
steps = [[i,int(i.split('-')[1])] for i in list(odb.steps.keys())]
steps = sorted(steps, key=operator.itemgetter(1,0))
steps = [i[0] for i in steps]
t     = 0.0
for step in steps:
    print('--------- Step %10s ----------' %step)
    N = len(odb.steps[step].frames)
    count = 0
    for frame in odb.steps[step].frames:
        if count<N-1:
            time = float(frame.description.split('=')[-1].strip())
            count+=1 
            print 'processing %i of %i frames' %(count,N)
            node = [node for node in frame.fieldOutputs['RF'].values if node.nodeLabel==ref][0]
            P = -node.data[1]
    
            node = [node for node in frame.fieldOutputs['U'].values if node.nodeLabel==ref][0]
            D = -node.data[1]*1e-3

            # indentation
            TPD_file.write('%s,%s,%s\n' %(t+time,P,D))
    t+=time

closeOdb(odb)
TPD_file.close()
