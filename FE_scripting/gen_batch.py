import os
import sys

dir = sys.argv[2]

file = open('batch.pbs','w')
file.write('# batch pbs file\n')
file.write('#PBS -N sim_%s\n' %sys.argv[1])
file.write('#PBS -l nodes=1:ppn=2\n')
file.write('#PBS -l mem=2gb\n')
file.write('#PBS -l walltime=150:00:00\n')
file.write('#PBS -q prometheus\n')
file.write('#PBS -j oe\n')
file.write('#PBS -o sim_%s.out\n' %sys.argv[1])
file.write('#PBS -m abe\n')
file.write('#PBS -M patxifernandezzelaia@gmail.com\n')
file.write('#PBS -V\n')
file.write('\n')
file.write('module load abaqus/6.16\n')
file.write('module load anaconda3/4.1.1\n')
file.write('\n')
file.write('cd "%s/%s"\n' %(dir,sys.argv[1]))
file.write('echo "cd = $PWD"\n')
file.write('if [ -f "sim_%s.lck" ];\n' %sys.argv[1])
file.write('then\n')
file.write('    echo "deleting sim_%s.lck file"\n' %sys.argv[1])
file.write('    rm sim_%s.lck\n' %sys.argv[1])
file.write('else\n')
file.write('    echo "sim_%s.lck not detected"\n' %sys.argv[1])
file.write('fi\n')
file.write('\n')
file.write('echo "------ Begining Abaqus FEA ------"\n')
file.write("abaqus cpus=2 memory='2 gb' job=sim_%s interactive input=LH_u.inp\n" %sys.argv[1])
file.write('echo "------ End of Abaqus FEA   ------"\n')
file.write('\n')
file.write('abaqus python odbreader.py %s\n' %sys.argv[1])
file.write('python postprocess.py %s\n' %sys.argv[1])
file.close()
