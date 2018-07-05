import os
import sys

home = r'/nv/hp22/pfz3/scratch/'

print('deleting lck files in the following directories:')
for i in sys.argv[1:]:
	print('              %s%s' %(home,i))
os.system('cd %s%s' %(home,sys.argv[1]))
print('changing directories...')
print(os.getcwd())
if '%s.lck' in os.listdir(os.getcwd()):
	print('Lock file being removed...') 
	os.system('rm %s.lck' %sys.argv[1])
else:
	print('Lock file not detected')
os.system('cd %s' %home)


