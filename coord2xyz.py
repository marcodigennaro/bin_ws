#!/home/mdi0316/anaconda3/bin/python

import sys 
bhor2ang=0.529
coord_file = sys.argv[1]
with open(coord_file, 'r') as f:
  read_lines = f.readlines()

for cc, rl in enumerate(read_lines):
  if rl.startswith('$coord'):
     start = cc
  elif rl.startswith('$'):
     end = cc
     break

xyz_file = coord_file.split('.')[0] + '.xyz'
with open( xyz_file , 'w+' ) as f:
  f.write(f'{end-start-1}\n#xyz file from turbomole output\n')
  for rl in read_lines[start+1:end]:
    x,y,z,e = rl.strip().split()
    f.write( f'{e} {bhor2ang*float(x)} {bhor2ang*float(y)} {bhor2ang*float(z)}\n' )

