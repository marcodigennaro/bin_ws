#!/home/mdi0316/anaconda3/bin/python

import os
import numpy as np
import shutil
import subprocess as sp

def find_log_in_dir( folder ):
    if not os.path.exists( folder ):
       return False
    else:
       logs_list = [ f for f in os.listdir(folder) if f.startswith('log') ]
       if len(logs_list) == 0:
          return False
       elif len(logs_list) == 1:
          return logs_list[0]
       else:
          print( too_many_log_file_present )
          return 0

main_dir='/data/mdi0316/WORK/TEST/NWCHEM/BF4'
template='bf4.nw'
template_lines = open( template, 'r' ).readlines()
slurm='/home/mdi0316/SCRIPTS/submit_nwchem.sh'

for ecut in range(5, 50, 2):
  new_dir = os.path.join( main_dir, 'ecut_convergence', 'ecut_{}'.format(ecut) )
  new_lab = 'bf4_ecut_{}'.format(ecut) 
  new_log = find_log_in_dir( new_dir )
  
  if new_log == False:
     os.makedirs( new_dir, exist_ok = True )
     os.chdir(new_dir)
     with open( '{}.nw'.format(new_lab), 'w+' ) as inp:
       for line in template_lines:
         new_line = line
         if new_line.strip().startswith( 'energy_cutoff' ):
            new_line = ' energy_cutoff   {}\n'.format(ecut)
         inp.write(new_line)
     shutil.copy( slurm, new_dir )
     sp.call( 'sbatch -J {} submit_nwchem.sh {}.nw '.format(new_lab, new_lab), shell=True )
  elif new_log:
     print(new_log)


os.chdir(main_dir)
ecut_list = [ int(ecut.replace('ecut_','')) for ecut in os.listdir( 'ecut_convergence' )] 
ecut_list.sort()

