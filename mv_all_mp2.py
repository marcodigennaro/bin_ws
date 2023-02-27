#!/home/mdi0316/anaconda3/bin/python

import os, shutil 

with open('all_mp2s', 'r') as a:
  lines = a.readlines()

for l in lines:
  if len(l.strip().split('/')) == 13:
     try:
       folder   =  '/'.join( l.strip().split('/')[1:-1]  )
       folder   = os.path.join( '/data/mdi0316/WORK', folder )
       dimer    = l.strip().split('/')[2] 
       opt_from = l.strip().split('/')[6]
       basis    = l.strip().split('/')[8]
       T        = l.strip().split('/')[9] 
       P        = l.strip().split('/')[10] 
       R        = l.strip().split('/')[11] 
       name     = l.strip().split('/')[12] 
       label    = 'gms_SCAN_{}_ENE_MP2_{}_{}_{}_{}'.format( dimer.lower(), basis, T, P, R )
       if name.startswith('log'):
          slurm_id = name.split('_')[-1]
          new_name = 'log.{}_{}'.format(label, slurm_id)
       else:
          if name.endswith('inp'):
             new_name = '{}.inp'.format(label)
          elif name.endswith('dat'):
             new_name = '{}.dat'.format(label)
          else:
             print( explosion )
       os.chdir( folder )
       shutil.move( name, new_name )
     except( FileNotFoundError ):
       print( folder )
       print( name )
       
