#!/home/mdi0316/anaconda3/bin/python

import os
import subprocess as sp

with open( 'to_del', 'r' ) as f:
  lines = f.readlines()

fold = os.getcwd()

for l in lines:
  print(l.strip())
  basis = l.strip().split('/')[7]
  funct = l.strip().split('/')[8]
  T = l.strip().split('/')[9] 
  P = l.strip().split('/')[10] 
  R = l.strip().split('/')[11] 

  old_inp = [ x for x in os.listdir(l.strip()) if x.endswith('inp') ][0]
  old_dat = [ x for x in os.listdir(l.strip()) if x.endswith('dat') ][0]
  old_log = [ x for x in os.listdir(l.strip()) if x.startswith('log') ][0]

  old_slurm = old_log.split('_')[-1]
  label = 'gms_SCAN_emim_bf4_ENE_DFT'
  new_inp = '{}_{}_{}_{}_{}_{}.inp'.format( label, basis, funct, T, P, R )
  new_dat = '{}_{}_{}_{}_{}_{}.dat'.format( label, basis, funct, T, P, R )
  new_log = 'log.{}_{}_{}_{}_{}_{}_{}'.format( label, basis, funct, T, P, R, old_slurm )


  os.chdir( l.strip() )
  sp.call( 'mv {} {}'.format(old_inp, new_inp), shell=True )
  sp.call( 'mv {} {}'.format(old_dat, new_dat), shell=True )
  sp.call( 'mv {} {}'.format(old_log, new_log), shell=True )
  os.chdir( fold )
