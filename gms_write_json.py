#!/home/mdi0316/anaconda3/bin/python

### common input start
import os, sys, re
import numpy as np
import pandas as pd
import shutil
import subprocess as sp
import datetime
import time
import math

scripts_dir = '/home/mdi0316/FUNCTIONS'
classes_dir = '/home/mdi0316/CLASSES'
zmat_converter_dir = '/home/mdi0316/CLASSES/zmatrix-master'
sys.path.insert(0, scripts_dir)
sys.path.insert(0, classes_dir)
sys.path.insert(0, zmat_converter_dir)

import GAMESS
import json
from Functions import running_jobs
run_ids, run_job_labels = running_jobs()

def recursive_items(dictionary, depth):
    depth += 1
    for key, value in dictionary.items():
        if type(value) is dict:
            yield from recursive_items(value, depth)
        else:
            yield ('Depth: {}'.format( depth ), key, value)

def read_input( filename ):
    inp_lines = open( filename, 'r' ).readlines()
    inp_dict = {}
    for count, line in enumerate(inp_lines):
        if line.strip().startswith('$'):
           tmp_card = line.replace('$','').strip()
           if not tmp_card == 'END':
              inp_dict[tmp_card] = {}
        else:
           if '=' in line:
              tmp_k, tmp_v = line.split('=')
              inp_dict[tmp_card][tmp_k.strip()] = tmp_v.strip()
           else:
              inp_dict['DATA'][count] = line.split() 
    return( inp_dict ) 

def main():

  run_dir = os.getcwd()
  obj_inp = sys.argv[1]
  inp_dict = read_input( obj_inp )
  
  optimize = False
  hessian = False
  tmp_coords = inp_dict['CONTRL']['COORD'] 
  tmp_runtyp = inp_dict['CONTRL']['RUNTYP']
  if tmp_runtyp == 'OPTIMIZE':
     optimize = True
  elif tmp_runtyp == 'HESSIANY':
     hessian = True


  if 'MPLEVL' in inp_dict['CONTRL'].keys():
     tmp_postscf = 'MP2'
     tmp_postscf_lab = 'MP2'
     mp2 = True
     dft = False
  elif 'DFTTYP' in inp_dict['CONTRL'].keys():
     tmp_postscf = 'DFTTYP'
     tmp_postscf_lab = 'DFT'
     dft = True
     mp2 = False
  else:
     tmp_postscf = 'NONE'
     tmp_postscf_lab = 'NONE'

  if tmp_coords == 'UNIQUE':
     natoms = len(inp_dict['DATA']) - 2
  elif tmp_coords == 'ZMT':
     natoms = int( (6.+float( inp_dict['CONTRL']['NZVAR'] ))/3 )

  obj_calc = GAMESS.GAMESS( inp_name = obj_inp, run_dir = run_dir, natoms = natoms, 
                            run_type = tmp_runtyp, post_scf = tmp_postscf, coordinates = tmp_coords )

  calc_exec, calc_exec_err = obj_calc.get_job_exec()
  calc_inp_dict, calc_out_dict, calc_scf, calc_geom = obj_calc.get_job_results()
  calc_gms_err = obj_calc.read_error()

  print( '\n === EXEC Stauts: {}, EXEC Error: {} === '.format(calc_exec, calc_exec_err))
  print( '\n     GMS Error: {}, SCF: {}, GEOM: {}    '.format( calc_gms_err, calc_scf, calc_geom))
  print( '\n' )
  print( obj_calc )

  if calc_exec == 'TERMINATED.NORMALLY':
    print( 'All keys in dict: {}'.format( calc_out_dict.keys()) )
    for k,v in calc_out_dict.items():
        if k in [ 'ZMAT', 'INTERNUCL.DISTANCES', 'MULL.CHARGES' ]:
           pass
        else:
          if isinstance(v, dict):
             print( '   all keys in {}: {}'.format( k, v.keys() ) )
          else:
             print( k, v )

    with open( 'gms.json', 'w+' ) as json_file:
         json.dump( calc_out_dict, json_file )
  else:
    calc_err = obj_calc.read_error()
    print( '\n === Calculation FAILED: {} === \n'.format(calc_err) )
  
  #if tmp_coords == 'ZMT' and tmp_runtyp == 'OPTIMIZE':
  #   print_zmat = input('Print zmat? (Y)')
  #   if print_zmat == 'Y':
  #      lowest_conf = obj_calc.get_lowest_opt_energy()
  #      for k,v in lowest_conf.items():
  #          print(k,v)
        
  print_ccc = 'Y' #input('Print atomic_coords_and_charges.csv? (Y)')
  print_ind = 'Y' #input('Print inter_nucl_dist.csv? (Y)')

  if print_ccc == 'Y':
     cart_df = pd.DataFrame(columns = [ 'Radius', 'cart.coords.', 'mull.charges' ] ) 
     cart_coords = calc_out_dict['CART.COORDS.']
     mull_charges = calc_out_dict['MULL.CHARGES'] 
  
     ## atomic coordinates and charges
     acac_df = pd.DataFrame( columns=['elem.', 'idx.', 'x', 'y', 'z', 'pop.', 'charge'] )
     
     for cc, mc in zip( cart_coords.values(), mull_charges.values() ):
           acac_df = acac_df.append( pd.Series( { 'elem.' : cc['elem.'], 'idx.' : cc['idx.'], 
                                                  'x': cc['x'], 'y': cc['y'], 'z': cc['z'], 
                                                  'pop.' : mc['pop.'], 'charge': mc['charge'] } ), ignore_index = True )
     acac_df.to_csv( 'atomic_coords_and_charges.csv' )
 
  natoms =  int( (int(inp_dict['CONTRL']['NZVAR']) + 6 )/ 3)
  print( 'NATOMS = {}'.format(natoms) )
  
  if print_ind == 'Y':
     ind_df = pd.DataFrame(columns = [ 'idx1', 'idx2', 'elem.1', 'elem.2', 'distance' ] ) 
     for at_idx_1 in range(1, natoms+1):
         for at_idx_2 in range(1, at_idx_1):
             [ (tmp_k, tmp_v) ]= [ (k,v) for k,v in calc_out_dict['INTERNUCL.DISTANCES'].items() if 
                                 int(v['at.1']['idx.1']) == at_idx_1 and int(v['at.2']['idx.2']) == at_idx_2 ] 
             tmp_dict = {  'idx1'     : tmp_v['at.1']['idx.1'], #'idx{}'.format(v['at.1']['idx.1']),
                           'idx2'     : tmp_v['at.2']['idx.2'], #'idx{}'.format(v['at.2']['idx.2']),
                           'elem.1'   : tmp_v['at.1']['elem.1'],
                           'elem.2'   : tmp_v['at.2']['elem.2'],
                           'distance' : tmp_v['dist.'] }
             ind_df = ind_df.append( pd.Series( tmp_dict ), ignore_index=True )

     ind_df.idx2 = pd.to_numeric( ind_df.idx2 )
     ind_df.idx1 = pd.to_numeric( ind_df.idx1 )
     ind_df = ind_df.sort_values('idx2')
     ind_df = ind_df.sort_values('idx1')
     ind_df = ind_df.reset_index( drop=True )
     ind_df.to_csv( 'inter_nucl_dist.csv' )

if __name__ == '__main__':
  main()

