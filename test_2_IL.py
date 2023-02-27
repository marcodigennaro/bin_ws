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

import filecmp
import warnings

import GAMESS
import IONIC_LIQUID
import SLURM

from converter import Converter
from Functions import running_jobs, print_tab
import numpy.linalg as LA
import scipy.constants as const
Ha2eV =  const.value('hartree-electron volt relationship') #27.211
Ang2Bohr = 1.8897259886

from IONIC_LIQUID import mono_dict, complete_R_LIST
from GAMESS import functionals_list, gbasis_list, polar_dict, ispher_1_list

work_dir = '/data/mdi0316/WORK'
scan_dir = os.path.join( work_dir, 'EMIM_BF4/RUNS/SCAN' )
mono_dir = os.path.join( work_dir, 'MONOMERS' )

import json
mono_dict_file = '/home/mdi0316/Inputfiles/GAMESS/monomers.json'
with open(mono_dict_file,'r') as json_file:  
    mono_dict = json.load(json_file)

run_ids, run_job_labels = running_jobs()

anion = 'BF4'
cation = 'EMIM'

gb = 'N311'
fun = 'B3LYP'

run_typ = 'OPTIMIZE'

def make_coordinates( T, P, R, anion, cation ):
    cation_nat = mono_dict[cation]['nat']
    anion_nat = mono_dict[anion]['nat']
    cation_cart_coords = mono_dict[cation]['OUT'][gb][fun]['DFT']['OPT']['CART.COORDS.']
    anion_cart_coords = mono_dict[anion]['OUT'][gb][fun]['DFT']['OPT']['CART.COORDS.']
    cation_mull_charges = mono_dict[cation]['OUT'][gb][fun]['DFT']['OPT']['MULL.CHARGES']
    anion_mull_charges = mono_dict[anion]['OUT'][gb][fun]['DFT']['OPT']['MULL.CHARGES']
    
    nat = cation_nat + anion_nat

    DX1 = float(R) * np.sin(np.deg2rad(T)) * np.cos(np.deg2rad(P))
    DY1 = float(R) * np.sin(np.deg2rad(T)) * np.sin(np.deg2rad(P))
    DZ1 = float(R) * np.cos(np.deg2rad(T))

    T2, P2, R2 = 90, 180, 20.
    DX2 = float(R2) * np.sin(np.deg2rad(T2)) * np.cos(np.deg2rad(P2))
    DY2 = float(R2) * np.sin(np.deg2rad(T2)) * np.sin(np.deg2rad(P2))
    DZ2 = float(R2) * np.cos(np.deg2rad(T2))


    label = '{}_{}_T_{}_P_{}_R_{}'.format(cation, anion, T, P, R)
    cart_file = os.path.join( '/data/mdi0316/WORK/COORDS' , '2couples_{}_cart.dat'.format(label) )
    zmat_file = os.path.join( '/data/mdi0316/WORK/COORDS' , '2couples_{}_zmat.dat'.format(label) )

    ## cartesian file
    with open(cart_file, 'w+') as df:
         df.write( '{}\n\n'.format(cation_nat+anion_nat) )
         for kk, coord in cation_cart_coords.items():
             elem = coord['elem.']
             x = coord['x']
             y = coord['y']
             z = coord['z']
             df.write( '{} {} {} {}\n'.format(elem.upper(), x, y, z) )
         for kk, coord in anion_cart_coords.items():
             elem = coord['elem.']
             x = DX1 + float(coord['x'])
             y = DY1 + float(coord['y'])
             z = DZ1 + float(coord['z'])
             df.write( '{} {} {} {}\n'.format(elem.upper(), str(x), str(y), str(z)) )
         for kk, coord in cation_cart_coords.items():
             elem = coord['elem.']
             x = DX2 + float(coord['x'])
             y = DY2 + float(coord['y'])
             z = DZ2 + float(coord['z'])
             df.write( '{} {} {} {}\n'.format(elem.upper(), x, y, z) )
         for kk, coord in anion_cart_coords.items():
             elem = coord['elem.']
             x = DX2 + DX1 + float(coord['x'])
             y = DY2 + DY1 + float(coord['y'])
             z = DZ2 + DZ1 + float(coord['z'])
             df.write( '{} {} {} {}\n'.format(elem.upper(), str(x), str(y), str(z)) )
    
    return( cart_file )

post_scf = 'NONE'
post_scf = 'DFTTYP'

def main():
    equil_dir = os.path.join( work_dir, '{}_{}/RUNS_2_IL/EQUILIBRIUM'.format(cation,anion) )

    tpr_dict = { 
                 '0' : { 'T' :   '5', 'P' :   '0' , 'R' : '10.0' },
                 '1' : { 'T' :  '90', 'P' :   '0' , 'R' : '10.0' },
                 '2' : { 'T' :  '90', 'P' :  '90' , 'R' : '10.0' },
                 '3' : { 'T' :  '90', 'P' : '180' , 'R' : '10.0' },
                 '4' : { 'T' :  '90', 'P' : '270' , 'R' : '10.0' },
                 '5' : { 'T' : '175', 'P' :   '0' , 'R' : '10.0' },
                 '6' : { 'T' :   '5', 'P' :   '0' , 'R' : '10.0' }
                  }

    for tpr_item in tpr_dict.values():
      T = tpr_item['T']
      P = tpr_item['P']
      R = tpr_item['R']

      root_dir  = os.path.join( equil_dir, 'T_{}'.format(T),'P_{}'.format(P),'R_{}'.format(R))
    
      if post_scf == 'NONE': 
         post_label = 'NONE'
      else:
         post_label = post_scf[:3]

      obj_label = '2couples_{}_{}_T_{}_P_{}_R_{}'.format(cation, anion, T, P, R)
     
      if len( [ ii for ii in run_job_labels if obj_label in ii ] ) != 0:
         print_tab( 1, '{} running'.format(obj_label) )

      else:
         obj_calc = GAMESS.GAMESS_calculation( obj_label, 
                                               root_dir = root_dir, natoms = 24, 
                                               runtyp = run_typ, post_scf = post_scf, 
                                               basis = gb, functional = fun )
         print( obj_calc )
         if os.path.isdir( obj_calc.run_dir ):
            pass
         else:
           cart_label = 'cart_{}'.format( obj_label )
           cart_calc = GAMESS.GAMESS_calculation( cart_label, 
                                                  root_dir = root_dir, natoms = 48, 
                                                  runtyp = run_typ, post_scf = post_scf, 
                                                  basis = gb, functional = fun, coordinates = 'UNIQUE', 
                                                  verbose = True 
                                                   )
           if os.path.isdir( cart_calc.run_dir ):
              
              cart_out_dict = cart_calc.get_out_dict()
              print( cart_out_dict.keys() )
           else:
              ## make cart_relaxation
              cart_file = make_coordinates( float(T), float(P), float(R), anion, cation )
              cart_calc.write_input_file_CART( cart_file )
              slurm_obj = SLURM.SLURM( cart_calc.run_dir , 'GAMESS', job_name = cart_calc.inp_name, job_queue = 'nodeshiq' )
              slurm_obj.write_batch()
              slurm_obj.submit_batch()
  #            print( afdsp)


if __name__ == '__main__':
  main()

