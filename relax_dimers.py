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
    cation_cart_coords  = mono_dict[cation]['OUT'][gb][fun]['DFT']['OPT']['FINAL']['CART.COORDS.']
    anion_cart_coords   = mono_dict[ anion]['OUT'][gb][fun]['DFT']['OPT']['FINAL']['CART.COORDS.']
    cation_mull_charges = mono_dict[cation]['OUT'][gb][fun]['DFT']['OPT']['FINAL']['MULL.CHARGES']
    anion_mull_charges  = mono_dict[ anion]['OUT'][gb][fun]['DFT']['OPT']['FINAL']['MULL.CHARGES']
    
    nat = cation_nat + anion_nat

    DX = float(R) * np.sin(np.deg2rad(T)) * np.cos(np.deg2rad(P))
    DY = float(R) * np.sin(np.deg2rad(T)) * np.sin(np.deg2rad(P))
    DZ = float(R) * np.cos(np.deg2rad(T))

    label = '{}_{}_T_{}_P_{}_R_{}'.format(cation, anion, T, P, R)
    cart_file = os.path.join( '/data/mdi0316/WORK/COORDS' , '{}_cart.dat'.format(label) )
    zmat_file = os.path.join( '/data/mdi0316/WORK/COORDS' , '{}_zmat.dat'.format(label) )

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
             x = DX + float(coord['x'])
             y = DY + float(coord['y'])
             z = DZ + float(coord['z'])
             df.write( '{} {} {} {}\n'.format(elem.upper(), str(x), str(y), str(z)) )
    
    ## zmatrix file
    c = Converter()
    c.run_cartesian( cart_file )
    shutil.move( 'zmatrix.dat', zmat_file )
    
    ## zmatrix dict
    zmat_dict = {}
    for count, line in enumerate(open( zmat_file, 'r' ).readlines()[2:]):
        if len(line.split()) == 1:
           elem = line.split()[0]
           zmat_dict[count] = { 'elem.' : elem , 'idx.' : count + 1 }
        elif len(line.split()) == 3:
           elem, ref_str, ref_val = line.split()
           zmat_dict[count] = { 'elem.' : elem , 'idx.' : count + 1, 
                                'STR'   : { 'ref' : ref_str, 'val' : ref_val }  }
        elif len(line.split()) == 5:
           elem, ref_str, ref_val, ben_str, ben_val = line.split()
           zmat_dict[count] = { 'elem.' : elem , 'idx.' : count + 1, 
                                'STR'   : { 'ref' : ref_str, 'val' : ref_val }, 
                                'BEN'   : { 'ref' : ben_str, 'val' : ben_val }  }
        elif len(line.split()) == 7:
           elem, ref_str, ref_val, ben_str, ben_val, tor_str, tor_val = line.split()
           zmat_dict[count] = { 'elem.' : elem , 'idx.' : count + 1, 
                                'STR'   : { 'ref' : ref_str, 'val' : ref_val }, 
                                'BEN'   : { 'ref' : ben_str, 'val' : ben_val },
                                'TOR'   : { 'ref' : tor_str, 'val' : tor_val }  }
    return( zmat_dict, cart_file )

post_scf = 'NONE'
post_scf = 'DFTTYP'

def main():
    equil_dir = os.path.join( work_dir, '{}_{}/RUNS/EQUILIBRIUM'.format(cation,anion) )
    projct = True

    #equil_dir = os.path.join( work_dir, '{}_{}/RUNS/EQUILIBRIUM_projct_f'.format(cation,anion) )
    #projct = False

    tpr_dict = { 
                 '0' : { 'T' :   '5', 'P' :   '0' , 'R' : '10.0' },
                 '1' : { 'T' :  '90', 'P' :   '0' , 'R' : '10.0' },
                 '2' : { 'T' :  '90', 'P' :  '90' , 'R' : '10.0' },
                 '3' : { 'T' :  '90', 'P' : '180' , 'R' : '10.0' },
                 '4' : { 'T' :  '90', 'P' : '270' , 'R' : '10.0' },
                 '5' : { 'T' : '175', 'P' :   '0' , 'R' : '10.0' },
                 '6' : { 'T' :   '5', 'P' :   '0' , 'R' : '10.0' },
                 '7' : { 'T' :   '5', 'P' :   '0' , 'R' :  '7.0' },
                 '8' : { 'T' :  '90', 'P' :   '0' , 'R' :  '7.0' },
                 '9' : { 'T' :  '90', 'P' :  '90' , 'R' :  '7.0' },
                 '0' : { 'T' :  '90', 'P' : '180' , 'R' :  '7.0' },
                '11' : { 'T' :  '90', 'P' : '270' , 'R' :  '7.0' },
                '12' : { 'T' : '175', 'P' :   '0' , 'R' :  '7.0' },
                '13' : { 'T' :   '5', 'P' :   '0' , 'R' :  '7.0' }
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

      obj_label = '{}_{}_T_{}_P_{}_R_{}'.format(cation, anion, T, P, R)
      cart_label = 'cart_{}'.format( obj_label )
      print( obj_label )
     
      if len( [ ii for ii in run_job_labels if obj_label in ii ] ) != 0:
         print_tab( 1, 'running' )

      else:
         obj_calc = GAMESS.GAMESS_calculation( obj_label, 
                                               root_dir = root_dir, natoms = 24, 
                                               runtyp = run_typ, post_scf = post_scf, 
                                               basis = gb, functional = fun, projct = projct )
         cart_calc = GAMESS.GAMESS_calculation( cart_label, 
                                                root_dir = root_dir, natoms = 24, 
                                                runtyp = run_typ, post_scf = post_scf, 
                                                #verbose = True 
                                                basis = gb, functional = fun, coordinates = 'UNIQUE', projct = projct ) 
         if os.path.exists( cart_calc.run_dir ):
            cart_exec = cart_calc.get_execution( run_job_labels )
            cart_out_dict = cart_calc.get_out_dict()
            print_tab( 1, "cart : {}".format(cart_exec) )
            print_tab( 1, "cart : {}".format(cart_out_dict['DFT']['OPT']['FINAL']['SCF']) )
            print_tab( 1, "cart : {}".format(cart_out_dict['DFT']['OPT']['FINAL']['GEOM.']) )
            if os.path.exists( obj_calc.run_dir ):
               obj_exec = obj_calc.get_execution( run_job_labels )
               obj_out_dict = obj_calc.get_out_dict()
               print_tab( 1, "zmat : {} ".format(obj_exec) )
               print_tab( 1, "zmat : {} ".format(obj_out_dict['DFT']['OPT']['FINAL']['SCF']) )
               print_tab( 1, "zmat : {} ".format(obj_out_dict['DFT']['OPT']['FINAL']['GEOM.']) )
            else:
               ## make zmat_relaxation from cartesian result
               relaxed_cart_coords = cart_out_dict['DFT']['OPT']['FINAL']['CART.COORDS.']
               relaxed_zmat_dict = cart_calc.cart_to_file( relaxed_cart_coords )
               obj_calc.write_input_file_ZMAT( relaxed_zmat_dict )
               slurm_obj = SLURM.SLURM( obj_calc.run_dir , 'GAMESS', job_name = obj_calc.inp_name, job_queue = 'nodeshiq' )
               slurm_obj.write_batch()
               slurm_obj.submit_batch()
         else:
            ## make cart_relaxation
            zmat_dict, cart_file = make_coordinates( float(T), float(P), float(R), anion, cation )
            cart_calc.write_input_file_CART( cart_file )
            slurm_obj = SLURM.SLURM( cart_calc.run_dir , 'GAMESS', job_name = cart_calc.inp_name, job_queue = 'nodesloq' )
            slurm_obj.write_batch()
            slurm_obj.submit_batch()


if __name__ == '__main__':
  main()

