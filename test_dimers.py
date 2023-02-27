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
from Functions import running_jobs
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

mono_dict_file = '/home/mdi0316/Inputfiles/GAMESS/monomers.json'
with open(mono_dict_file,'r') as json_file:
    mono_dict = json.load(json_file)

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

    DX = float(R) * np.sin(np.deg2rad(T)) * np.cos(np.deg2rad(P))
    DY = float(R) * np.sin(np.deg2rad(T)) * np.sin(np.deg2rad(P))
    DZ = float(R) * np.cos(np.deg2rad(T))

    ## cartesian file
    label = '{}_{}_T_{}_P_{}_R_{}'.format(cation, anion, T, P, R)
    cart_file = os.path.join( '/data/mdi0316/WORK/COORDS' , '{}_cart.dat'.format(label) )
    zmat_file = os.path.join( '/data/mdi0316/WORK/COORDS' , '{}_zmat.dat'.format(label) )
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
#post_scf = 'DFTTYP'

def main():
    scan_dir = os.path.join( work_dir, '{}_{}/RUNS/SCAN'.format(cation,anion) )
    T,P,R = 90, 270, 15.0
    root_dir  = os.path.join( scan_dir, 'T_{}'.format(T),'P_{}'.format(P),'R_{}'.format(R))
    #opt_inp = '{}_{}_{}_{}_{}_R_{}_T_{}_P_{}.inp'.format(cation, anion, gb, fun, post_typ[:3], run_typ[:3], R, T, P)
    if post_scf == 'NONE': 
       post_label = 'HF'
    else:
       post_label = post_scf[:3]
    obj_inp = '{}_{}_{}_T_{}_P_{}_R_{}_OPT.inp'.format(cation, anion, post_label, T, P, R)
    obj_calc = GAMESS.GAMESS_calculation( obj_inp, root_dir = root_dir, natoms = 24, 
                                         runtyp = run_typ, post_scf = post_scf, 
                                         basis = gb, functional = fun )
    if os.path.isdir( obj_calc.run_dir ):
       if obj_calc.out_file!= None:
          obj_exec = obj_calc.get_execution( run_job_labels )
          if obj_exec == 'NORMALLY':
            print( make_DFT )
            #print( 'Reading: {}'.format(obj_calc.out_file) )
            obj_calc_dict = obj_calc.get_out_dict()['DFT']['OPT']
            all_charges = obj_calc_dict['MULL.CHARGES'] 
            cat_charges = sum([ v['charge'] for (k,v) in all_charges.items() if k < 19 ] )
            ani_charges = sum([ v['charge'] for (k,v) in all_charges.items() if k >= 19 ] )
            opt_res_dict[ ('OPT','TOT.EN.') ] = obj_calc_dict['TOT.EN.']
            opt_res_dict[ ('OPT','INT.EN.') ] = obj_calc_dict['INT.EN.'] 
            opt_res_dict[ ('OPT','MULL.CH.CAT.') ] = cat_charges
            opt_res_dict[ ('OPT','MULL.CH.ANI.') ] = ani_charges
            
    else:
       zmat_dict, cart_file = make_coordinates( T, P, R, anion, cation )
       obj_calc.write_input_file_ZMAT( zmat_dict )
       slurm_obj = SLURM.SLURM( obj_calc.run_dir, 'GAMESS', job_name = obj_inp, job_queue = 'nodeshiq' )
       slurm_obj.write_batch()
       slurm_obj.submit_batch()

       zmat_run_dir_list = obj_calc.run_dir.split('/')
       zmat_inp_file_list = obj_calc.inp_file.split('/')
       zmat_run_dir_list[-1] =  'CART_' + zmat_run_dir_list[-1]
       cart_run_dir  = '/'.join( zmat_run_dir_list )
       cart_inp_file = 'cart_' + zmat_inp_file_list[-1]
       
       obj_calc.write_input_file_CART( cart_file )
       slurm_obj = SLURM.SLURM( cart_run_dir , 'GAMESS', job_name = cart_inp_file, job_queue = 'nodeshiq' )
       slurm_obj.write_batch()
       slurm_obj.submit_batch()


if __name__ == '__main__':
  main()

