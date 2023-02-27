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
mono_dict_file = 'monomers.json'
with open(mono_dict_file,'r') as json_file:  
    mono_dict = json.load(json_file)

run_ids, run_job_labels = running_jobs()

template = open( 'emim_bf4.inp', 'r' ).readlines()

def write_new_and_run( root_dir, inp_label, runtyp, fun, gb, template, rad, queue ):
    calc_type = 'DFTTYP'
    os.makedirs( root_dir, exist_ok=True )
    os.chdir( root_dir )
    inp_file = '{}.inp'.format(inp_label)
    with open( inp_file, 'w+' ) as nf:
      for line in template:
        if   line == '  RUNTYP=ENERGY\n':
             new_line = '  RUNTYP={}\n'.format(runtyp)
             #if calc_type == 'DFTTYP':
             #   new_line += '  NPRINT=-5\n'
        elif line == '  DFTTYP=B3LYP\n':
             if calc_type == 'DFTTYP':
                new_line = '  DFTTYP={}\n'.format(fun)
                if fun == 'B2PLYP':
                   new_line += '  NUMGRD=.T.\n'
             elif calc_type == 'MP2':
                new_line = '  MPLEVL=2\n'
             if gb in ispher_1_list:
                new_line += '  ISPHER=1\n'
        elif line == '  GBASIS=N311\n':
             new_line = '  GBASIS={}\n'.format(gb)
             polar_key = [ k for (k,v) in polar_dict.items() if gb in v  ]
             if polar_key == []:
                pass
             else:
                new_line += '  POLAR={}\n'.format(polar_key[0])
        elif line == 'bn20      7.0000000\n':
             new_line = 'bn20      {}\n'.format(rad)
        else:
             new_line = line
        nf.write( new_line )
    shutil.copy( '/home/mdi0316/scripts/submit_gamess.sh' , 'submit_gamess.sh' )
    sp.call( 'sbatch -p {} -J {} submit_gamess.sh {}'.format(queue, inp_file, inp_file),shell=True)

coords_columns = [ ('Coordinates', ii) for ii in [ 'Radius' ]]
opt_columns = [ ('OPT', ii) for ii in [ 'TOT.EN.', 'INT.EN.', 'MULL.CH.CAT.','MULL.CH.ANI.']] 
mp2_columns = [ ('MP2', ii) for ii in [ 'SCF.EN.', 'MP2.EN.', 'SCF.INT.EN.', 'MP2.INT.EN.', 'MULL.CH.CAT.','MULL.CH.ANI.' ]]
eda_columns = [ ('EDA', ii) for ii in [ 'STATUS', 'EXEC.', 'SCF', 'ES.','EX.','REP.','POL.','INT.EN.' ]]
ene_columns = [ ('ENE', ii) for ii in [ 'STATUS', 'EXEC.', 'SCF', 'TOT.EN.', 'INT.EN.' ]]

def main():
    #for gb in gbasis_list:
    for gb in ['APCseg-1']:
      for fun in functionals_list:
      #for fun in ['B3LYP', 'PBE0']:
          if fun == 'B2PLYP':
             pass
          else:
            try:
               print( '>>> {} {} <<<'.format(gb, fun) )
               data_dir = os.path.join( scan_dir, 'DATA')
               gb_fun_label = '{}_{}'.format(gb,fun)
               opt_label = 'emim_bf4_DFT_{}_OPT'.format(gb_fun_label)
               opt_csv = os.path.join( '{}/opt_{}.csv'.format(data_dir,gb_fun_label) )
               mp2_csv = os.path.join( '{}/mp2_{}.csv'.format(data_dir,gb_fun_label) )
    
               opt_df = pd.DataFrame( columns = pd.MultiIndex.from_tuples(coords_columns + opt_columns), dtype=object )
               mp2_df = pd.DataFrame( columns = pd.MultiIndex.from_tuples(coords_columns + mp2_columns), dtype=object )


               emim_opt_dict = mono_dict['EMIM']['OUT'][gb][fun]['DFT']['OPT']['FINAL']
               bf4_opt_dict  = mono_dict[ 'BF4']['OUT'][gb][fun]['DFT']['OPT']['FINAL']
               emim_mp2_dict = mono_dict['EMIM']['OUT'][gb][fun]['MP2']['ENE']
               bf4_mp2_dict  = mono_dict[ 'BF4']['OUT'][gb][fun]['MP2']['ENE']
               #try:
               error_list = [ vv for (kk,vv) in list(emim_opt_dict.items()) + list(emim_mp2_dict.items()) +
                                                list( bf4_opt_dict.items()) + list( bf4_mp2_dict.items())
                                                if 'ERR' == kk ]
               if len( error_list ) > 0:
                  print( 'emim', [ vv for (kk,vv) in list(emim_opt_dict.items()) + list(emim_mp2_dict.items()) if kk == 'ERR' ] )
                  print(  'bf4', [ vv for (kk,vv) in list( bf4_opt_dict.items()) + list( bf4_mp2_dict.items()) if kk == 'ERR' ] )
                         
                  zero_opt_etot = 0
                  zero_mp2_etot = 0
                  zero_mp2_emp2 = 0
               else: 
                  emim_opt_etot = float(emim_opt_dict['TOT.EN.'])
                  bf4_opt_etot  = float( bf4_opt_dict['TOT.EN.'])
                  emim_mp2_etot = float(emim_mp2_dict['TOT.EN.'])
                  bf4_mp2_etot  = float( bf4_mp2_dict['TOT.EN.'])
                  emim_mp2_emp2 = float(emim_mp2_dict['MP2.EN.'])
                  bf4_mp2_emp2  = float( bf4_mp2_dict['MP2.EN.'])
                  zero_opt_etot = emim_opt_etot + bf4_opt_etot
                  zero_mp2_etot = emim_mp2_etot + bf4_mp2_etot
                  zero_mp2_emp2 = emim_mp2_emp2 + bf4_mp2_emp2
                  #except(KeyError):
               
           
               queue = 'nodesloq'
               radius_list = [ '10.0' , '5.0', '5.5', '7.0', '15.0']
               radius_list = [ '2.3', '2.5', '2.8', '3.0', '3.3', '3.5', '3.8', '4.0', '4.5', '5.0', '5.5', '7.0', '10.0', '15.0', '20.0']
               if gb in ['APCseg-1', 'N311']:
                  if fun in ['B3LYP']:
                     radius_list = complete_R_LIST
                     queue = 'nodeshiq'
               for rad in radius_list: 
                 #try:
                       opt_res_dict = { ('Coordinates','Radius') : float(rad) }
                       mp2_res_dict = { ('Coordinates','Radius') : float(rad) }
                       os.chdir( scan_dir )
                       opt_lab = '{}_{}'.format(opt_label, rad)
                       root_dir = os.path.join( scan_dir, 'T_90','P_90','R_{}'.format(rad))
                       opt_calc = GAMESS.GAMESS_calculation( opt_lab,
                                                             root_dir, zero_en = zero_opt_etot, 
                                                             natoms = 24, runtyp = 'OPTIMIZE', 
                                                             post_scf = 'DFTTYP', basis = gb, functional = fun )
                       if os.path.isdir( opt_calc.run_dir ):
                          if opt_calc.out_file!= None:
                             opt_exec = opt_calc.get_execution( run_job_labels )
                             if opt_exec == 'NORMALLY':
                               #print( 'Reading: {}'.format(opt_calc.out_file) )
                               opt_calc_dict = opt_calc.get_out_dict()['DFT']['OPT']['FINAL']
                               opt_calc_geom = opt_calc_dict['GEOM.']
                               opt_calc_scf  = opt_calc_dict['SCF']
                               if opt_calc_geom == 'LOCATED' and opt_calc_scf == 'CONVERGED':
                                  all_charges = opt_calc_dict['MULL.CHARGES'] 
                                  cat_charges = sum([ v['charge'] for (k,v) in all_charges.items() if k < 19 ] )
                                  ani_charges = sum([ v['charge'] for (k,v) in all_charges.items() if k >= 19 ] )
                                  opt_res_dict[ ('OPT','TOT.EN.') ] = opt_calc_dict['TOT.EN.']
                                  opt_res_dict[ ('OPT','INT.EN.') ] = opt_calc_dict['INT.EN.'] 
                                  opt_res_dict[ ('OPT','MULL.CH.CAT.') ] = cat_charges
                                  opt_res_dict[ ('OPT','MULL.CH.ANI.') ] = ani_charges
                                  
                                  ## Moller-Plesset analysis
                                  mp2_label = 'emim_bf4_MP2_{}_{}_ENE'.format(gb,fun)
                                  mp2_calc = GAMESS.GAMESS_calculation( mp2_label,
                                                                        root_dir, zero_en = zero_mp2_etot,
                                                                        natoms = 24, runtyp='ENERGY', 
                                                                        post_scf='MP2', basis=gb, functional=fun )

                                  if os.path.isdir( mp2_calc.run_dir ):
                                     if mp2_calc.out_file != None:
                                        mp2_exec = mp2_calc.get_execution( run_job_labels )
                                        if mp2_exec == 'NORMALLY':
                                           mp2_calc_dict = mp2_calc.get_out_dict()['MP2']['ENE']
                                           mp2_charges = mp2_calc_dict['MULL.CHARGES'] 
                                           cat_charges = sum([ v['charge'] for (k,v) in mp2_charges.items() if k < 19 ] )
                                           ani_charges = sum([ v['charge'] for (k,v) in mp2_charges.items() if k >= 19 ] )
                                           #print( 'Reading: {}'.format(mp2_calc.out_file) )
                                           mp2_etot = mp2_calc_dict['TOT.EN.']
                                           mp2_emp2 = mp2_calc_dict['MP2.EN.']
                                           mp2_eint = mp2_etot - zero_mp2_etot
                                           mp2_Dmp2 = mp2_emp2 - zero_mp2_emp2
                                           mp2_res_dict[ ('MP2','SCF.EN.') ] = mp2_etot 
                                           mp2_res_dict[ ('MP2','MP2.EN.') ] = mp2_emp2
                                           mp2_res_dict[ ('MP2','SCF.INT.EN.') ] = mp2_eint
                                           mp2_res_dict[ ('MP2','MP2.INT.EN.') ] = mp2_Dmp2
                                           mp2_res_dict[ ('MP2','MULL.CH.CAT.') ] = cat_charges
                                           mp2_res_dict[ ('MP2','MULL.CH.ANI.') ] = ani_charges
                                  
                                  else:
                                     dft_opt_zmat = opt_calc_dict['ZMAT'] 
                                     mp2_calc.write_input_file_ZMAT( dft_opt_zmat, msg='dfttyp.optimized.zmat' )
                                     os.chdir( mp2_calc.run_dir )
                                     slurm_obj = SLURM.SLURM( mp2_calc.run_dir, 'GAMESS', job_name = mp2_calc.inp_name, 
                                                              job_queue = 'nodeshiq' )
                                     slurm_obj.write_batch()
                                     slurm_obj.submit_batch()

                                  opt_df = opt_df.append( pd.Series( opt_res_dict ), ignore_index=True )
                                  mp2_df = mp2_df.append( pd.Series( mp2_res_dict ), ignore_index=True )
                               else:
                                  print( 'NOT Converged', opt_calc_geom, opt_calc_scf, opt_calc.run_dir )
                             else:
                               print( 'ABNORMALLY ', opt_calc.run_dir, opt_calc.read_error() )
                       else:
                          write_new_and_run( opt_calc.run_dir, opt_lab, 'OPTIMIZE', fun, gb, template, rad, queue = queue )

               opt_df.to_csv( opt_csv )
               mp2_df.to_csv( mp2_csv )
            except(IndexError,NameError,KeyError,TypeError):
               pass
    
if __name__ == '__main__':
  main()
