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
from Functions import print_tab, running_jobs, find_last_log
from IONIC_LIQUID import mono_dict
from GAMESS import functionals_list,gbasis_list, ispher_1_list, polar_dict

import numpy.linalg as LA
import scipy.constants as const
Ha2eV =  const.value('hartree-electron volt relationship') #27.211
Ang2Bohr = 1.8897259886

run_ids, run_job_labels = running_jobs()

runtyp = 'OPTIMIZE'
 
from bs_conv_dimers import work_dir, mono_dir

def write_mono_and_run( root_dir, inp_label, calc_typ, run_typ, funct, basis, template, radius, queue ):
    os.makedirs( root_dir, exist_ok=True )
    os.chdir( root_dir )
    inp_file = '{}.inp'.format(inp_label)
    with open( inp_file, 'w+' ) as nf:
      for line in template:
        if   line == '  RUNTYP=ENERGY\n':
             new_line = '  RUNTYP={}\n'.format(run_typ)
        elif line == '  DFTTYP=B3LYP\n':
             if calc_type == 'DFTTYP':
                new_line = '  DFTTYP={}\n'.format(funct)
                #if fun == 'B2PLYP':
                #   new_line += '  NUMGRD=.T.\n'
             elif calc_type == 'MP2':
                new_line = '  MPLEVL=2\n'
             elif calc_type == 'CCSD(T)':
                new_line = '  CCTYP=CCSD(T)\n'
             if basis in ispher_1_list:
                new_line += '  ISPHER=1\n'
        elif line == '  GBASIS=N311\n':
             new_line = '  GBASIS={}\n'.format(basis)
             polar_key = [ k for (k,v) in polar_dict.items() if basis in v  ]
             if polar_key == []:
                pass
             else:
                new_line += '  POLAR={}\n'.format(polar_key[0])
        else:
             new_line = line
        nf.write( new_line )
    shutil.copy( '/home/mdi0316/scripts/submit_gamess.sh' , 'submit_gamess.sh' )
    sp.call( 'sbatch -p {} -J {} submit_gamess.sh {}'.format(queue, inp_file, inp_file),shell=True)

mono_dict = {
   'BF4' : {'nat' :  5, 'charge': 1},
  'EMIM' : {'nat' : 19, 'charge':-1}
}

post_dict = {
   'DFT_opt'     : { 'runtyp' : 'OPTIMIZE', 'post_typ' : 'DFTTYP'  , 'read' : None},
#   'DFT_ene'     : { 'runtyp' : 'ENERGY',   'post_typ' : 'DFTTYP'  , 'read' : 'DFT_opt' },
#   'MP2_opt'     : { 'runtyp' : 'OPTIMIZE', 'post_typ' : 'MP2'     , 'read' : 'DFT_opt' },
   'MP2_ene'     : { 'runtyp' : 'ENERGY',   'post_typ' : 'MP2'     , 'read' : 'DFT_opt' },
#   'CCSD(T)_opt' : { 'runtyp' : 'OPTIMIZE', 'post_typ' : 'CCSDT'   , 'read' : 'DFT_opt' },
   'CCSD(T)_ene' : { 'runtyp' : 'ENERGY',   'post_typ' : 'CCSDT'   , 'read' : 'DFT_opt' }
}


def get_output(gms_obj):
    post_lab = gms_obj.post_scf
    if gms_obj.post_scf == 'DFTTYP':
       post_lab = 'DFT'
    run_lab  = gms_obj.runtyp[:3] 
    gms_scf  = False
    gms_geom = False
    gms_zmat = False
    if os.path.isdir( gms_obj.run_dir ): 
       gms_exec = gms_obj.get_execution( run_job_labels )
       if gms_exec == 'Running':
          pass
       elif gms_exec == 'NORMALLY':
          if gms_obj.runtyp == 'OPTIMIZE':
             gms_dict = gms_obj.get_out_dict()[post_lab][run_lab]['FINAL']
             gms_geom = gms_dict['GEOM.']
             gms_zmat = gms_dict['ZMAT']
          else:
             gms_dict = gms_obj.get_out_dict()[post_lab][run_lab] 
          gms_scf = gms_dict['SCF']
       else:
          print( "#W@",  gms_obj.run_dir, gms_obj.read_error() )
          gms.
    else:
       gms_exec = 'NO.FOLDER'
    return( gms_exec, gms_scf, gms_geom, gms_zmat )


def main():
    for mono_label, mono_v in mono_dict.items():
       mono_natom  = mono_v['nat']
       mono_charge = mono_v['charge']
       template = open( os.path.join( work_dir, '{}.inp'.format(mono_label.lower())) ,'r' ).readlines()
       root_dir = os.path.join( mono_dir, mono_label )
       for basis in gbasis_list:
         for funct in functionals_list:

             os.chdir(root_dir)
             print_tab( 0, '{} {} {}'.format(mono_label, basis, funct) )

             for k,v in post_dict.items():
                 print_tab(1, k)
                 tmp_run  = v['runtyp']
                 tmp_post = v['post_typ']
                 tmp_read = v['read']

                 tmp_obj = GAMESS.GAMESS_calculation( mono_label.lower(), root_dir=root_dir, 
                                                      natoms = mono_natom, icharge = mono_charge, 
                                                      runtyp = tmp_run, post_scf = tmp_post, 
                                                      basis = basis, functional = funct) #, verbose=True  )
                 tmp_exec, tmp_scf, tmp_geom, tmp_zmat = get_output( tmp_obj )
                 print_tab( 2, tmp_exec )
        
                 if tmp_exec == 'NO.FOLDER':

                   if tmp_read == None:
                      print( write_from_template, dft_opt )
                   else:
                      read_run  = post_dict[tmp_read]['runtyp']
                      read_post = post_dict[tmp_read]['post_typ']
                      read_obj = GAMESS.GAMESS_calculation( mono_label.lower(), root_dir=root_dir, 
                                                            natoms = mono_natom, icharge = mono_charge, 
                                                            runtyp = read_run, post_scf = read_post, 
                                                            basis = basis, functional = funct )#, verbose=True  )
                      read_exec, read_scf, read_geom, read_zmat = get_output( read_obj )
                      print_tab( 2, [ 'Reading from: {}'.format(read_obj.run_dir) , 'Read status = {}'.format( read_exec )] )
 
                      if [read_exec, read_scf, read_geom] == ['NORMALLY','CONVERGED','LOCATED']:
                         tmp_obj.write_input_file_ZMAT( read_zmat )
                         slurm_obj = SLURM.SLURM( tmp_obj.run_dir, 'GAMESS', job_name=tmp_obj.inp_name, job_queue = 'nodeshiq' )
                         slurm_obj.write_batch()
                         slurm_obj.submit_batch()

       
             

if __name__ == '__main__':
   main()
