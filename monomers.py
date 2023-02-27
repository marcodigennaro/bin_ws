#!/home/mdi0316/anaconda3/bin/python

### common input start
import os, sys, re
import numpy as np
import shutil
import pandas as pd

import getpass
user = getpass.getuser()

scripts_dir = '/home/{}/FUNCTIONS'.format(user)
classes_dir = '/home/{}/CLASSES'.format(user)
zmat_converter_dir = '/home/{}/CLASSES/zmatrix-master'.format(user)

sys.path.insert(0, scripts_dir)
sys.path.insert(0, classes_dir)
sys.path.insert(0, zmat_converter_dir)

import json
import subprocess as sp

import GAMESS
import SLURM

from Functions import print_tab, running_jobs, running_label
from GAMESS import full_functionals_list, full_gbasis_list
from dimers import print_converged_zmat, print_internucl_dist, print_dft_results, print_mp2_results, read_object, read_zmat, print_ccsdt_results
from make_scan_list import *

work_dir = '/data/{}/WORK'.format(user)

mono_runs_dir = os.path.join( work_dir, 'MONOMERS/RUNS' )
mono_csv_dir = os.path.join( work_dir, 'MONOMERS/CSV' )

os.makedirs( mono_runs_dir, exist_ok = True )
os.makedirs( mono_csv_dir, exist_ok = True )

temp_dir = '/home/{}/Inputfiles/GAMESS/MONOMERS/AVOGADRO/'.format(user)

mono_json = os.path.join( work_dir, 'monomers_{}.json'.format(user) )
print( 'Reading: {}'.format( mono_json ))
with open(mono_json,'r') as json_file:
   mono_dict = json.load(json_file)

verbose=False
verbose=True

def locate_mono_line(df, mono_label, basis, funct ):

    line_complete = False
    if df.empty:
       mono_line = pd.DataFrame()
       mono_idx = float('nan')
       print_tab( 3, 'new DataFrame' )
    else:
       mono_line = df[ (df['LABEL'] == mono_label) & (df['BASIS'] == basis) & (df['FUNCT'] == funct) ]
       if mono_line.empty:
          mono_idx = float('nan')
          print_tab( 3, '{} {} missing'.format(basis, funct) )
       else:
          mono_idx  = mono_line.index.values[0]
          null_columns = mono_line.columns[ mono_line.isnull().any() ]
          if null_columns.empty:
             if all( [ x in mono_line.columns for x in ['OPT.DFT.TOT.EN.', 'DFT.TOT.EN.', 'CC.CCSD[T]', 'MP2.TOT.EN.' ] ] ):
                print_tab( 3, 'idx. = {}: complete'.format(mono_idx) )
                line_complete = True
             else:
                print_tab( 3, 'idx. = {}: incomplete'.format(mono_idx) )
          else:
             print_tab( 3, 'idx. = {}: empty columns found'.format(mono_idx) )

    return line_complete, mono_idx

def template_to_dict( template_file ):
    with open( template_file, 'r' ) as tmpl_file:
       tmpl_lines = tmpl_file.readlines()

    tmpl_dict = {}
    for tmpl_l in tmpl_lines:
      if tmpl_l.split() == []:
         pass
      elif '=' in tmpl_l:
         pass
      else:
         print( tmpl_l.split() )

def main():

    global DIMER
    DIMER = False
    
    # read from monomers_mdi0316.json
    # run /home/mdi0316/bin/initialize_mono_dict.py

    mono_label = sys.argv[1]   
    mono_v = mono_dict[mono_label]
    mono_natom  = mono_v['nat']
    mono_charge = mono_v['charge']
    mono_mult   = mono_v['mult']
    mono_scftyp = mono_v['scftyp']

    mono_zmat_template = os.path.join( temp_dir, '{}.inp'.format(mono_label.lower())) 
    mono_csv_file = os.path.join( mono_csv_dir, '{}.csv'.format(mono_label) )
    
    print_tab( 0, '>>>> {} <<<<'.format(mono_label) )
    print_tab( 1, 'Template: {}'.format(mono_zmat_template ))
    print_tab( 1, 'csv file: {}'.format(mono_csv_file ))
    
    full_gbasis_list = ['N311']
    full_functionals_list = ['B3LYP']
    if mono_label in ['EMIM', 'BF4', 'PF6']:
       full_gbasis_list = GAMESS.full_gbasis_list
       full_functionals_list = GAMESS.full_functionals_list

       #full_gbasis_list, full_functionals_list = extended_basis_list, extended_funct_list
       #full_gbasis_list, full_functionals_list = ['N311', 'STO', 'APCseg-1', 'CCQ'], ['B3LYP', 'M11', 'PBE0', 'wB97x-D' ]

    if os.path.exists( mono_csv_file ):
       mono_df = pd.read_csv( mono_csv_file, index_col = 0 )
    else:
       mono_df = pd.DataFrame()

    for basis in full_gbasis_list:
        print_tab( 1, '=== {} ==='.format(basis) )
        for funct in full_functionals_list:
            print_tab( 2, '--- {} ---'.format(funct) )
        
            line_complete, mono_tmp_idx = locate_mono_line( mono_df, mono_label, basis, funct )

            if line_complete == False: #mono_tmp_line.empty:

               root_dir = os.path.join( mono_runs_dir, mono_label, basis, funct )
               os.makedirs( root_dir, exist_ok=True )
               os.chdir(root_dir)
               
               ## MAKE OPT
               mono_tmp_dict = { 'LABEL' : mono_label, 'BASIS' : basis, 'FUNCT' : funct }
               opt_label = '{}_{}_{}'.format( mono_label.lower(), basis, funct )
               opt_obj = GAMESS.GAMESS( inp_label = opt_label, root_dir = root_dir, natoms = mono_natom, icharge = mono_charge, 
                  	       	        run_type = 'OPTIMIZE', post_scf = 'DFTTYP', mult  = mono_mult, scftyp = mono_scftyp,
                  		        basis = basis, functional = funct )
               opt_exec, opt_exec_err, opt_gms_err, opt_scf, opt_geom, opt_time = read_object( opt_obj, dimer = False, 
                                                                    read_template = mono_zmat_template, read_msg = mono_label )
               mono_tmp_dict['OPT.EXEC.']       = opt_exec
               mono_tmp_dict['OPT.EXEC.ERR.']   = opt_exec_err
               mono_tmp_dict['OPT.GMS.ERR.']    = opt_gms_err
               mono_tmp_dict['OPT.SCF']         = opt_scf
               mono_tmp_dict['OPT.GEOM']        = opt_geom
               mono_tmp_dict['OPT.TIME']        = opt_time

               if opt_exec == 'Running':
                  print( opt_exec )
               elif [ opt_exec, opt_exec_err, opt_gms_err, opt_scf ] == [ 'TERMINATED.NORMALLY', False, False, 'CONVERGED' ]:
                  opt_out_dict = opt_obj.get_job_results()[1]
                  mono_tmp_dict['OPT.DFT.TOT.EN.'] = opt_out_dict['TOT.EN.']
                  opt_dict = read_zmat( opt_obj.zmat_file )
 
                  for post_lab in [ 'DFTTYP', 'MP2', 'CCSDT' ]:
                      pst_label = '{}_{}_{}'.format( mono_label.lower(), basis, funct )
                      post_obj = GAMESS.GAMESS( inp_label = pst_label, root_dir = root_dir, 
                                                natoms = mono_natom, icharge = mono_charge, 
                                                run_type = 'ENERGY', post_scf = post_lab,
                                                mult  = mono_mult, scftyp = mono_scftyp,
                                                basis = basis, functional = funct )
                     
                      post_exec, post_exec_err, post_gms_err, post_scf, post_geom, post_time = \
                                       read_object( post_obj, dimer = False, 
                                       read_dict = opt_dict, read_msg = '{}.{}'.format(mono_label, post_lab)  )
                      if [ post_exec, post_exec_err, post_gms_err ] == [ 'TERMINATED.NORMALLY', False, False ]:
                         post_out_dict = post_obj.get_job_results()[1]
                         if post_lab == 'DFTTYP':
                            dft_line = print_dft_results( post_obj, post_out_dict, dimer=False, distance=False, zero_dft_ener=False ) 
                            mono_tmp_dict[ 'DFT.TOT.EN.' ]    = dft_line['TOT.EN.'].values[0]
                            mono_tmp_dict[ 'DFT.BANDGAP']     = dft_line['BANDGAP'].values[0]
                            mono_tmp_dict[ 'DFT.BASIS.DIM.' ] = dft_line['BASIS.DIM.'].values[0]
                         elif post_lab == 'MP2':
                            mp2_line = print_mp2_results( post_obj, post_out_dict, dimer=False, distance=False, zero_mp2_ener=False ) 
                            mono_tmp_dict[ 'MP2.TOT.EN.' ] = mp2_line['EN.MP2'].values[0]
                            mono_tmp_dict[ 'MP2.BANDGAP' ] = mp2_line['BANDGAP'].values[0]
                         elif post_lab == 'CCSDT':
                            ccsdt_line = print_ccsdt_results( post_obj, post_out_dict, dimer=False ) 
                            mono_tmp_dict['CC.REF.EN']  = ccsdt_line['REF.EN'].values[0] 
                            mono_tmp_dict['CC.MBPT(2)'] = ccsdt_line['MBPT(2)'].values[0]
                            mono_tmp_dict['CC.CCSD']    = ccsdt_line['CCSD'].values[0] 
                            mono_tmp_dict['CC.CCSD(T)'] = ccsdt_line['CCSD(T)'].values[0]
                            mono_tmp_dict['CC.CCSD[T]'] = ccsdt_line['CCSD[T]'].values[0]

                  if all( [ x in mono_tmp_dict.keys() for x in ['OPT.DFT.TOT.EN.', 'DFT.TOT.EN.'] ] ):
                     if abs(mono_tmp_dict['OPT.DFT.TOT.EN.'] - mono_tmp_dict['DFT.TOT.EN.'] ) > .001:
                        print( mono_tmp_dict['OPT.DFT.TOT.EN.'] )
                        print( mono_tmp_dict['DFT.TOT.EN.'] ) 
                        print( 'hhhhhhhhhhhhh', opt_obj.run_dir )
                        #raise NameError( "OPT and DFT energies are different" ) 

               else: 
                  print( 'unknown status for result: ', opt_obj.run_dir )
                  print(  opt_exec, opt_exec_err, opt_gms_err, opt_scf )
               #new_mono_line = pd.DataFrame( [mono_tmp_dict] )
               #null_columns = mono_line.columns[ new_mono_line.isnull().any() ]
               #if null_columns.empty:
               #   mono_df.loc[mono_idx] = new_mono_line
               if np.isnan(mono_tmp_idx):
                  print_tab( 3, 'New DataFrame/New Line' )
                  mono_df = mono_df.append( [mono_tmp_dict], ignore_index = True )
                  mono_df.to_csv( mono_csv_file )
               else:
                  print_tab( 3, 'replacing line' )
                  mono_df.at[int(mono_tmp_idx)] = pd.Series( mono_tmp_dict )

    print_tab( 0, 'printing {}'.format( mono_csv_file  ))
    mono_df.to_csv( mono_csv_file )

if __name__ == '__main__':
   main()
