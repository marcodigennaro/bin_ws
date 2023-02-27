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
from Functions import print_tab, running_jobs, find_last_log, calculate_center_of_mass, calculate_center_of_charge, Coulomb_Energy
from IONIC_LIQUID import R_LIST, T_LIST, P_LIST
import numpy.linalg as LA
import scipy.constants as const
Ha2eV =  const.value('hartree-electron volt relationship') #27.211
Ang2Bohr = 1.8897259886

from IONIC_LIQUID import mono_dict

from pyionic_liquid import inp_dir, work_dir, mono_dir, get_gms_info, read_df_line, write_and_run_new, resubmit_jobs, POST_SCF, BASIS, FUNCTIONAL

## define calculation labels and dictionaries ##
global IL_LABEL
IL_LABEL  = sys.argv[1].split('/')[0]

## define IL from CLASS 
IL_dir  = os.path.join( work_dir, IL_LABEL )
cat_label, ani_label = IL_LABEL.split('_')

cat_atnum, cat_type, cat_inp, cat_composition, cat_log_dict = mono_dict[cat_label].values()
ani_atnum, ani_type, ani_inp, ani_composition, ani_log_dict = mono_dict[ani_label].values()
cat_log_file = cat_log_dict[POST_SCF][FUNCTIONAL][BASIS]
ani_log_file = ani_log_dict[POST_SCF][FUNCTIONAL][BASIS]

global N_ATOMS 
N_ATOMS = cat_atnum + ani_atnum

cat_dir = os.path.join( mono_dir, cat_label, POST_SCF, BASIS, FUNCTIONAL )
ani_dir = os.path.join( mono_dir, ani_label, POST_SCF, BASIS, FUNCTIONAL )

cation = GAMESS.GAMESS_calculation( cat_dir, cat_inp, cat_log_file,
                                    runtyp = 'OPTIMIZE', post_scf = POST_SCF, basis = BASIS, functional = FUNCTIONAL,
                                    natoms = cat_atnum )

anion  = GAMESS.GAMESS_calculation( ani_dir, ani_inp, ani_log_file, 
                                    runtyp = 'OPTIMIZE', post_scf = POST_SCF, basis = BASIS, functional = FUNCTIONAL,
                                    natoms = ani_atnum )

global ZERO_ENERGY
cat_ener = cation.get_out_dict()['OPT']['TOT.EN.']
ani_ener = anion.get_out_dict()['OPT']['TOT.EN.']
cat_atoms_list = cation.get_atoms_list()
ani_atoms_list = anion.get_atoms_list( shift = cat_atnum )
ZERO_ENERGY = cat_ener + ani_ener  # -20924.916677969093 eV

IL = IONIC_LIQUID.IONIC_LIQUID( IL_LABEL, IL_dir , zero_energy=ZERO_ENERGY, cation_label=cat_label, anion_label=ani_label,
                                post_scf = POST_SCF, basis = BASIS, functional = FUNCTIONAL )
for fold in [ IL.scan_dir, IL.csv_dir ]:
    os.makedirs(fold, exist_ok = True)

input_csv  = IL.inp_scan_csv
result_csv = IL.opt_scan_csv
scan_type = 'SCAN'
runtyp  = 'OPTIMIZE'

if len( sys.argv ) == 3:
   scan_type = sys.argv[2]
   if scan_type == 'SCAN':
      runtyp  = sys.argv[3]
      if runtyp == 'ENERGY':
         result_csv = IL.ene_scan_csv
      elif runtyp == 'EDA':
         result_csv = IL.eda_scan_csv
   elif scan_type == 'EQUILIBRIUM':
      input_csv  = IL.equil_inp_csv
      result_csv = IL.equil_out_csv
   elif scan_type == 'MIN.SURF.':
      input_csv  = IL.min_surf_inp_csv
      result_csv = IL.min_surf_out_csv
   elif scan_type == 'GPR':
      gpr_idx  = str(sys.argv[3])
      input_csv = os.path.join( IL.csv_dir, 'gpr_{}_inp.csv'.format(gpr_idx) )
      result_csv = IL.opt_gpr_csv
   elif scan_type == 'RANDOM':
      input_csv  = IL.inp_rndm_csv
      result_csv = IL.opt_rndm_csv
   else:
      print( "possible values: 'EQUILIBRIUM', 'MIN.SURF.', 'GPR', 'RANDOM'")
      print( unrecognized_scan_type )

run_ids, run_job_labels = running_jobs()

def resubmit(resumbit_obj, error, close_R_tuple=None):
    gms_inp_name = resumbit_obj.inp_name
    os.chdir( resumbit_obj.run_dir )
    # save failed files 
    skip = False
    if os.path.exists( 'FAILED' ):
       print_tab( 2, ['WARNING: FAILED folder exists', resumbit_obj.run_dir] )
       all_failed = [ int(item) for item in os.listdir('FAILED') ]
       last_one = os.path.join( resumbit_obj.run_dir, 'FAILED', str(max(all_failed)), gms_inp_name )
       if filecmp.cmp( last_one, resumbit_obj.inp_file ):
          print_tab( 2, 'WARNING: no change has been made to the input file' )
          skip = True
    
    if not skip:
       #print_tab( 2 , [ error, gms_inp_name, resumbit_obj.run_dir ] )
       os.makedirs( 'FAILED', exist_ok=True )
       fail_fold = os.path.join( resumbit_obj.run_dir, 'FAILED',  str(len(os.listdir( 'FAILED' )) ) )
       os.makedirs( fail_fold )
       shutil.copy2( resumbit_obj.inp_file, fail_fold )
       shutil.copy2( resumbit_obj.opt_file, fail_fold )
       shutil.copy2( resumbit_obj.dat_file, fail_fold )
       with open( os.path.join(fail_fold, 'fail.txt'), 'w+' ) as fail_txt:
            now = datetime.datetime.now()
            fail_txt.write( 'FAILED on {}\nERROR: {}'.format(now, error) )
       
       ## write new files 
       if   error == 'TOO MANY STEPS TAKEN':
            ## new_Zmat from same directory
            new_Zmat = resumbit_obj.get_last_geometry()[2]
       elif error in ['SCF.Failed', '-ABNORMALLY-']:
            ## new_Zmat from closest converged R directory
            [tmp_R, tmp_T, tmp_P], tmp_closest_R = close_R_tuple
            new_obj = get_gms_info( tmp_r = tmp_closest_R, tmp_t = tmp_T, tmp_p = tmp_P )
            new_Zmat = new_obj.get_equilibrium_geometry()[2]
            print_tab( 2, [ 'Failed RTP    = {}, {}, {}'.format( tmp_R, tmp_T, tmp_P ),
                            'CONVERGED RTP = {}, {}, {}'.format( tmp_closest_R, tmp_T, tmp_P ) ] )
            ## set RTP back to tmp_values
            new_Zmat[51][5] = float(tmp_R)*Ang2Bohr
            new_Zmat[51][6] = tmp_R
            new_Zmat[52][6] = tmp_T
            new_Zmat[53][6] = tmp_P
       else:
            print( unknown_error )   
       ## resubmit
       resumbit_obj.write_input_from_Zmat(new_Zmat)
       sp.call( 'sbatch -J {} submit_gamess.sh {}'.format(gms_inp_name, gms_inp_name), shell=True )


## !!ALL energies in EV!!
def main():

    global R
    global T
    global P

    ## READ previously written results
    print_tab( 0, 'input_csv       {}'.format(input_csv) )
    print_tab( 0, 'result_csv      {}'.format(result_csv) )
    print_tab( 0, 'scan_type       {}'.format(scan_type) )
    if not os.path.exists( input_csv ):
       IL.write_inp_csv( scan_type = scan_type )

    input_df = pd.read_csv( input_csv, header = [0,1], dtype=object ) #index_col=0
    
    print( input_df )
    for index, row in input_df.iterrows():
        R, T, P  = str(row[('Coordinates','Radius')]), str(row[('Coordinates','Theta')]), str(row[('Coordinates','Phi')])
        run_obj = get_gms_info( R, T, P, runtyp , scan_type ) 
        run_exec = run_obj.get_execution(run_job_labels)
        print( T,P,R, run_exec )
        if run_exec == 'Running':
           run_stat = 'Running'
           print_tab(3 , run_stat )
        elif run_exec == 'Missing':
           run_exec = resubmit_jobs( run_obj.run_dir )
        elif run_exec in ['unknown.exec.', '-ABNORMALLY-']:
           run_stat = run_obj.read_error() # 'FAIL.'
        elif run_exec == 'NORMALLY':
           run_stat = 'unknown.stat.'
           gamess_out_dict = run_obj.get_out_dict()['OPT']
           print( gamess_out_dict['SCF'] )
           print( gamess_out_dict['GEOM.'] )
           run_obj.read_warning()
        

if __name__ == '__main__':
  main()
