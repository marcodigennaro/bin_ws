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

dft_dict = { 'N311' : { 'gbasis' : 'N311' , 'ngauss' : 6 } }

"""
Theta = BEND
Phi   = TORSION (or dihedral)
bn20  = Radius [~7]
bcn20 = Theta  [0,180]
dih20 = Phi    [0,360]
"""

inp_dir  = '/home/mdi0316/Inputfiles/GAMESS'
work_dir = '/data/mdi0316/WORK'
mono_dir = os.path.join( work_dir, 'MONOMERS' )

## define calculation labels and dictionaries ##
global IL_LABEL
IL_LABEL  = sys.argv[1].split('/')[0]
global POST_SCF
global BASIS
global FUNCTIONAL

POST_SCF = 'DFTTYP'
#POST_SCF = 'NONE'

#BASIS = 'APCseg-1'
BASIS = 'N311'
#BASIS = 'STO'

FUNCTIONAL = 'B3LYP'

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
#try:
cat_ener = cation.get_out_dict()['OPT']['TOT.EN.']
ani_ener = anion.get_out_dict()['OPT']['TOT.EN.']
cat_atoms_list = cation.get_atoms_list()
ani_atoms_list = anion.get_atoms_list( shift = cat_atnum )
ZERO_ENERGY = cat_ener + ani_ener  # -20924.916677969093 eV
#except(FileNotFoundError):
#  warnings.warn( 'Cation/Anion not completed' )
#  ZERO_ENERGY = 0

IL = IONIC_LIQUID.IONIC_LIQUID( IL_LABEL, IL_dir , zero_energy=ZERO_ENERGY, cation_label=cat_label, anion_label=ani_label,
                                post_scf = POST_SCF, basis = BASIS, functional = FUNCTIONAL )
for fold in [ IL.scan_dir, IL.csv_dir ]:
    os.makedirs(fold, exist_ok = True)

#print(cation)
#print(anion)
#print(IL)
## default
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

## read running jobs
run_ids, run_job_labels = running_jobs()

#def resubmit(resumbit_obj, error, close_R_tuple=None):
#    gms_inp_name = resumbit_obj.inp_name
#    os.chdir( resumbit_obj.run_dir )
#    # save failed files 
#    skip = False
#    if os.path.exists( 'FAILED' ):
#       print_tab( 2, ['WARNING: FAILED folder exists', resumbit_obj.run_dir] )
#       all_failed = [ int(item) for item in os.listdir('FAILED') ]
#       last_one = os.path.join( resumbit_obj.run_dir, 'FAILED', str(max(all_failed)), gms_inp_name )
#       if filecmp.cmp( last_one, resumbit_obj.inp_file ):
#          print_tab( 2, 'WARNING: no change has been made to the input file' )
#          skip = True
#    
#    if not skip:
#       #print_tab( 2 , [ error, gms_inp_name, resumbit_obj.run_dir ] )
#       os.makedirs( 'FAILED', exist_ok=True )
#       fail_fold = os.path.join( resumbit_obj.run_dir, 'FAILED',  str(len(os.listdir( 'FAILED' )) ) )
#       os.makedirs( fail_fold )
#       shutil.copy2( resumbit_obj.inp_file, fail_fold )
#       shutil.copy2( resumbit_obj.opt_file, fail_fold )
#       shutil.copy2( resumbit_obj.dat_file, fail_fold )
#       with open( os.path.join(fail_fold, 'fail.txt'), 'w+' ) as fail_txt:
#            now = datetime.datetime.now()
#            fail_txt.write( 'FAILED on {}\nERROR: {}'.format(now, error) )
#       
#       ## write new files 
#       if   error == 'TOO MANY STEPS TAKEN':
#            ## new_Zmat from same directory
#            new_Zmat = resumbit_obj.get_last_geometry()[2]
#       elif error in ['SCF.Failed', '-ABNORMALLY-']:
#            ## new_Zmat from closest converged R directory
#            [tmp_R, tmp_T, tmp_P], tmp_closest_R = close_R_tuple
#            new_obj = get_gms_info( tmp_r = tmp_closest_R, tmp_t = tmp_T, tmp_p = tmp_P )
#            new_Zmat = new_obj.get_equilibrium_geometry()[2]
#            print_tab( 2, [ 'Failed RTP    = {}, {}, {}'.format( tmp_R, tmp_T, tmp_P ),
#                            'CONVERGED RTP = {}, {}, {}'.format( tmp_closest_R, tmp_T, tmp_P ) ] )
#            ## set RTP back to tmp_values
#            new_Zmat[51][5] = float(tmp_R)*Ang2Bohr
#            new_Zmat[51][6] = tmp_R
#            new_Zmat[52][6] = tmp_T
#            new_Zmat[53][6] = tmp_P
#       else:
#            print( unknown_error )   
#       ## resubmit
#       resumbit_obj.write_input_from_Zmat(new_Zmat)
#       sp.call( 'sbatch -J {} submit_gamess.sh {}'.format(gms_inp_name, gms_inp_name), shell=True )

def get_gms_info( tmp_r, tmp_t, tmp_p, runtyp, scan_type ):
    tmp_r = tmp_r.strip()
    tmp_t = tmp_t.strip()
    tmp_p = tmp_p.strip()
    if scan_type == 'SCAN':
       run_label = 'SCAN_{}'.format(runtyp[:3])
       head_dir  = IL.scan_dir
    elif scan_type == 'GPR':
       run_label = 'GPR_{}'.format(runtyp[:3])
       tmp_r = '{:4.4f}'.format(float(tmp_r))
       tmp_t = '{:4.4f}'.format(float(tmp_t))
       tmp_p = '{:4.4f}'.format(float(tmp_p))
       head_dir = IL.gpr_dir
    elif scan_type == 'RANDOM':
       tmp_r = '{:4.4f}'.format(float(tmp_r))
       tmp_t = '{:4.4f}'.format(float(tmp_t))
       tmp_p = '{:4.4f}'.format(float(tmp_p))
       run_label = 'RND_{}'.format(runtyp[:3])
       head_dir = IL.rnd_dir
    elif scan_type == 'MIN.SURF.':
       run_label = 'MIN_{}'.format(runtyp[:3])
       head_dir = IL.min_surf_dir
    elif scan_type == 'EQUILIBRIUM':
       run_label = 'EQU_{}'.format(runtyp[:3])
       head_dir = IL.equil_dir
       
    gms_dir = os.path.join( head_dir, POST_SCF, BASIS, FUNCTIONAL, 
                            'T_{}'.format(tmp_t), 'P_{}'.format(tmp_p), 'R_{}'.format(tmp_r) )
    run_dir = os.path.join( gms_dir, runtyp )
    gms_inp = '{}_{}_{}_{}_{}_T_{}_P_{}_R_{}.inp'.format(IL_LABEL,POST_SCF[:3],BASIS,FUNCTIONAL,run_label,tmp_t,tmp_p,tmp_r)
    gms_log = None ## new calc
    if os.path.exists( run_dir ):
       all_logs = [ff for ff in os.listdir(run_dir) if ff.startswith('log') ]
       if len( all_logs ) == 0:
          gms_log = None
       elif len( all_logs ) == 1:
          gms_log = all_logs[-1]
       else:
          print( run_dir )
          print( stop_more_than_one_log )
    
    gms_obj = GAMESS.GAMESS_calculation( gms_dir, gms_inp, gms_log,
                                         zero_en = ZERO_ENERGY, 
                                         runtyp = runtyp, post_scf = POST_SCF, basis = BASIS, functional = FUNCTIONAL,
                                         natoms = N_ATOMS )

    return( gms_obj )

def read_df_line( R_read, T_read, P_read, df_read, run_label, print_opt = False ):
    df_line = df_read.loc[ df_read[('Coordinates', 'Radius')] == str(R_read) ].loc[ 
                           df_read[('Coordinates', 'Theta')]  == str(T_read) ].loc[ 
                           df_read[('Coordinates', 'Phi')]    == str(P_read) ]
    if len(df_line) > 1:
       print( error_more_than_one_line_in_df )
    if df_line.empty:
      df_index = None
      new_line = True
      prev_stat = 'New.Line'
    else:
      df_index = df_line.index.values[0]
      new_line = False
      prev_stat = df_line[(run_label,'STATUS')].values[0]
    if print_opt:
       print_tab( 1, 'R={:5s},T={:3s},P={:3s}: {}:{}'.format(R_read, T_read, P_read, run_label, prev_stat) )
    return( df_line, prev_stat, new_line, df_index )

def write_and_run_new( run_dir, inp_name, runtyp, scan_type = 'SCAN',
                       run_zmat = None, nat_cat = None, nat_ani = None ):

    warnings.warn( 'will skip new files till ATOMATE is correctly installed' )
    ## define new zmatrix from relaxed anion and cation according to R,T,P if not explicitely in input
    #if run_zmat == None:
    #   run_zmat = IL.get_cation_anion_gamess_zmat_dict(radius=R, theta=T, phi=P)

    if scan_type == 'SCAN':
       run_ifreeze = IL.ifreeze_coords
    elif scan_type == 'MIN.SURF.':
       run_ifreeze = ','.join( IL.ifreeze_coords.split(',')[1:] )
    elif scan_type == 'EQUILIBRIUM':
       if not os.path.exists( run_dir ):
          os.makedirs( run_dir )
       run_ifreeze = None
       ## read relaxed anion/cation and create guess coordinates BEGINS
       cat_cart_coords = IL.CATION.cart_dict_to_dat()
       ani_cart_coords = IL.ANION.cart_dict_to_dat()
       IL_guess_cart_coords = os.path.join( run_dir, 'guess_cart_coords.dat' )
       IL_guess_zmat_coords = os.path.join( run_dir, 'guess_zmat_coords.dat' )

       Z = float(R)*np.sin(np.deg2rad(float(T)))
       XY = Z/np.tan(np.deg2rad(float(T)))
       X = XY*np.cos(np.deg2rad(float(P)))
       Y = XY*np.sin(np.deg2rad(float(P)))
       
       if not os.path.exists( IL_guess_cart_coords ):
         with open( IL_guess_cart_coords, 'w+' ) as IL_cc:
              IL_cc.write( '{}\n\n'.format(IL.cation_nat+IL.anion_nat))
              for cat_line in open( cat_cart_coords, 'r' ).readlines()[2:]:
                  IL_cc.write(cat_line)                                   
              for ani_line in open( ani_cart_coords, 'r' ).readlines()[2:]:
                  sym,x,y,z = ani_line.split()
                  IL_cc.write('{:2s}  {:1.8f}  {:1.8f}  {:1.8f}\n'.format( sym, float(x) + X, float(y) + Y, float(z) + Z ) )
       a = Converter()
       a.run_cartesian( input_file = IL_guess_cart_coords, output_file=IL_guess_zmat_coords )

       run_zmat = {}
       run_zmat[0] = { 'elem.':open(cat_cart_coords,'r').readlines()[2].split()[0], 'idx.':'1' } 
       IL_zm_lines = open( IL_guess_zmat_coords, 'r' ).readlines()
       for zm_idx, zm_line in enumerate(IL_zm_lines[3:]):
           if len(zm_line.split()) == 3:
              elem, str_ref, str_val = zm_line.split()
              line_dict = { 'elem.':elem, 'idx.':zm_idx+2, 
                            'STR' : {'ref':str_ref, 'label': 'r{}'.format(zm_idx), 'val':float(str_val)} }
           elif len(zm_line.split()) == 5:
              elem, str_ref, str_val, ben_ref, ben_val = zm_line.split()
              line_dict = {'elem.':elem, 'idx.':zm_idx+2,
			    'STR' : {'ref':str_ref, 'label': 'r{}'.format(zm_idx), 'val':float(str_val)},
                            'BEN' : {'ref':ben_ref, 'label': 'b{}'.format(zm_idx), 'val':float(ben_val)} }
           elif len(zm_line.split()) == 7:
              elem, str_ref, str_val, ben_ref, ben_val, tor_ref, tor_val = zm_line.split()
              line_dict = {'elem.':elem, 'idx.':zm_idx+2,
			    'STR' : {'ref':str_ref, 'label': 'r{}'.format(zm_idx), 'val':float(str_val)},
                            'BEN' : {'ref':ben_ref, 'label': 'b{}'.format(zm_idx), 'val':float(ben_val)},
                            'TOR' : {'ref':tor_ref, 'label': 't{}'.format(zm_idx), 'val':float(tor_val)} }
                                                         
                                                         
           else:
              print(error_unkonwn_coordinate_type)
           run_zmat[zm_idx+1] = line_dict
       ## read relaxed anion/cation and create guess coordinates ENDS
       
    ## define new GAMESS object and run
    new_gms = GAMESS.GAMESS_calculation( run_dir, inp_name, '',
                                         zero_en = ZERO_ENERGY, 
                                         runtyp = runtyp, post_scf = POST_SCF, basis = BASIS, functional = FUNCTIONAL,
                                         natoms = N_ATOMS, ifreeze = run_ifreeze, nat_cat = nat_cat, nat_ani = nat_ani )
    new_gms.write_input_file(run_zmat, msg=inp_name)

    ## submit
    queue = 'nodesloq'
    if runtyp in ['EDA', 'ENERGY'] or scan_type in ['EQUILIBRIUM', 'MIN.SURF.']:
       queue = 'nodeshiq'
    slurm_obj = SLURM.SLURM( new_gms.run_dir, 'GAMESS', job_name = new_gms.inp_name, job_queue = queue )
    slurm_obj.write_batch()
    slurm_obj.submit_batch()
    return( 'Running' )
    return( 'NEW' )

def resubmit_jobs( folder ):
    os.chdir( folder )
    inp_name = [ ff for ff in os.listdir( os.getcwd() ) if ff.endswith('inp') ][0]
    warnings.warn( 'missing log, will resubmit' )
    shutil.copy( '/home/mdi0316/scripts/submit_gamess.sh', './' )
    sp.call( 'sbatch --partition=nodesloq -J {} submit_gamess.sh {}'.format(inp_name, inp_name), shell=True)
    return( 'Running' )

### common input end

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
    
    coords_columns = [ ('Coordinates', ii) for ii in [ 'Radius', 'Theta', 'Phi' ]]
    equil_columns  = [   ('OPT', ii) for ii in [ 'STATUS','EXEC.','SCF','GEOM.','TOT.EN.','INT.EN.','Rlx.R','Rlx.T','Rlx.P' ]]
    min_surf_columns = [ ('OPT', ii) for ii in [ 'STATUS','EXEC.','SCF','GEOM.','TOT.EN.','INT.EN.','Rlx.R' ]]
    opt_columns = [      ('OPT', ii) for ii in [ 'STATUS','EXEC.','SCF','GEOM.','TOT.EN.','INT.EN.', 
                                                 'MULL.CHAR.CAT.','MULL.CHAR.AN.','D.COM.', 'D.COC.', 
                                                 'OLD.COUL.EN.','NEW.COUL.EN.' ]]
    eda_columns = [      ('EDA', ii) for ii in [ 'STATUS','EXEC.','SCF','ES.','EX.','REP.','POL.','INT.EN.' ]]
    ene_columns = [      ('ENE', ii) for ii in [ 'STATUS','EXEC.','SCF','TOT.EN.','INT.EN.' ]]
    
    scan_dict = { 'OPTIMIZE' : { 'csv' : IL.opt_scan_csv, 'res_cols' : opt_columns},
                       'EDA' : { 'csv' : IL.eda_scan_csv, 'res_cols' : eda_columns},
                    'ENERGY' : { 'csv' : IL.ene_scan_csv, 'res_cols' : ene_columns} }
    opt_csv = scan_dict[runtyp]['csv']
    opt_res_cols = scan_dict[runtyp]['res_cols']
    opt_label = runtyp[:3]

    reduced_coordinates = False
    if scan_type == 'SCAN':
       results_columns = opt_res_cols
       reduced_coordinates = False #True
       reduced_coordinates = True
       if reduced_coordinates:
          #input_df = input_df.sample(n=100)
          reduced_R_list = [str(r) for r in R_LIST]
          reduced_T_list = [str(t) for t in T_LIST] 
          reduced_P_list = [str(p) for p in P_LIST] 
          reduced_R_list = [ '3.0', '4.0', '5.0', '6.0', '7.0', '8.0', '10.0', '12.0', '15.0', '20.0' ]
          #reduced_R_list = [ '5.0' ] #'5.0', '5.5', '6.0', '6.5', '7.0', '10.0', '15.0', '20.0' ]
          reduced_T_list = [ '90' ]# 
          reduced_P_list = [ '90' ]# 
    elif scan_type == 'EQUILIBRIUM':
       results_columns = equil_columns
    elif scan_type == 'MIN.SURF.':
       results_columns = min_surf_columns
    elif scan_type == 'GPR':
       reduced_coordinates = True
       pass
    elif scan_type == 'RANDOM':
       reduced_coordinates = True

    if reduced_coordinates:
       input_df = input_df.loc[input_df[('Coordinates','Radius')].isin( reduced_R_list ) ] 
       input_df = input_df.loc[input_df[('Coordinates','Theta')].isin(  reduced_T_list ) ] 
       input_df = input_df.loc[input_df[('Coordinates','Phi')].isin(    reduced_P_list ) ] 
       print_tab( 0, 'reduced R_list: {}'.format(reduced_R_list) )
       print_tab( 0, 'reduced R_list: {}'.format(reduced_T_list) )
       print_tab( 0, 'reduced P_list: {}'.format(reduced_P_list) )

    ## read out CSV
    if os.path.exists( result_csv ):
       result_df = pd.read_csv( result_csv, header = [0,1], dtype=object )
    else:
      result_df = pd.DataFrame( columns = pd.MultiIndex.from_tuples(coords_columns + results_columns), dtype=object )
      print_tab(1, 'NEW file: {}'.format( result_csv ) )
    ## READ previously written results ends
 
    print_tab( 1, ['scan_type = {}, runtyp = {}, opt_label = {}'.format( scan_type, runtyp, opt_label), 
                   'reduced_coordinates = {}, result_csv = {}'.format( reduced_coordinates, result_csv )] )

    print( input_df )
    for index, row in input_df.iterrows():
        R, T, P  = str(row[('Coordinates','Radius')]), str(row[('Coordinates','Theta')]), str(row[('Coordinates','Phi')])
        run_dict = {('Coordinates','Radius'):str(R), ('Coordinates','Theta'):str(T), ('Coordinates','Phi'):str(P) } 
        prev_out_line, prev_out_stat, new_out_line, df_out_index = read_df_line(R,T,P,result_df,opt_label,print_opt=True )
        run_obj = get_gms_info( R, T, P, runtyp , scan_type ) 
        print(run_obj)
        
        proceed = True
        
        ### system is running
        #if run_obj.inp_name in run_labels:
        #   run_exec = 'Running'
        #   print_tab(1, 'Running')
        #   proceed = False

        ## read previous result R,T,P for optimize case
        if runtyp in ['EDA', 'ENERGY']:
           opt_obj = get_gms_info( R, T, P, 'OPTIMIZE', 'SCAN' ) 
           print( 'read opt' ) 
           prev_opt_line, prev_opt_stat, new_opt_line, df_opt_index = read_df_line(R,T,P,opt_df,'OPT')
           if prev_opt_stat == 'FAIL.':
              ## Do not calculate EDA and ENERGY if OPT did not converge
              print_tab(1, 'skipping')
              proceed = False

        ## Has the calculation already succeede/failed
        if prev_out_stat in ['FAIL.', 'SUCC.']:
           proceed = False

        #print( 'proceed:', proceed ) 
        if proceed:
           if os.path.exists( run_obj.run_dir ):
              run_exec = run_obj.get_execution(run_job_labels)
              print( run_exec )
              ## check whether input and output have the same parameters
              # run_obj.read_RTP_io( RADIUS_STR, THETA_STR, PHI_STR, CHECK_LINE )
              ## READ OUTPUT BEGINS
              if run_exec == 'Running':
                 run_stat = 'Running'
                 print_tab(3 , run_stat )
              elif run_exec == 'Missing':
                 ## folder exists but no out file
                 run_exec = resubmit_jobs( run_obj.run_dir )
              elif run_exec in ['unknown.exec.', '-ABNORMALLY-']:
                 run_stat = run_obj.read_error() # 'FAIL.'
              elif run_exec == 'NORMALLY':
                 run_stat = 'unknown.stat.'
                 gamess_out_dict = run_obj.get_out_dict()[opt_label]

                 if runtyp == 'OPTIMIZE':
                    opt_scf  = gamess_out_dict['SCF'] 
                    opt_geom = gamess_out_dict['GEOM.']
                    run_dict[('OPT','SCF' )]    = opt_scf
                    run_dict[('OPT','GEOM.')]   = opt_geom
                    run_dict[('OPT','TOT.EN.')] = gamess_out_dict['TOT.EN.'] 
                    run_dict[('OPT','INT.EN.')] = gamess_out_dict['INT.EN.']
                    
                    cart_dict = gamess_out_dict['CART.COORDS.'] 
                    atoms_dict = gamess_out_dict['ATOMS'] 
                    mulliken_dict = gamess_out_dict['MULL.CHARGES']
                    internuc_dict = gamess_out_dict['INTERNUCL.DISTANCES']
                     
                    ## calculate Distance center of masses (com) and charges (coc)
                    Z_cat, xyz_cat, charge_cat = [], [], []
                    for at_index in range(0, cat_atnum):
                        Z_cat.append( float(atoms_dict[at_index]['Z']) )
                        xyz_cat.append( [ float(cart_dict[at_index][idx]) for idx in ['x','y','z'] ] )
                        charge_cat.append( float(mulliken_dict[at_index]['charge']) )
                    com_cat = calculate_center_of_mass( Z_cat, xyz_cat )
                    coc_cat = calculate_center_of_charge( Z_cat, xyz_cat, charge_cat )
                    
                    Z_ani, xyz_ani, charge_ani = [], [], []
                    for at_index in range(cat_atnum, cat_atnum+ani_atnum):
                        Z_ani.append( float(atoms_dict[at_index]['Z']) )
                        xyz_ani.append( [ float(cart_dict[at_index][idx]) for idx in ['x','y','z'] ] )
                        charge_ani.append( float(mulliken_dict[at_index]['charge']) )
                    com_ani = calculate_center_of_mass( Z_ani, xyz_ani )
                    coc_ani = calculate_center_of_charge( Z_ani, xyz_ani, charge_ani )

                    distance_center_of_masses = LA.norm( com_ani - com_cat )
                    distance_center_of_charges = LA.norm( coc_ani - coc_cat )

                    run_dict[('OPT','D.COM.')] = distance_center_of_masses  
                    run_dict[('OPT','D.COC.')] = distance_center_of_charges 

                    ## OLD Mulliken charges
                    cat_charges_dict = { k:v for (k,v) in mulliken_dict.items() if k < cat_atnum }
                    ani_charges_dict = { k:v for (k,v) in mulliken_dict.items() if k >= cat_atnum }
                    cat_charges_sum  = sum([v['charge'] for v in cat_charges_dict.values() ])
                    ani_charges_sum  = sum([v['charge'] for v in ani_charges_dict.values() ])
                    run_dict[('OPT','MULL.CHAR.CAT.')] = cat_charges_sum 
                    run_dict[('OPT','MULL.CHAR.AN.')]  = ani_charges_sum  
                    run_dict[('OPT','OLD.COUL.EN.')] = Coulomb_Energy(distance_center_of_masses, cat_charges_sum, ani_charges_sum)

                    ## NEW Mulliken charges
                    tot_int = 0
                    tot_coul_energy = 0
                    for int_key, int_dict in internuc_dict.items():
                      ii_idx = int(int_dict['at.1']['idx.1'])
                      ii_sym = int_dict['at.1']['elem.1']
                      jj_idx = int(int_dict['at.2']['idx.2'])
                      jj_sym = int_dict['at.2']['elem.2']
                      ii_label = '{}{}'.format( ii_sym, ii_idx ) 
                      jj_label = '{}{}'.format( jj_sym, jj_idx ) 
                      print( cat_atoms_list )
                      print( ani_atoms_list )
                      if ii_label in cat_atoms_list and jj_label in ani_atoms_list:
                         ii_dict = [v for v in mulliken_dict.values() if int(v['idx'])==int(ii_idx) and v['elem.']==ii_sym ][0]
                         jj_dict = [v for v in mulliken_dict.values() if int(v['idx'])==int(jj_idx) and v['elem.']==jj_sym ][0]
                         ii_charge = ii_dict['charge']
                         jj_charge = jj_dict['charge']
                         ii_jj_dist = int_dict['dist.']
                         ii_jj_coul_energy = Coulomb_Energy( ii_jj_dist, q1 = ii_charge, q2 = jj_charge )
                         tot_coul_energy += ii_jj_coul_energy
                         tot_int += 1

                    run_dict[('OPT','NEW.COUL.EN.')] = tot_coul_energy
                    #run_dict[('OPT','NEW.COUL.EN.')] = run_obj.internal_coulomb_energy()
                    if run_exec == 'NORMALLY' and opt_scf == 'CONVERGED' and opt_geom == 'LOCATED':
                       run_stat = 'SUCC.'
                    else:
                       run_stat = 'FAIL.'

                 elif runtyp == 'EDA':
                    run_dict[ ('EDA','INT.EN.')] = gamess_out_dict['ALL.BS']['INT.EN.']
                    run_dict[ ('EDA','ES.')]     = gamess_out_dict['ALL.BS']['ES.']
                    run_dict[ ('EDA','EX.')]     = gamess_out_dict['ALL.BS']['EX.']
                    run_dict[ ('EDA','REP.')]    = gamess_out_dict['ALL.BS']['REP.']
                    run_dict[ ('EDA','POL.')]    = gamess_out_dict['ALL.BS']['POL.']
                    if run_exec == 'NORMALLY':
                       run_stat = 'SUCC.'
                    else:
                       run_stat = 'FAIL.'

                 elif runtyp == 'ENERGY':
                    ene_scf  = gamess_out_dict['SCF'] 
                    run_dict[('ENE','SCF' )] = ene_scf
                    run_dict[('ENE','TOT.EN.')] = gamess_out_dict['TOT.EN.'] 
                    run_dict[('ENE','INT.EN.')] = gamess_out_dict['INT.EN.']
                    if run_exec == 'NORMALLY' and ene_scf == 'CONVERGED':
                       run_stat = 'SUCC.'
                    else:
                       run_stat = 'FAIL.'

              run_dict[(opt_label,'STATUS')] = run_stat
              ## READ OUTPUT ENDS

           else:
              print( scan_type )
              ## folder does not exists
              if runtyp == 'OPTIMIZE':
                 run_exec = write_and_run_new( run_obj.gms_dir, run_obj.inp_name, runtyp, scan_type = scan_type )
              else:
                 prev_opt_line, prev_opt_stat, new_opt_line, df_opt_index = read_df_line( R, T, P, opt_df, 'OPT' )
                 if prev_opt_stat == 'SUCC.':
                    opt_zmat =  opt_obj.get_out_dict()['OPT']['ZMAT'] 
                    run_exec = write_and_run_new( run_obj.gms_dir, run_obj.inp_name, runtyp, scan_type = scan_type,
                                                  run_zmat = opt_zmat, nat_cat = cat_atnum, nat_ani = ani_atnum )
           
                 else:
                    run_exec = 'OPT.FAIL.'
           run_dict[(opt_label,'EXEC.')]  = run_exec
           if new_out_line:
             result_df = result_df.append( pd.Series( run_dict ), ignore_index=True )
           else:
             result_df.replace( result_df.loc[df_out_index], pd.Series( run_dict ) )

           ### copy reduced grid results from scan csv file
           #if scan_type == 'GPR':
           #   if gpr_idx == '0':
           #      link_r = '{:3.1f}'.format(float(R))
           #      link_p = int(P)
           #      link_t = int(T)
           #      run_obj = get_gms_info(tmp_r=link_r,tmp_p=link_p,tmp_t=link_t) 
           ### copy reduced grid results from scan csv file

    print( 'writing {} reuslts to {}'.format(runtyp, result_csv) )
    result_df.to_csv( result_csv, index=False )

     ### subtract ElectroStatic Contribution
     #tot_df = pd.concat([opt_df, eda_df], axis=1)
     #print( tot_df )
     #tot_df[('OPT','NON.COUL.EN.')] = tot_df[('OPT', 'TOT.EN.')] - tot_df[('EDA','ES.')]
     #print( tot_df[('OPT','NON.COUL.EN.')] )
     #print( here2 )
     #opt_df[('OPT','NON.COUL.EN.')] = tot_df[('OPT','NON.COUL.EN.')] 
     #print( opt_df[('OPT','NON.COUL.EN.')] ) #= opt_df[('OPT', 'TOT.EN.')] - eda_df[('EDA','ES.')]
     #opt_df.to_csv( result_csv, index=False )

    ## Check error in OPT DF
    if os.path.exists( IL.opt_scan_csv ):
        all_opt_df = pd.read_csv( IL.opt_scan_csv, header = [0,1], dtype=object )
        suc_opt_df = all_opt_df.loc[ all_opt_df[('OPT','STATUS')] == 'SUCC.' ]
        run_opt_df = all_opt_df.loc[ all_opt_df[('OPT','EXEC.')]  == 'Running' ]
        abn_opt_df = all_opt_df.loc[ all_opt_df[('OPT','EXEC.')] == '-ABNORMALLY-' ]
        #unknown_opt_df   = all_opt_df.loc[ all_opt_df[('OPT','STATUS')] == 'unknown.stat.' ]

        print( '{:4d} Total configurations'.format(len(all_opt_df)))      
        print( '{:4d} Succeeded       jobs'.format(len(suc_opt_df)))
        print( '{:4d} Running         jobs'.format(len(run_opt_df)))
        print( '{:4d} Abnormally      jobs'.format(len(abn_opt_df)))
        #print( '{:4d} Unknown status  jobs'.format(len(unknown_opt_df)))
        ## fix errors

        ##for fail_df in [abnormal_opt_df, unknown_opt_df]:
        #for index, fail_row in abn_opt_df.iterrows():
        #    fail_stat = fail_row[('OPT', 'STATUS')]
        #    fail_exex = fail_row[('OPT', 'EXEC.')]
        #    fail_R = fail_row[('Coordinates', 'Radius')]
        #    fail_T = fail_row[('Coordinates', 'Theta')]
        #    fail_P = fail_row[('Coordinates', 'Phi')]
        #    fail_out_line, fail_out_stat, fail_out_line, df_fail_index = read_df_line( fail_R, fail_T, fail_P, abn_opt_df, 'OPT', print_opt=True )
        #    fail_obj = get_gms_info( fail_R, fail_T, fail_P, 'OPTIMIZE', 'SCAN' )
        #    ## double check for execution
        #    fail_exec = fail_obj.get_execution(run_job_labels)
        #    if fail_exec == '-ABNORMALLY-':
        #       fail_stat = fail_obj.read_error()
        #       if fail_stat in ['reduce QMTTOL']:
        #          # find closest succeeded configuration with same T,P
        #          fixed_angles_df = suc_opt_df.loc[ suc_opt_df[('Coordinates', 'Theta')] == fail_T 
        #                                     ].loc[ suc_opt_df[('Coordinates','Phi')] == fail_P ]
        #          suc_R_list = fixed_angles_df[('Coordinates','Radius')].values
        #          delta = 1e2
        #          for suc_R in suc_R_list:
        #              tmp_delta = abs(float(suc_R) - float(fail_R))
        #              if tmp_delta < delta:
        #                 closest_R = suc_R
        #                 delta = tmp_delta
        #          print( fail_R, closest_R, delta )
        #       else:
        #          print( fail_stat )
        #    else:
        #      print( porca_merda_df_file_corrotto )
               
if __name__ == '__main__':
  main()
