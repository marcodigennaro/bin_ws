#!/home/mdi0316/anaconda3/bin/python

### common input start
import os, sys, re
import shutil
import numpy as np
from numpy import linalg as LA
import pandas as pd
import subprocess as sp

import getpass
user = getpass.getuser()

scripts_dir = '/home/{}/FUNCTIONS'.format(user)
classes_dir = '/home/{}/CLASSES'.format(user)
zmat_converter_dir = '/home/{}/CLASSES/zmatrix-master'.format(user)

sys.path.insert(0, scripts_dir)
sys.path.insert(0, classes_dir)
sys.path.insert(0, zmat_converter_dir)

import json
import math
import ast
from collections import defaultdict
from mendeleev import element

import GAMESS
import SLURM

import IONIC_LIQUID as IL

from Functions import print_tab, running_jobs, compose_zmatrices, running_label, center_of_charge, center_of_mass, Coulomb_Energy
from GAMESS import functionals_list, gbasis_list, full_R_list, full_T_list, full_P_list, full_functionals_list, full_gbasis_list
 
#from monomers import change_all_file_names 

if user == 'mdi0316':
   work_dir = '/data/{}/WORK'.format(user)
else:
   work_dir = '/data/scratch-no-backup/{}/WORK'.format(user)

dimers_dir = os.path.join( work_dir, 'DIMERS' )
os.makedirs( dimers_dir, exist_ok = True )
temp_dir = '/home/{}/Inputfiles/GAMESS/MONOMERS/AVOGADRO/'.format(user)
mono_json = os.path.join( work_dir, 'monomers_{}.json'.format(user) )

with open(mono_json,'r') as json_file:
  mono_dict = json.load(json_file)

global VERBOSE
VERBOSE = True 
VERBOSE = False

make_equil = True
make_equil = False

make_scan = False
make_scan = True

if len(sys.argv) == 1:
  DIMER_LABEL = 'EMIM_BF4'
else:
  DIMER_LABEL = sys.argv[1]
  
CAT_LABEL, ANI_LABEL = DIMER_LABEL.split('_')
print(CAT_LABEL)
print(ANI_LABEL)
 
R_EQUIL_LIST = [ '3.0', '4.0', '5.0', '6.0', '7.0', '8.0', '9.0', '10.0' ]
R_EQUIL_LIST = [ '10.0', '9.0' ] #, '8.0', '7.0', '6.0', '5.0' ]

global T_SCAN_LIST
T_SCAN_LIST = [ '5', '45', '90', '135', '175' ]
P_SCAN_LIST = [ '0', '45', '90', '135', '180', '225', '270', '315' ]
R_SCAN_LIST  = [ '2.0', #'2.1', '2.2', '2.3', '2.4', 
                 '2.5', '2.6', '2.7', '2.8', '2.9',
                 '3.0', '3.1', '3.2', '3.3', '3.4', 
                 '3.5', '3.6', '3.7', '3.8', '3.9',
                 #'3.75', '3.85', '3.95',
                 #'4.05', 
                 '4.0', '4.1', '4.2', '4.3', '4.4', 
                 '4.5', '4.6', '4.7', '4.8', '4.9', 
                 '5.0', '5.5', '6.0', '6.5', '7.0', 
                 '7.5', '8.0', '8.5', '9.0', '9.5', 
                 '10.0', '11.0', '12.0', '13.0', '14.0', '15.0' ]

print(T_SCAN_LIST)
#reduced_T_SCAN_LIST = [ '90' ]
#reduced_P_SCAN_LIST = [ '90' ]
reduced_T_SCAN_LIST = [ '5', '45', '90', '135', '175' ]
reduced_P_SCAN_LIST = [ '0', '45', '90', '135', '180', '225', '270', '315' ]

reduced_R_SCAN_LIST = [ '2.0', '2.5', '3.0', '3.5', '4.0', '4.5', '5.0', '5.5', '6.0', '6.5', '7.0',  
                        '7.5', '8.0', '8.5', '9.0', '9.5', '10.0', '11.0', '12.0', '13.0', '15.0' ]

#full_gbasis_list      = [ 'APCseg-1', 'STO',  'N311' ]

full_gbasis_list      = [ 'N311' ]
full_functionals_list = [ 'B3LYP']

#full_functionals_list = [ 'PBE0', 'B3LYP', 'M11' , 'wB97x-D' ]

if DIMER_LABEL not in ['EMIM_BF4', 'EMIM_PF6']:
   full_gbasis_list = ['N311']
   full_functionals_list = ['B3LYP']
   T_SCAN_LIST = [ '5', '45', '90', '135', '175' ]
   P_SCAN_LIST = [ '0', '45', '90', '135', '180', '225', '270', '315' ]
   
   #T_SCAN_LIST = ['90']
   #P_SCAN_LIST = ['90']

print( 'R_SCAN_LIST: {}'.format( R_SCAN_LIST ))
print( 'T_SCAN_LIST: {}'.format( T_SCAN_LIST ))
print( 'P_SCAN_LIST: {}'.format( P_SCAN_LIST ))
print( 'full_gbasis_list: {}'.format( full_gbasis_list ))
print( 'full_functionals_list: {}'.format( full_functionals_list ))

global READ_FROM
READ_FROM = 'ISOLATED'


RTP_columns = ['Radius', 'Theta', 'Phi']

dft_columns = RTP_columns + ['TOT.EN.', 'INT.EN.', 'DISP.EN.', 'CHARG.ANI.', 'CHARG.CAT.', 'COUL.EN.', 'DIST.COM', 
                              'BASIS.DIM.', 'INERT.MOM', 'COM', 'COC', 'Run.Time']
mp2_columns = RTP_columns + ['TOT.EN.', 'MP2.EN.', 'DISP.EN.', 'CHARG.ANI.', 'CHARG.CAT.', 'COUL.EN.', 'DIST.COM']
ccs_columns = RTP_columns + ['TOT.EN.', 'CCSD(T).EN.', 'DISP.EN.', 'CHARG.ANI.', 'CHARG.CAT.', 'COUL.EN.', 'DIST.COM']
err_columns = RTP_columns + ['ERROR']
frc_columns = RTP_columns + ['f_x', 'f_y', 'f_z']

def write_pd_series( R, T, P, scan_out_dict, scan_inp_dict, post_proc=False, equil=False ):

    mull_charges = scan_out_dict['MULL.CHARGES'] 
    cart_coords  = scan_out_dict['CART.COORDS.']
    cart_dict = { 'Radius' : R, 'Theta' : T, 'Phi' : P, 'cart.coords.' : cart_coords, 'mull.charges' : mull_charges }
    ##
    com = center_of_charge( cart_coords, mull_charges ) 
    coc = center_of_mass( cart_coords ) 
    cat_coords = dict(( str(k), cart_coords[k]) for k in range(CAT_NAT) )
    ani_coords = dict(( str(k), cart_coords[k]) for k in range(CAT_NAT, CAT_NAT + ANI_NAT ) )
    cat_com = center_of_mass( cat_coords )
    ani_com = center_of_mass( ani_coords )
    dcom = LA.norm( np.array(cat_com)-np.array(ani_com ))

    pd_dict = { 'Radius' : R, 'Theta' : T, 'Phi' : P }

    ## common
    pd_dict['TOT.EN.'] = scan_out_dict['TOT.EN.']
    pd_dict['CHARG.CAT.'] = scan_out_dict['CHARG.CAT.']
    pd_dict['CHARG.ANI.'] = scan_out_dict['CHARG.ANI.'] 
    pd_dict['DIST.COM'] = dcom
    pd_dict['COUL.EN.'] = Coulomb_Energy(float(pd_dict['DIST.COM']), float(pd_dict['CHARG.CAT.']), float(pd_dict['CHARG.ANI.']))

    if not post_proc:
       pd_dict['INT.EN.'] = scan_out_dict['INT.EN.']
       pd_dict['DISP.EN.'] = scan_out_dict['INT.EN.'] - pd_dict['COUL.EN.']
       pd_dict['BASIS.DIM.'] = scan_out_dict['BASIS.DIM.']
       pd_dict['INERT.MOM']  = scan_out_dict['INERT.MOM.']
       pd_dict['COM'] = com
       pd_dict['COC'] = coc

    if post_proc == 'MP2':
       pd_dict[ 'MP2.EN.'] = scan_out_dict['MP2']['MP2.EN.'] - ZERO_MP2_ENER
       pd_dict['DISP.EN.'] = pd_dict['MP2.EN.'] - pd_dict['COUL.EN.']

    elif post_proc in ['EDA', 'HES']:
      print( 'write_out' )
    else:
      if equil:
        pd_dict['Relax.Radius'] = scan_out_dict['FINAL']['ZMAT'][19]['STR']['val']

    pd_dict['Run.Time'] = scan_out_dict['TIME']

    pd_series = pd.Series( pd_dict ) #, dtype=object )
    cart_series = pd.Series( cart_dict ) #, dtype=object )
    return( pd_series, cart_series )

def get_gms_object( basis, funct, T, P, R, equil = False, opt_method = 'QA', post_scf = 'DFTTYP', run_type = 'OPTIMIZE' ):
    #tmp_ifreeze = '52,53,54' # EMIM_BF4
    cat_zmat_dim = 3*CAT_NAT-6
    if equil:
       TPR_CONF = IL.DIMER_EQUIL_CONF( DIMER_LABEL, basis, funct, T=T, P=P, R=R )
       TPR_label = 'EQUIL_{}_T_{}_P_{}_R_{}_{}_{}'.format(DIMER_LABEL.lower(), T, P, R, basis, funct )
       tmp_ifreeze = '{},{}'.format( cat_zmat_dim + 1,  cat_zmat_dim + 2 )
    else:
       TPR_CONF = IL.DIMER_SCAN_CONF( DIMER_LABEL, basis, funct, T=T, P=P, R=R )
       TPR_label = 'SCAN_from_ISOLATED_{}_T_{}_P_{}_R_{}_{}_{}'.format(DIMER_LABEL.lower(), T, P, R, basis, funct )
       tmp_ifreeze = '{},{},{}'.format( cat_zmat_dim + 1,  cat_zmat_dim + 2, cat_zmat_dim + 3 ) 
 
    gms_obj = GAMESS.GAMESS( inp_label = TPR_label, root_dir = TPR_CONF.R_dir, 
                             natoms = CAT_NAT + ANI_NAT, nat_cat = CAT_NAT, nat_ani = ANI_NAT, 
                             icharge = 0, zero_energy = ZERO_DFT_ENER,
                             run_type = run_type, post_scf = post_scf,
                             basis = basis, functional = funct,
                             ifreeze = tmp_ifreeze, opt_method = opt_method )

    if VERBOSE:
       print_tab( 3, gms_obj.run_dir )
    return( TPR_label, gms_obj )

def read_dimer( basis, funct ):

    DIMER = IL.DIMER( DIMER_LABEL, basis, funct )

    CAT_NAT = DIMER.cat_dict['nat']
    ANI_NAT = DIMER.ani_dict['nat']

    CATION = IL.MONOMER( CAT_LABEL, basis, funct )
    ANION  = IL.MONOMER( ANI_LABEL, basis, funct )
      
    try:
      cat_zmat   = CATION.mono_dict['OUT'][basis][funct]['DFT']['FINAL']['ZMAT']
      CAT_DFT_EN = CATION.mono_dict['OUT'][basis][funct]['DFT']['TOT.EN.'] 

      ani_zmat   =  ANION.mono_dict['OUT'][basis][funct]['DFT']['FINAL']['ZMAT']
      ANI_DFT_EN =  ANION.mono_dict['OUT'][basis][funct]['DFT']['TOT.EN.']

      ZERO_DFT_ENER = CAT_DFT_EN + ANI_DFT_EN

      ANI_MP2_EN = ANION.mono_dict['OUT'][basis][funct]['MP2']['MP2.EN.']
      CAT_MP2_EN = CATION.mono_dict['OUT'][basis][funct]['MP2']['MP2.EN.']
      ZERO_MP2_ENER = CAT_MP2_EN + ANI_MP2_EN
    except(KeyError):
      print( 'Missing Monomer' )
      cat_zmat, ani_zmat, ZERO_DFT_ENER, ZERO_MP2_ENER = False, False, False, False    
    
    return( DIMER, CAT_NAT, ANI_NAT, cat_zmat, ani_zmat, ZERO_DFT_ENER, ZERO_MP2_ENER )

def post_process( pp_label, pp_obj, dft_zmat, R, T, P, pp_jq = 'nodeshiq' ):
    pp_series = None
    running = running_label( pp_obj.inp_name )
    if running:
       print_tab( 4, '{}: Running'.format(pp_label) )
    else:
       if os.path.exists( pp_obj.inp_file ):
          pp_exec, pp_exec_err = pp_obj.get_job_exec()
          if [ pp_exec, pp_exec_err ] == ['MISSING.OUTPUT.FILE', False]:
             pp_obj.resubmit()
          elif [ pp_exec, pp_exec_err ] == ['TERMINATED.NORMALLY', False]:
             pp_inp_dict, pp_out_dict, pp_scf, pp_geom = pp_obj.get_job_results()
             if pp_scf == 'CONVERGED':
                print_tab( 4, '{}: ok (new)'.format(pp_label) )
                pp_series = write_pd_series( R, T, P, pp_out_dict, pp_inp_dict, post_proc = pp_label )[0]
             else:
                print_tab( 4, '{}: not ok (new)'.format(pp_label) )
                pp_obj.fix_error()
          else:
             print_tab( 4, '{}: EXEC: {}, ERR: {}'.format(pp_label, pp_exec, pp_exec_err ) )
             pp_obj.fix_error()
       else:
          print_tab( 4, 'Submitting {}'.format(pp_label) )
          pp_obj.run_new( zmat_dict = dft_zmat, job_queue = pp_jq )

    return pp_series


def main():

    global T_SCAN_LIST
    global P_SCAN_LIST
    global R_SCAN_LIST
    global CAT_NAT
    global ANI_NAT
    global ZERO_DFT_ENER
    global ZERO_MP2_ENER

    print_tab( 0, '>>>> {} <<<<'.format(DIMER_LABEL) )

    for tmp_basis in full_gbasis_list:
      for tmp_funct in full_functionals_list:
         print_tab( 1, '=== {} ==='.format(tmp_basis) )
         print_tab( 2, '=== {} ==='.format(tmp_funct) )

         DIMER, CAT_NAT, ANI_NAT, cat_zmat, ani_zmat, ZERO_DFT_ENER, ZERO_MP2_ENER = read_dimer( tmp_basis, tmp_funct )

         os.makedirs( DIMER.runs_dir, exist_ok=True )
         os.makedirs( DIMER.csv_dir,  exist_ok=True )

         if [cat_zmat, ani_zmat, ZERO_DFT_ENER, ZERO_MP2_ENER] == [False, False, False, False]:
           proceed = False

         else:
           proceed = True

           eq_err_df = pd.DataFrame( columns = err_columns )
           eq_dft_df = pd.DataFrame( columns = RTP_columns + [ 'Relax.Radius' ] )
           eq_mp2_df = pd.DataFrame( columns = RTP_columns ) 
           err_df = pd.DataFrame( columns = err_columns )
           dft_df = pd.DataFrame( columns = dft_columns )
           mp2_df = pd.DataFrame( columns = mp2_columns )
           ccs_df = pd.DataFrame( columns = ccs_columns )
           eda_df = pd.DataFrame( columns = RTP_columns )
           hes_df = pd.DataFrame( columns = RTP_columns )
           crd_df = pd.DataFrame( columns = [ 'Radius', 'cart.coords.', 'mull.charges' ] ) 
           frc_df = pd.DataFrame( columns = frc_columns )
              
           if os.path.exists( DIMER.equil_err_csv ):
              eq_err_df = pd.read_csv( DIMER.equil_err_csv, index_col=0, dtype=object )
           if os.path.exists( DIMER.equil_dft_csv ):
              eq_dft_df = pd.read_csv( DIMER.equil_dft_csv, index_col=0, dtype=object )
           if os.path.exists( DIMER.equil_mp2_csv ):
              eq_mp2_df = pd.read_csv( DIMER.equil_mp2_csv, index_col=0, dtype=object )
           if os.path.exists( DIMER.scan_err_csv ):
              err_df = pd.read_csv( DIMER.scan_err_csv, index_col=0, dtype=object )
           if os.path.exists( DIMER.scan_dft_csv ):
              dft_df = pd.read_csv( DIMER.scan_dft_csv, index_col=0, dtype=object )
           if os.path.exists( DIMER.scan_mp2_csv ):
              mp2_df = pd.read_csv( DIMER.scan_mp2_csv, index_col=0, dtype=object )
           if os.path.exists( DIMER.scan_ccs_csv ):
              ccs_df = pd.read_csv( DIMER.scan_ccs_csv, index_col=0, dtype=object )
           if os.path.exists( DIMER.scan_eda_csv ):
              eda_df = pd.read_csv( DIMER.scan_eda_csv, index_col=0, dtype=object )
           if os.path.exists( DIMER.scan_hes_csv ):
              hes_df = pd.read_csv( DIMER.scan_hes_csv, index_col=0, dtype=object )
           if os.path.exists( DIMER.scan_crd_csv ):
              crd_df = pd.read_csv( DIMER.scan_crd_csv, index_col=0, dtype=object )
           if os.path.exists( DIMER.scan_frc_csv ):
              frc_df = pd.read_csv( DIMER.scan_frc_csv, index_col=0, dtype=object )

           ##############
           if make_equil:
              print_tab( 3, 'EQUILIBRIUM' )
    
              for T in equil_T_list:
                for P in equil_P_list:
                  for R in R_EQUIL_LIST:
                    
                    err_line = eq_err_df.loc[ eq_err_df['Radius']==str(R) ].loc[ eq_err_df[ 'Theta']==str(T) ].loc[ eq_err_df['Phi']==str(P) ]
                    dft_line = eq_dft_df.loc[ eq_dft_df['Radius']==str(R) ].loc[ eq_dft_df[ 'Theta']==str(T) ].loc[ eq_dft_df['Phi']==str(P) ]
                    mp2_line = eq_mp2_df.loc[ eq_mp2_df['Radius']==str(R) ].loc[ eq_mp2_df[ 'Theta']==str(T) ].loc[ eq_mp2_df['Phi']==str(P) ]
                    if err_line.empty:
                      if dft_line.empty or mp2_line.empty:
                         run_mp2 = False
                         run_hes = False
                         print_tab( 3, 'T = {}, P = {}, R = {}'.format(T,P,R) )
                         eq_label, eq_opt_obj = get_gms_object( tmp_basis, tmp_funct, T, P, R, equil = True )
                         eq_label, eq_mp2_obj = get_gms_object( tmp_basis, tmp_funct, T, P, R, equil = True, post_scf = 'MP2', run_type = 'ENERGY')
                         eq_label, eq_hes_obj = get_gms_object( tmp_basis, tmp_funct, T, P, R, equil = True, post_scf = 'DFTTYP', run_type = 'HESSIAN')
                         running = running_label( eq_opt_obj.inp_name )
                         if running:
                            print_tab( 4, 'Running' )
                            break ## skip till this is finished
                         else:
                            if os.path.exists( eq_opt_obj.inp_file ):
                               eq_exec, eq_exec_err = eq_opt_obj.get_job_exec()
                               if [ eq_exec, eq_exec_err ] == ['MISSING.OUTPUT.FILE', False]:
                                  eq_obj.resubmit()
                               if [eq_exec, eq_exec_err] == ['TERMINATED.NORMALLY', False]:
                                  eq_inp_dict, eq_out_dict, eq_scf, eq_geom = eq_opt_obj.get_job_results()
                                  if [ eq_scf, eq_geom ] == ['CONVERGED', 'LOCATED']:
                                     print_tab( 4, 'OPT.EQ. ok' )
                                     run_mp2 = True
                                     run_ccs = True
                                     if dft_line.empty:
                                        eq_series = write_pd_series( R, T, P, eq_out_dict, eq_inp_dict, equil=True )[0]
                                        eq_dft_df = eq_dft_df.append( eq_series , ignore_index=True )
                                  else:
                                     print_tab( 4, 'OPT.EQ. not ok' )
                                     eq_opt_obj.fix_error()
                                     eq_opt_err = eq_opt_obj.read_error()
                                     #if eq_opt_err in [ 'atoms.too.close', 'Stationary.Point.Location.failed' ]:
                                     #   failed_serie = pd.Series({ 'Radius':R, 'Theta':T, 'Phi':P, 'TOT.EN.': eq_opt_err} )
                                     #   eq_dft_df = eq_dft_df.append( failed_serie, ignore_index=True )
                                     #   eq_mp2_df = eq_mp2_df.append( failed_serie, ignore_index=True )
                               else:
                                  print_tab( 4, 'OPT.EQ. FAILED' )
                                  eq_opt_obj.fix_error()
                            else:
                               os.makedirs( eq_opt_obj.run_dir, exist_ok=True )
                               comp_zmat = compose_zmatrices( cat_zmat, ani_zmat, radius=R , theta=T, phi=P )
                               #eq_opt_obj.run_new( zmat_dict = comp_zmat, msg='equilibrium', job_queue='nodesloq' )
                               eq_opt_obj.run_new( zmat_dict = comp_zmat, msg='equilibrium', job_queue='nodeshiq' )
                               break ## skip till this is finished

                            ## MP2 STARTS
                            if run_mp2:
                               eq_dft_zmat = eq_out_dict['FINAL']['ZMAT']
                               eq_mp2_series = post_process( 'MP2', eq_mp2_obj, eq_dft_zmat, R, T, P )
                               if isinstance(eq_mp2_series, pd.core.series.Series) :
                                  eq_mp2_df = eq_mp2_df.append( eq_mp2_series, ignore_index=True )
                            ## MP2 ENDS

                            ## HES STARTS
                            run_hes = False
                            if run_hes:
                               eq_dft_zmat = eq_out_dict['FINAL']['ZMAT']
                               eq_hes_series = post_process( 'HES', eq_hes_obj, eq_dft_zmat, R, T, P )
                               #if isinstance(eq_hes_series, pd.core.series.Series) :
                               #   eq_hes_df = eq_hes_df.append( eq_hes_series, ignore_index=True )
                            ## HES ENDS

              eq_err_df.to_csv( DIMER.equil_err_csv )
              eq_dft_df.to_csv( DIMER.equil_dft_csv )
              eq_mp2_df.to_csv( DIMER.equil_mp2_csv )


           ##############
           if make_scan:
              print_tab( 3, 'SCAN' )

              if [tmp_basis, tmp_funct] != ['N311', 'B3LYP'] or DIMER not in ['EMIM_BF4', 'EMIM_PF6']:
		 #T_SCAN_LIST = [ '5', '45', '90', '135', '175']
		 #P_SCAN_LIST = [ '0', '45', '90', '135', '180', '225', '270', '315']
	         #T_SCAN_LIST = T_SCAN_LIST
		 #P_SCAN_LIST = P_SCAN_LIST
                 T_SCAN_LIST = reduced_T_SCAN_LIST
                 R_SCAN_LIST = reduced_R_SCAN_LIST
                 P_SCAN_LIST = reduced_P_SCAN_LIST


              # ## find equilibrium with highest initial distance
              # if READ_FROM == 'HCER':
              #    eq_df = pd.read_csv( DIMER.equil_dft_csv, index_col=0, dtype=object )
              #    try:
              #      hcer = max( [ float(r) for r in eq_df['Radius'].values ] ) #highest_converged_equil_r
              #      print_tab( 2, '--- Highest converged equilibrium radius = {}'.format(hcer) )
              #    except(ValueError):
              #      print_tab( 2, '--- No converged equilibrium radius ' )
              #      print_tab( 2, '--- skip functional ' )
              #      proceed = False
              #   
              # if READ_FROM == 'ISOLATED':
              #    #pass
              #    proceed = True
              #    if [cat_zmat, ani_zmat, ZERO_DFT_ENER, ZERO_MP2_ENER] == [False, False, False, False]:
              #       proceed = False
              # elif READ_FROM == 'HCER':
              #    print_tab( 2, '--- Reading zmat from equilibrium radius = {}'.format(hcer) )
              #    guess_label, guess_obj = get_gms_object( tmp_basis, tmp_funct, T, P, hcer, equil = True )
              #    guess_exec, guess_exec_err = guess_obj.get_job_exec()
              #    print_tab( 4, 'GUESS: {}, {}'.format( guess_exec, guess_exec_err ))
              #    if [ guess_exec, guess_exec_err ] == ['TERMINATED.NORMALLY', False]:
              #       guess_inp_dict, guess_out_dict, guess_scf, guess_geom = guess_obj.get_job_results()
              #       print_tab( 4, 'GUESS: {}, {}'.format( guess_scf, guess_geom ))
              #       runs_dir = '/data/mdi0316/WORK/DIMERS/EMIM_BF4/RUNS/'
              #       read_hcer_dir = os.path.join(  runs_dir, 'SCAN_from_{}'.format(hcer), tmp_basis, tmp_funct )
              #       write_hcer_dir = os.path.join( runs_dir, 'SCAN_from_HCER', tmp_basis, tmp_funct )
              #       if not os.path.exists( write_hcer_dir ):
              #          sp.call( 'mkdir -p {}'.format( write_hcer_dir ), shell=True )
              #       print( 'rsync  {}/ {}'.format( read_hcer_dir, write_hcer_dir  ))
              #       exit() 
              #       sp.call( 'rsync  {}/ {}'.format( read_hcer_dir, write_hcer_dir ), shell=True )
              #       exit() 
              # else:
              #    print_tab( 2, '--- Reading zmat from equilibrium radius = {}'.format(READ_FROM) )
              #    if float(READ_FROM) != float(hcer):
              #       print_tab( 2, '--- skipping ' )
              #       proceed = False
              #    else:
              #       guess_label, guess_obj = get_gms_object( tmp_basis, tmp_funct, T, P, READ_FROM, equil = True )
              #       running = running_label( guess_obj.inp_name )
              #       if running:
              #          print_tab( 2, '--- equilibrium radius is Running ' )
              #          proceed = False
              #       else:
              #          guess_exec, guess_exec_err = guess_obj.get_job_exec()
              #          print_tab( 4, 'GUESS: {}, {}'.format( guess_exec, guess_exec_err ))
              #          if [ guess_exec, guess_exec_err ] == ['TERMINATED.NORMALLY', False]:
              #             guess_inp_dict, guess_out_dict, guess_scf, guess_geom = guess_obj.get_job_results()
              #             print_tab( 4, 'GUESS: {}, {}'.format( guess_scf, guess_geom ))
              #             if [ guess_scf, guess_geom ] == ['CONVERGED', 'LOCATED']:
              #                proceed = True
              #             else:
              #                proceed = False
              #          else:
              #             proceed = False
              # ## find equilibrium with highest initial distance
               
              for T in T_SCAN_LIST:
                print_tab(3, '==========================================')
                for P in P_SCAN_LIST:
                  print_tab(3, '------------------------------------------')
                  for R in R_SCAN_LIST:
                    #if tmp_basis == 'APCseg-1' and tmp_funct == 'PBE0' and R == '4.2':
                    #   pass
                    #elif [T,P] == ['5', '0'] and float(R) in [2.0, 2.9]:
                    #   pass
                    #elif [T,P] == ['5', '90'] and float(R) <= 3.8:
                    #   pass
                    #elif [T,P] == ['5', '180'] and float(R) <= 3.1:
                    #   pass
                    #elif [T,P] == ['45', '0'] and ( float(R) <= 3.8 or  4.3 <= float(R) <= 5.0 ) :
                    #   pass
                    #elif [T,P] == ['45', '90'] and float(R) <= 2.6 :
                    #   pass
                    #elif [T,P] == ['90', '0'] and  ( 2.9 <= float(R) <= 3.9 or  4.2 <= float(R) <= 4.6 ) :
                    #   pass
                    #elif [T,P] == ['135', '0'] and float(R) <= 3.5 and float(R) != 3.2 :
                    #   pass
                    #elif [T,P] == ['135', '90'] and float(R) in [6.5]:
                    #   pass
                    #elif [T,P] == ['175', '0'] and float(R) <= 5.5:
                    #   pass
                    #elif [T,P] == ['175', '90'] and float(R) <= 4.3:
                    #   pass
                    #elif [T,P] == ['175', '180'] and float(R) <= 3.2:
                    #   pass
                    #else:
                    print_tab( 3, 'T = {}, P = {}, R = {}'.format(T,P,R) )
                    run_dft = False
                    run_mp2 = False
                    run_ccs = False
                    run_hes = False
                    run_eda = False
                    dft_line = dft_df.loc[dft_df['Radius']==str(R)].loc[dft_df['Theta']==str(T)].loc[dft_df['Phi']==str(P)]
                    mp2_line = mp2_df.loc[mp2_df['Radius']==str(R)].loc[mp2_df['Theta']==str(T)].loc[mp2_df['Phi']==str(P)]
                    ccs_line = ccs_df.loc[ccs_df['Radius']==str(R)].loc[ccs_df['Theta']==str(T)].loc[ccs_df['Phi']==str(P)]
                    eda_line = eda_df.loc[eda_df['Radius']==str(R)].loc[eda_df['Theta']==str(T)].loc[eda_df['Phi']==str(P)]
                    hes_line = hes_df.loc[hes_df['Radius']==str(R)].loc[hes_df['Theta']==str(T)].loc[hes_df['Phi']==str(P)]
                    err_line = err_df.loc[err_df['Radius']==str(R)].loc[err_df['Theta']==str(T)].loc[err_df['Phi']==str(P)]

                    if err_line.empty: 
                       #if dft_line.empty or mp2_line.empty or hes_line.empty: 
                       if dft_line.empty or mp2_line.empty: 
                          scan_label, scan_obj = get_gms_object( tmp_basis, tmp_funct, T, P, R )
                          scan_label, mp2_obj  = get_gms_object(tmp_basis,tmp_funct,T,P,R,post_scf='MP2',run_type='ENERGY')
                          scan_label, ccs_obj  = get_gms_object(tmp_basis,tmp_funct,T,P,R,post_scf='CCSDT',run_type='ENERGY')
                          scan_label, hes_obj  = get_gms_object(tmp_basis,tmp_funct,T,P,R,post_scf='MP2',run_type='HESSIAN')
                          scan_label, eda_obj  = get_gms_object(tmp_basis,tmp_funct,T,P,R,post_scf='DFTTYP',run_type='EDA')

                          abort_dir = os.path.join( scan_obj.run_dir, 'FAILED', 'ABORTED' )
                          if os.path.exists( abort_dir ):
                             print_tab( 4, 'ABORTED' )
                          else:
                             if VERBOSE:
                                print_tab( 3, '{}, {}'.format(scan_obj.run_dir, scan_obj.inp_name ) )
                             #running = running_label( scan_label )
                             run_check_scan = running_label( scan_obj.inp_name )
                             run_check_mp2  = running_label( mp2_obj.inp_name  )
                             run_check_ccs  = running_label( ccs_obj.inp_name  )
                             run_check_hes  = running_label( hes_obj.inp_name  )
                             run_check_eda  = running_label( eda_obj.inp_name  )
                             if run_check_scan:
                                print_tab( 4, 'Running scan' )
                             elif run_mp2: 
                                print_check_tab( 4, 'Running mp2' )
                             elif run_hes: 
                                print_check_tab( 4, 'Running hes' )
                             elif run_eda: 
                                print_check_tab( 4, 'Running eda' )
                             else:
                                
                                if os.path.exists( scan_obj.inp_file ):
                                   scan_exec, scan_exec_err = scan_obj.get_job_exec()
                                   if [ scan_exec, scan_exec_err ] == ['MISSING.OUTPUT.FILE', False]:
                                      scan_obj.resubmit()
                                   elif [ scan_exec, scan_exec_err ] == ['TERMINATED.NORMALLY', False]:
                                      scan_inp_dict, scan_out_dict, scan_scf, scan_geom = scan_obj.get_job_results()
                                      if [ scan_scf, scan_geom ] == ['CONVERGED', 'LOCATED']:
                                         print_tab( 4, 'OPT: ok (new)' )
                                         run_mp2 = True
                                         run_ccs = True
                                         if dft_line.empty:
                                            dft_series, cart_series = write_pd_series( R, T, P, scan_out_dict, scan_inp_dict )
                                            dft_df = dft_df.append( dft_series  , ignore_index=True )
                                            crd_df = crd_df.append( cart_series , ignore_index=True )
          
                                            frc_line = pd.read_csv( scan_obj.force_file, index_col = 0 )
                                            frc_line['Radius'] = R
                                            frc_line['Theta'] = T
                                            frc_line['Phi'] = P
                                            frc_df = pd.concat( [frc_df, frc_line], sort=True ) 
                                      else:
                                         print_tab( 4, 'OPT: not ok' )
                                         scan_obj.fix_error()
                                         scan_obj_err = scan_obj.read_error()
                                         if scan_obj_err in [ 'atoms.too.close', 'Stationary.Point.Location.failed', 'Gradient.out.of.range' ]:
                                            failed_serie = pd.Series({ 'Radius' : R, 'Theta' : T, 'Phi': P, 
                                                                       'ERROR':scan_obj_err}, dtype=object )
                                            print_tab( 4, 'OPT: EXEC: {}, ERR: {}'.format( scan_exec, scan_obj_err ) )
                                            err_df = err_df.append( failed_serie, ignore_index=True )
                                   else:
                                      print_tab( 4, 'OPT: EXEC: {}, ERR: {}'.format( scan_exec, scan_exec_err ) )
                                      scan_obj.fix_error()
                                else:
                                  run_dft = True
                       else:
                          print_tab( 4, 'DFT + MP2 + EDA = OK' )

                       ## RUN NEW DFT
                       if run_dft:
                          if READ_FROM == 'ISOLATED':
                             guess_zmat = compose_zmatrices( cat_zmat, ani_zmat, radius=R , theta=T, phi=P )
                             msg = 'zmat.from.isolated.ions'
                          else:
                             guess_zmat = guessout_dict['FINAL']['ZMAT']
                             guess_zmat[19]['STR']['val'] = R
                             msg = 'zmat.from.equil.R={}'.format(READ_FROM)
                          scan_obj.run_new( zmat_dict = guess_zmat, msg=msg, job_queue='nodesloq' )
                          #scan_obj.run_new( zmat_dict = guess_zmat, msg=msg, job_queue='nodeshiq' )

                       ## MP2 STARTS
                       if run_mp2 and mp2_line.empty and not run_check_mp2:
                          scan_inp_dict, scan_out_dict, scan_scf, scan_geom = scan_obj.get_job_results()
                          dft_zmat = scan_out_dict['FINAL']['ZMAT']
                          dft_com  = dft_df.loc[ dft_df['Radius'] == R ]['DIST.COM'].values
                          mp2_series = post_process( 'MP2', mp2_obj, dft_zmat, R, T, P )
                          if isinstance(mp2_series, pd.core.series.Series) :
                             mp2_df = mp2_df.append(mp2_series,ignore_index=True)
                       ## MP2 ENDS

                       ## CCSDT STARTS
                       if run_ccs and ccs_line.empty and not run_check_ccs and [T,P] == ['90', '90'] and R == '5.0':
                          try:
                            scan_inp_dict, scan_out_dict, scan_scf, scan_geom = scan_obj.get_job_results()
                            ccs_zmat = scan_out_dict['FINAL']['ZMAT']
                            dft_com  = dft_df.loc[ dft_df['Radius'] == R ]['DIST.COM'].values
                            ccs_series = post_process( 'CCSDT', ccs_obj, dft_zmat, R, T, P )
                            if isinstance(ccs_series, pd.core.series.Series) :
                               ccs_df = ccs_df.append(ccs_series,ignore_index=True)
                          except:
                            pass
                       ## CCSDT ENDS

                       ## HES STARTS
                       run_hes = False
                       if [T, P] == ['90', '90']:
                          if run_hes and hes_line.empty and not run_check_hes:
                             scan_inp_dict, scan_out_dict, scan_scf, scan_geom = scan_obj.get_job_results()
                             dft_zmat = scan_out_dict['FINAL']['ZMAT']
                             dft_com  = dft_df.loc[ dft_df['Radius'] == R ]['DIST.COM'].values
                             hes_series = post_process( 'HES', hes_obj, dft_zmat, R, T, P )
                             if isinstance(hes_series, pd.core.series.Series) :
                                hes_df = hes_df.append(hes_series,ignore_index=True)
                             hes_obj.get_out_dict()
                       ## HES ENDS

                       ## EDA STARTS
                       run_eda = False
                       if run_eda and eda_line.empty and not run_check_eda:
                          scan_inp_dict, scan_out_dict, scan_scf, scan_geom = scan_obj.get_job_results()
                          dft_zmat = scan_out_dict['FINAL']['ZMAT']
                          dft_com  = dft_df.loc[ dft_df['Radius'] == R ]['DIST.COM'].values
                          ### DFT-D METHODS ARE NOT SUPPORTED. ###
                          eda_series = post_process( 'EDA', eda_obj, dft_zmat, R, T, P )
                          if isinstance(eda_series, pd.core.series.Series) :
                             eda_df = eda_df.append(eda_series,ignore_index=True)
                       ## EDA ENDS

#                      ## BSSE STARTS
#                      ## BSSE ENDS
                    else:
                      print_tab( 4, err_line['ERROR'].values[0] )


              for df_obj, csv_obj in zip( [dft_df, mp2_df, eda_df, hes_df, err_df, crd_df, frc_df], 
                                          [DIMER.scan_dft_csv, DIMER.scan_mp2_csv, DIMER.scan_eda_csv, 
                                           DIMER.scan_hes_csv, DIMER.scan_err_csv, DIMER.scan_crd_csv, DIMER.scan_frc_csv] ):
                  df_obj.sort_values(by=['Radius'], inplace=True )
                  df_obj.reset_index()
                  df_obj.to_csv( csv_obj )

#              dft_df.sort_values(by=['Radius'], inplace=True ).reset_index()
#              mp2_df.sort_values(by=['Radius'], inplace=True ).reset_index()
#              eda_df.sort_values(by=['Radius'], inplace=True ).reset_index()
#              hes_df.sort_values(by=['Radius'], inplace=True ).reset_index()
#              err_df.sort_values(by=['Radius'], inplace=True ).reset_index()
#              crd_df.sort_values(by=['Radius'], inplace=True ).reset_index()
#              frc_df.sort_values(by=['Radius'], inplace=True ).reset_index()
#
#              dft_df.to_csv( DIMER.scan_dft_csv )
#              mp2_df.to_csv( DIMER.scan_mp2_csv )
#              eda_df.to_csv( DIMER.scan_eda_csv )
#              hes_df.to_csv( DIMER.scan_hes_csv )
#              err_df.to_csv( DIMER.scan_err_csv )
#              crd_df.to_csv( DIMER.scan_crd_csv )
#              frc_df.to_csv( DIMER.scan_frc_csv )
 
#              ## FIXING ERROR
#              for T in T_SCAN_LIST:
#                print_tab(3, '==========================================')
#                for P in P_SCAN_LIST:
#                  print_tab(3, '------------------------------------------')
#                  tmp_err_df = err_df.loc[ err_df['Theta'] == T ].loc[ err_df['Phi'] == P ]
#                  tmp_dft_df = dft_df.loc[ dft_df['Theta'] == T ].loc[ dft_df['Phi'] == P ]
#
#                  tmp_err_r_list = list(tmp_err_df['Radius'].astype(float).values)
#                  tmp_dft_r_list = list(tmp_dft_df['Radius'].astype(float).values)
#                  tmp_err_r_list.sort(reverse = True)
#                  tmp_dft_r_list.sort(reverse = True)
#                  for tmp_err_r in tmp_err_r_list:
#                      closest_r = 100
#                      delta = 100
#                      for tmp_r in tmp_dft_r_list:
#                          tmp_delta = abs( tmp_r - tmp_err_r )
#                          if tmp_delta < delta:
#                             delta = tmp_delta
#                             closest_r = tmp_r
#                      tmp_err_label, tmp_err_obj = get_gms_object( tmp_basis, tmp_funct, T, P, tmp_err_r )
#                      tmp_dft_label, tmp_dft_obj = get_gms_object( tmp_basis, tmp_funct, T, P, tmp_dft_r )
#                      print( tmp_err_obj )
#                      print( tmp_dft_obj )
#                      break

 
if __name__ == '__main__':
   main()
