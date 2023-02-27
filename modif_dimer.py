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
print(user)
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

tp_dict = {
            '0'  : { 'T' :   '5', 'P' :   '0' },
            '1'  : { 'T' :  '90', 'P' :   '0' },
            '2'  : { 'T' :  '90', 'P' :  '90' },
            '3'  : { 'T' :  '90', 'P' : '180' },
            '4'  : { 'T' :  '90', 'P' : '270' },
            '5'  : { 'T' : '175', 'P' :   '0' },
            '6'  : { 'T' :  '45', 'P' :   '0' },
            '7'  : { 'T' :  '45', 'P' :  '90' },
            '8'  : { 'T' :  '45', 'P' : '180' },
            '9'  : { 'T' :  '45', 'P' : '270' },
      }

global VERBOSE
VERBOSE = False
global R_SCAN_LIST
global R_EQUIL_LIST

global CAT_LABEL
global ANI_LABEL
global DIMER_LABEL
global HCER
DIMER_LABEL = 'EMIM_BF4' #sys.argv[1]   
CAT_LABEL, ANI_LABEL = DIMER_LABEL.split('_')
CAT_LABEL, ANI_LABEL = 'EMIM', 'BF4'
 
for rm_basis in ['PCseg-2', 'APCseg-2']:
  try:
    gbasis_list.remove(rm_basis)
  except(ValueError):
    pass

R_SCAN_LIST  = full_R_list 
equil_T_list = full_T_list
equil_P_list = full_P_list

R_SCAN_LIST  = [ #'2.0', '2.1', '2.2', '2.3', '2.4', 
                 '2.5', '2.6', '2.7', '2.8', '2.9',
                 '3.0', '3.1', '3.2', '3.3', '3.4', 
                 '3.5', '3.6', '3.7', '3.8', '3.9',
                 '4.0', '4.1', '4.2', '4.3', '4.4', 
                 '4.5', '4.6', '4.7', '4.8', '4.9', 
                 '5.0', '5.5', '6.0', '6.5', '7.0', 
                 '7.5', 
                 #'8.0', '8.5', '9.0', '9.5', 
                 #'10.0', '11.0', '12.0', '13.0', '15.0' 
                 ]
reduced_R_list = [ 
                   '2.5', '2.8', 
                   '3.0', '3.3', '3.6', '3.9', 
                   '4.2', '4.5', '4.8', 
                   '5.1'
#                  '5.5', '6.0', '6.5', '7.0' 
                    ]

T_list = [ '5', '90', '175' ]
P_list = [ '0', '90', '180', '270' ] # ['90'] 
T_list = [ '90' ]
P_list = [ '90' ] 

global READ_FROM

READ_FROM = sys.argv[1]
print( READ_FROM )

#gbasis_list      = [ 'APCseg-1' ] 
gbasis_list      = [ 'APCseg-1', 'STO',  'N311' ]
functionals_list = [ 'PBE0', 'B3LYP', 'M11' , 'wB97x-D' ]
#df_columns = [('Coordinates','Radius'), ('Coordinates','Theta'), ('Coordinates', 'Phi'), ('LAMMPS','INT.EN.') ]
#df = pd.DataFrame( columns = pd.MultiIndex.from_tuples( df_columns ), dtype=object )

common_columns = [ 'Radius', 'Theta', 'Phi', 'TOT.EN.', 'INT.EN.', 'CHARG.ANI.', 'CHARG.CAT.', 'COUL.EN.',
                   'BASIS.DIM.', 'INERT.MOM', 'COM', 'COC', 'Run.Time']

def write_pd_series( R, T, P, scan_out_dict, scan_inp_dict, post_proc=False, equil=False ):

    pd_dict = { 'Radius' : R, 'Theta' : T, 'Phi' : P, 'BASIS.DIM.': scan_out_dict['BASIS.DIM.'], 
                'INERT.MOM' : scan_out_dict['INERT.MOM.'], 'Run.Time'  : scan_out_dict['TIME'] } 
    if post_proc == 'MP2':
      cart_dict = None
      pd_dict[ 'MP2.EN.'] = scan_out_dict['MP2']['MP2.EN.'] - ZERO_MP2_ENER
    elif post_proc == 'EDA':
      cart_dict = None
      print( write_out )
    else:
      if equil:
        pd_dict['Relax.Radius'] = scan_out_dict['FINAL']['ZMAT'][19]['STR']['val']
      mull_charges = scan_out_dict['MULL.CHARGES'] 
      #mull_charges = scan_out_dict['FINAL']['MULL.CHARGES'] 
      cart_coords  = scan_out_dict['FINAL']['CART.COORDS.']
      com = center_of_charge( cart_coords, mull_charges ) 
      coc = center_of_mass( cart_coords ) 
      cat_coords = dict(( str(k), cart_coords[k]) for k in range(CAT_NAT) )
      ani_coords = dict(( str(k), cart_coords[k]) for k in range(CAT_NAT, CAT_NAT + ANI_NAT ) )
      cat_com = center_of_mass( cat_coords )
      ani_com = center_of_mass( ani_coords )
      dcom = LA.norm( np.array(cat_com)-np.array(ani_com ))
      pd_dict[ 'COM' ] = com
      pd_dict[ 'COC' ] = coc
      pd_dict[ 'DIST.COM' ] = dcom
      cart_dict = { 'Radius' : R, 'cart.coords.' : cart_coords, 'mull.charges' : mull_charges }
      pd_dict[ 'CHARG.CAT.' ] = scan_out_dict['CHARG.CAT.']
      pd_dict[ 'CHARG.ANI.' ] = scan_out_dict['CHARG.ANI.'] 
      pd_dict['COUL.EN.'] = Coulomb_Energy(float(pd_dict['DIST.COM']), float(pd_dict['CHARG.CAT.']), float(pd_dict['CHARG.ANI.']))

    pd_dict[ 'TOT.EN.'] = scan_out_dict['TOT.EN.']
    pd_dict[ 'INT.EN.'] = scan_out_dict['INT.EN.']


    pd_series = pd.Series( pd_dict,  dtype=object )
    cart_series = pd.Series( cart_dict, dtype=object )
    return( pd_series, cart_series )

def get_gms_object( basis, funct, T, P, R, equil = False, opt_method = 'QA', post_scf = 'DFTTYP', run_type = 'OPTIMIZE', 
                    read_from='ISOLATED' ):
    if equil:
       TPR_CONF = IL.DIMER_EQUIL_CONF( DIMER_LABEL, basis, funct, T=T, P=P, R=R )
       TPR_label = 'EQUIL_{}_T_{}_P_{}_R_{}_{}_{}'.format(DIMER_LABEL.lower(), T, P, R, basis, funct )
       tmp_ifreeze = '53,54'
    else:
       TPR_CONF = IL.DIMER_SCAN_CONF( DIMER_LABEL, basis, funct, T=T, P=P, R=R )
       TPR_label = 'SCAN_{}_T_{}_P_{}_R_{}_{}_{}'.format(DIMER_LABEL.lower(), T, P, R, basis, funct )
       tmp_ifreeze = '52,53,54'
 
    #if full_relax:
#    print(1, TPR_label )
    TPR_CONF.R_dir = TPR_CONF.R_dir.replace('SCAN','SCAN_from_{}'.format(read_from))
    TPR_label = TPR_label.replace('SCAN','SCAN_from_{}'.format(read_from))
#    print(2, TPR_label )

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
    os.makedirs( DIMER.csv_dir, exist_ok=True)

    CAT_NAT = DIMER.cat_dict['nat']
    ANI_NAT = DIMER.ani_dict['nat']

    CATION = IL.MONOMER( CAT_LABEL, basis, funct )
    ANION  = IL.MONOMER( ANI_LABEL, basis, funct )
      
    cat_zmat   = CATION.mono_dict['OUT'][basis][funct]['DFT']['ZMAT']
    CAT_DFT_EN = CATION.mono_dict['OUT'][basis][funct]['DFT']['TOT.EN.'] 

    ani_zmat   = ANION.mono_dict['OUT'][basis][funct]['DFT']['ZMAT']
    ANI_DFT_EN = ANION.mono_dict['OUT'][basis][funct]['DFT']['TOT.EN.']

    ZERO_DFT_ENER = CAT_DFT_EN + ANI_DFT_EN

    try:
      ANI_MP2_EN = ANION.mono_dict['OUT'][basis][funct]['MP2']['MP2.EN.']
      CAT_MP2_EN = CATION.mono_dict['OUT'][basis][funct]['MP2']['MP2.EN.']
      ZERO_MP2_ENER = CAT_MP2_EN + ANI_MP2_EN
    except(KeyError):
      print_tab( 3, 'Missing (MP2) monomer' )
      cat_zmat, ani_zmat, ZERO_DFT_ENER, ZERO_MP2_ENER = False, False, False, False
    
    return( DIMER, CAT_NAT, ANI_NAT, cat_zmat, ani_zmat, ZERO_DFT_ENER, ZERO_MP2_ENER )

def post_process( pp_label, pp_obj, dft_zmat, R, T, P, pp_jq = 'nodesloq' ):
    pp_running = running_label( pp_obj.inp_name )
    pp_series = None
    if pp_running:
       print_tab( 4, 'Running ({})'.format(pp_label) )
    else:
       if os.path.exists( pp_obj.inp_file ):
          print_tab( 4, '{} ok (new)'.format(pp_label) )
          pp_exec, pp_exec_err = pp_obj.get_job_exec()
          pp_inp_dict, pp_out_dict, pp_scf, ppn_geom = pp_obj.get_job_results()
          pp_series = write_pd_series( R, T, P, pp_out_dict, pp_inp_dict, post_proc = pp_label )[0]
       else:
          print_tab( 4, 'Submitting {}'.format(pp_label) )
          pp_obj.run_new( zmat_dict = dft_zmat, job_queue = pp_jq )

    return pp_series

#def read_all_csv( DIMER, equil = False ):
#
#    if equil:
#       ## read CSV files for EQUIL 
#       if os.path.exists( DIMER.equil_dft_csv ):
#          eq_dft_df = pd.read_csv( DIMER.equil_dft_csv, index_col=0, dtype=object )
#       else:
#          eq_dft_df = pd.DataFrame( columns = common_columns + [ 'Relax.Radius' ] )
#   
#       if os.path.exists( DIMER.equil_mp2_csv ):
#          eq_mp2_df = pd.read_csv( DIMER.equil_mp2_csv, index_col=0, dtype=object )
#       else:
#          eq_mp2_df = pd.DataFrame( columns = common_columns ) 
#   
#       return_list = [eq_dft_df, eq_mp2_df]
#       ## read CSV files for EQUIL 
#
#    else:
#       ## read CSV files for SCAN 
#       if os.path.exists( DIMER.scan_dft_csv ):
#          dft_df = pd.read_csv( DIMER.scan_dft_csv, index_col=0, dtype=object )
#       else:
#          dft_columns = common_columns 
#          dft_df  = pd.DataFrame( columns = dft_columns )
#       
#       if os.path.exists( DIMER.scan_mp2_csv ):
#          mp2_df = pd.read_csv( DIMER.scan_mp2_csv, index_col=0, dtype=object )
#       else:
#          mp2_df  = pd.DataFrame( columns = common_columns + ['MP2.EN.'] )
#       
#       if os.path.exists( DIMER.scan_eda_csv ):
#          eda_df = pd.read_csv( DIMER.scan_eda_csv, index_col=0, dtype=object )
#       else:
#          eda_df  = pd.DataFrame( columns = ['Radius', 'Theta', 'Phi'] )
#       
#       if os.path.exists( DIMER.scan_coords_csv ):
#          cart_df = pd.read_csv( DIMER.scan_coords_csv, index_col=0, dtype=object )
#       else:
#          cart_df = pd.DataFrame(columns = [ 'Radius', 'cart.coords.', 'mull.charges' ] ) 
#
#       return_list = [dft_df, mp2_df, eda_df, cart_df]
#       ## read CSV files for SCAN
#
#    return return_list
   
def main():

    global CAT_NAT
    global ANI_NAT
    global ZERO_DFT_ENER
    global ZERO_MP2_ENER
    ## MAKE EQUIL (fixed T, P)
    print_tab( 0, '>>>> {} <<<<'.format(DIMER_LABEL) )
    for tmp_basis in gbasis_list: #['STO', 'N311']: 
        for tmp_funct in functionals_list: #['PBE0', 'B3LYP']:
           print_tab( 1, '=== {} ===='.format(tmp_basis) )
           print_tab( 2, '== {} ==='.format(tmp_funct) )
    
           DIMER, CAT_NAT, ANI_NAT, cat_zmat, ani_zmat, ZERO_DFT_ENER, ZERO_MP2_ENER = read_dimer( tmp_basis, tmp_funct )

           if [cat_zmat, ani_zmat, ZERO_DFT_ENER, ZERO_MP2_ENER] == [False, False, False, False]:
              proceed = False
           else:
              proceed = True

           if proceed:

              ## read CSV files for EQUIL 
              if os.path.exists( DIMER.equil_dft_csv ):
                 eq_dft_df = pd.read_csv( DIMER.equil_dft_csv, index_col=0, dtype=object )
              else:
                 eq_dft_df = pd.DataFrame( columns = common_columns + [ 'Relax.Radius' ] )
   
              if os.path.exists( DIMER.equil_mp2_csv ):
                 eq_mp2_df = pd.read_csv( DIMER.equil_mp2_csv, index_col=0, dtype=object )
              else:
                 eq_mp2_df = pd.DataFrame( columns = common_columns ) 
   

              for T in T_list:
                for P in P_list:

                  run_equil = False
                  if READ_FROM == 'HCER':
                     run_equil = True
                     HCER = False
                  ##########################################################################################################
                  if run_equil:
                     print_tab( 2, 'EQUILIBRIUM starts' )
                     find_lcer = True ## stop when lowest converged equilibrium radius reached
                     for R in ['10.0', '9.0', '8.0', '7.0']: # R_EQUIL_LIST:
                       if find_lcer:
                          dft_line = eq_dft_df.loc[ eq_dft_df['Radius']==str(R) ].loc[ 
                                                    eq_dft_df[ 'Theta']==str(T) ].loc[ eq_dft_df['Phi']==str(P) ]
                          mp2_line = eq_mp2_df.loc[ eq_mp2_df['Radius']==str(R) ].loc[ 
                                                    eq_mp2_df[ 'Theta']==str(T) ].loc[ eq_mp2_df['Phi']==str(P) ]
                          if dft_line.empty or mp2_line.empty:
                             run_mp2 = False
                             run_hes = False
                             print_tab( 3, 'T = {}, P = {}, R = {}'.format(T,P,R) )
                             TPR_eq_label, eq_obj = get_gms_object( tmp_basis, tmp_funct, T, P, R, equil = True )
                             running = running_label( eq_obj.inp_name )
                             if running:
                                print_tab( 4, 'Running' )
                                find_lcer = False ## will not go below this radius
                             else:
                                if os.path.exists( eq_obj.inp_file ):
                                   eq_exec, eq_exec_err = eq_obj.get_job_exec()
                                   if [eq_exec, eq_exec_err] == ['TERMINATED.NORMALLY', False]:
                                      eq_inp_dict, eq_out_dict, eq_scf, eq_geom = eq_obj.get_job_results()
                                      if [ eq_scf, eq_geom ] == ['CONVERGED', 'LOCATED']:
                                         print_tab( 4, 'OPT.EQ. ok' )
                                         run_mp2 = True
                                         run_hes = True
                                         if dft_line.empty:
                                            eq_series = write_pd_series( R, T, P, eq_out_dict, eq_inp_dict, equil=True )[0]
                                            eq_dft_df = eq_dft_df.append( eq_series , ignore_index=True )
                                         find_lcer = False ## will not go below this radius
                                         HCER = R
                                      else:
                                         print_tab( 4, 'OPT.EQ. not ok' )
                                         eq_obj.fix_error()
                                   else:
                                      print_tab( 4, 'OPT.EQ. FAILED' )
                                      eq_obj.fix_error()
                                else:
                                   os.makedirs( eq_obj.run_dir, exist_ok=True )
                                   comp_zmat = compose_zmatrices( cat_zmat, ani_zmat, radius=R , theta=T, phi=P )
                                   eq_obj.run_new( zmat_dict = comp_zmat, msg='equilibrium', job_queue='nodeshiq' )
                                   #eq_obj.run_new( zmat_dict = comp_zmat, msg='equilibrium', job_queue='nodeshiq' )
                                   find_lcer = False ## will skip further radii

                             ## MP2 STARTS
                             if run_mp2:
                                eq_mp2_label, eq_mp2_obj = get_gms_object( tmp_basis, tmp_funct, T, P, R, equil = True, 
                                                                           post_scf = 'MP2', run_type = 'ENERGY' )
                                eq_dft_zmat = eq_out_dict['FINAL']['ZMAT']
                                eq_mp2_series = post_process( 'MP2', eq_mp2_obj, eq_dft_zmat, R, T, P )
                                if isinstance(eq_mp2_series, pd.core.series.Series) :
                                   eq_mp2_df = eq_mp2_df.append( eq_mp2_series, ignore_index=True )
                             ## MP2 ENDS

                             ## HES STARTS
                             run_hes = False
                             if run_hes:
                                eq_hes_label, eq_hes_obj = get_gms_object( tmp_basis, tmp_funct, T, P, R, equil = True, 
                                                                           post_scf = 'DFTTYP', run_type = 'HESSIAN' )
                                eq_dft_zmat = eq_out_dict['FINAL']['ZMAT']
                                eq_hes_series = post_process( 'HES', eq_hes_obj, eq_dft_zmat, R, T, P )
                                if isinstance(eq_hes_series, pd.core.series.Series) :
                                   eq_hes_df = eq_hes_df.append( eq_hes_series, ignore_index=True )
                             ## HES ENDS
                     print_tab( 3, '--- Highest converged equilibrium radius = {}'.format(HCER) )
                     print_tab( 2, 'EQUILIBRIUM ends' )

                  ##########################################################################################################

                  ## find equilibrium with highest initial distance
                  new_csv_dir = os.path.join( DIMER.csv_dir, 'READ_FROM_{}'.format(READ_FROM) )
                  if READ_FROM == 'HCER':
                     if HCER:
                        print_tab( 3, '--- Reading zmat from equilibrium radius = {}'.format(HCER) )
                        runs_dir = '/data/mdi0316/WORK/DIMERS/EMIM_BF4/RUNS/'
                        read_HCER_dir = os.path.join(  runs_dir, 'SCAN_from_{}'.format(HCER), tmp_basis, tmp_funct )
                        write_HCER_dir = os.path.join( runs_dir, 'SCAN_from_HCER', tmp_basis, tmp_funct )
                        if not os.path.exists( write_HCER_dir ):
                              sp.call( 'mkdir -p {}'.format( write_HCER_dir ), shell=True )
                        sp.call( 'rsync -r {}/ {} --update'.format( read_HCER_dir, write_HCER_dir ), shell=True )
                        run_scan = True
                        old_radius = HCER 
                  elif READ_FROM == 'ISOLATED':
                        run_scan = True
                        old_radius = 'ISOLATED' 

                  os.makedirs( new_csv_dir, exist_ok=True)
                  DIMER.scan_dft_csv = os.path.join( new_csv_dir, 'scan_dft.csv' )
                  DIMER.scan_mp2_csv = os.path.join( new_csv_dir, 'scan_mp2.csv' )
                  DIMER.scan_coords_csv = os.path.join( new_csv_dir, 'scan_coords.csv' )


                  ## read CSV files for SCAN 
                  if os.path.exists( DIMER.scan_dft_csv ):
                     dft_df = pd.read_csv( DIMER.scan_dft_csv, index_col=0, dtype=object )
                  else:
                     dft_columns = common_columns 
                     dft_df  = pd.DataFrame( columns = dft_columns )
                  
                  if os.path.exists( DIMER.scan_mp2_csv ):
                     mp2_df = pd.read_csv( DIMER.scan_mp2_csv, index_col=0, dtype=object )
                  else:
                     mp2_df  = pd.DataFrame( columns = common_columns + ['MP2.EN.'] )
                  
                  if os.path.exists( DIMER.scan_eda_csv ):
                     eda_df = pd.read_csv( DIMER.scan_eda_csv, index_col=0, dtype=object )
                  else:
                     eda_df  = pd.DataFrame( columns = ['Radius', 'Theta', 'Phi'] )
                  
                  if os.path.exists( DIMER.scan_coords_csv ):
                     cart_df = pd.read_csv( DIMER.scan_coords_csv, index_col=0, dtype=object )
                  else:
                     cart_df = pd.DataFrame(columns = [ 'Radius', 'cart.coords.', 'mull.charges' ] ) 
                  ## read CSV files for SCAN 

                  if run_scan: 
                     print_tab( 2, 'SCAN starts' )
                     for R in R_SCAN_LIST:
                     #for R in loop_R_list:
                         print_tab( 3, '   R = {}'.format(R) )
                         run_dft = False
                         run_mp2 = False
                         run_eda = False
                         dft_line = dft_df.loc[dft_df['Radius']==str(R)].loc[dft_df['Theta']==str(T)].loc[dft_df['Phi']==str(P)]
                         mp2_line = mp2_df.loc[mp2_df['Radius']==str(R)].loc[mp2_df['Theta']==str(T)].loc[mp2_df['Phi']==str(P)]
                         eda_line = eda_df.loc[eda_df['Radius']==str(R)].loc[eda_df['Theta']==str(T)].loc[eda_df['Phi']==str(P)]
                         if dft_line.empty or mp2_line.empty: 
                         #if dft_line.empty or mp2_line.empty or eda_line.empty: 
                            #TPR_CONF = IL.DIMER_SCAN_CONF( DIMER_LABEL, tmp_basis, tmp_funct, T=T, P=P, R=R )
                            #print(TPR_CONF) #= IL.DIMER_SCAN_CONF( DIMER_LABEL, tmp_basis, tmp_funct, T=T, P=P, R=R )
                            TPR_scan_label, scan_obj = get_gms_object( tmp_basis, tmp_funct, T, P, R, read_from=old_radius )

                            #old_label = 'SCAN_{}_T_{}_P_{}_R_{}_{}_{}'.format(DIMER_LABEL.lower(), T, P, R, tmp_basis, tmp_funct )
                            #old_root_dir = os.path.join( work_dir, 'DIMERS/EMIM_BF4/RUNS/SCAN_from_{}'.format(HCER), 
                            #                             tmp_basis, tmp_funct, 'T_{}'.format(T), 'P_{}'.format(P), 'R_{}'.format(R) ) 
                            #old_obj = GAMESS.GAMESS( inp_label = old_label, root_dir = old_root_dir, 
                            #                         natoms = CAT_NAT + ANI_NAT, nat_cat = CAT_NAT, nat_ani = ANI_NAT, 
                            #                         icharge = 0, zero_energy = ZERO_DFT_ENER,
                            #                         run_type = 'OPTIMIZE', post_scf = 'DFTTYP',
                            #                         basis = tmp_basis, functional = tmp_funct,
                            #                         ifreeze = '52,53,54', opt_method = 'QA' )
                            #scan_obj = old_obj
                            if VERBOSE:
                               print_tab( 3, '{}, {}'.format(scan_obj.run_dir, scan_obj.inp_name ) )
                            running = running_label( scan_obj.inp_name )
                            if running:
                               print_tab( 4, 'Running' )
                            else:
                               if os.path.exists( scan_obj.inp_file ):
                                  scan_exec, scan_exec_err = scan_obj.get_job_exec()
                                  if [ scan_exec, scan_exec_err ] == ['TERMINATED.NORMALLY', False]:
                                     scan_inp_dict, scan_out_dict, scan_scf, scan_geom = scan_obj.get_job_results()
                                     if [ scan_scf, scan_geom ] == ['CONVERGED', 'LOCATED']:
                                        print_tab( 4, 'OPT.SCAN ok (new)' )
                                        run_mp2 = True
                                        run_eda = True
                                        if dft_line.empty:
                                           dft_series, cart_series = write_pd_series( R, T, P, scan_out_dict, scan_inp_dict )
                                           dft_df  =  dft_df.append( dft_series , ignore_index=True )
                                           cart_df = cart_df.append( cart_series , ignore_index=True )
                                     else:
                                        print_tab( 4, 'OPT.SCAN not ok' )
                                        scan_obj.fix_error()
                                  else:
                                     scan_obj.fix_error()
                               else:
                                 run_dft = True
                         else:
                            print_tab( 4, 'DFT + MP2 = OK' )
                            if VERBOSE:
                               print_tab( 4, dft_df.loc[ dft_df['Radius'] == str(R) ] )

                         ## RUN NEW DFT
                         if run_dft:
                            if READ_FROM == 'ISOLATED':
                               guess_zmat = compose_zmatrices( cat_zmat, ani_zmat, radius=R , theta=T, phi=P )
                               msg = 'zmat.from.isolated.ions'
                            else:
                               guess_label, guess_obj = get_gms_object( tmp_basis, tmp_funct, T, P, HCER, equil = True)
                               guess_inp_dict, guess_out_dict, guess_scf, guess_geom = guess_obj.get_job_results()
                               guess_zmat = guess_out_dict['FINAL']['ZMAT']
                               guess_zmat[19]['STR']['val'] = R
                               msg = 'zmat.from.equil.R={}'.format(HCER)
                               print( scan_obj )
                               print( guess_obj )
                            
                            scan_obj.run_new( zmat_dict = guess_zmat, msg=msg, job_queue='nodesloq' )
                            #scan_obj.run_new( zmat_dict = guess_zmat, msg=msg, job_queue='nodeshiq' )
                            exit()

                         ## MP2 STARTS
                         run_mp2 = False
                         if mp2_line.empty:
                            if run_mp2:
                               mp2_label, mp2_obj = get_gms_object(tmp_basis,tmp_funct,T,P,R,post_scf='MP2',run_type='ENERGY')
                               dft_zmat = scan_out_dict['FINAL']['ZMAT']
                               mp2_series = post_process( 'MP2', mp2_obj, dft_zmat, R, T, P )
                               if isinstance(mp2_series, pd.core.series.Series) :
                                  mp2_df = mp2_df.append(mp2_series,ignore_index=True)
#                                  run_eda = True
                         ## MP2 ENDS

                         ## EDA STARTS
                         run_eda = False
                         if eda_line.empty:
                            if run_eda:
                               eda_label, eda_obj = get_gms_object(tmp_basis,tmp_funct,T,P,R,post_scf='DFTTYP',run_type='EDA')
                               dft_zmat = scan_out_dict['FINAL']['ZMAT']
                               dft_com  = dft_df.loc[ dft_df['Radius'] == R ]['DIST.COM'].values
                               eda_series = post_process( 'EDA', eda_obj, dft_zmat, R, T, P )
                               if isinstance(eda_series, pd.core.series.Series) :
                                  eda_df = eda_df.append(eda_series,ignore_index=True)
                               exit()
                         ## MP2 ENDS
                     print_tab( 2, 'SCAN ends' )

             
                print_tab( 2, ['--- Writing output to {}\n'.format(DIMER.csv_dir) ] )
                eq_dft_df.to_csv( DIMER.equil_dft_csv )
                eq_mp2_df.to_csv( DIMER.equil_mp2_csv )

                dft_df.sort_values(by=['Radius'], inplace=True)
                dft_df.to_csv( DIMER.scan_dft_csv )
                mp2_df.sort_values(by=['Radius'], inplace=True)
                mp2_df.to_csv( DIMER.scan_mp2_csv )
                cart_df.sort_values(by=['Radius'], inplace=True)
                cart_df.to_csv( DIMER.scan_coords_csv )

 
if __name__ == '__main__':
   main()
