#!/home/mdi0316/anaconda3/bin/python

### common input start
import os, sys, re
import shutil
import numpy as np
from numpy import linalg as LA
import pandas as pd
import subprocess as sp
import csv 
import time

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

from Functions import print_tab, running_jobs, compose_zmatrices, running_label, center_of_charge, center_of_mass, Coulomb_Energy, angle_between
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

fix_errors = False
fix_errors = True

R_EQUIL_LIST = [ '3.0', '4.0', '5.0', '6.0', '7.0', '8.0', '9.0', '10.0' ]
R_EQUIL_LIST = [ '10.0', '9.0' ] #, '8.0', '7.0', '6.0', '5.0' ]

global T_SCAN_LIST

full_gbasis_list = ['N311']
full_functionals_list = ['B3LYP']

T_SCAN_LIST = [ '90' ]
P_SCAN_LIST = [ '90' ]
R_SCAN_LIST = [ 
                '2.0', '2.5', #'3.0', '3.5',
#                '4.0', '4.5', '5.0', '5.5', 
#                '6.0', '6.5', '7.0', '7.5', 
#                '8.0', '8.5', '9.0', '9.5', 
#                '10.0', '11.0', '12.0', '13.0', '15.0' 
                ]

#if DIMER_LABEL == 'EMIM_BF4':
#   full_gbasis_list      = [ 'APCseg-1', 'STO',  'N311' ]
#   full_functionals_list = [ 'PBE0', 'B3LYP', 'M11' , 'wB97x-D' ]

print( 'full_gbasis_list: {}'.format( full_gbasis_list ))
print( 'full_functionals_list: {}'.format( full_functionals_list ))

print_time = True

RTP_columns  = ['Radius', 'Theta', 'Phi']
#'INERT.MOM', 'COM', 'COC', 'DIST.COM'

dft_columns = RTP_columns + ['TOT.EN.', 'INT.EN.', 'DISP.EN.', 'CHARG.ANI.', 'CHARG.CAT.', 'COUL.EN.', 'DIST.COM', 'BASIS.DIM.', 'Run.Time']
mp2_columns = RTP_columns + ['TOT.EN.', 'MP2.EN.', 'DISP.EN.', 'CHARG.ANI.', 'CHARG.CAT.', 'COUL.EN.', 'DIST.COM' ]
ccs_columns = RTP_columns + ['TOT.EN.', 'CCSD(T).EN.', 'DISP.EN.', 'CHARG.ANI.', 'CHARG.CAT.', 'COUL.EN.', 'DIST.COM']
err_columns = RTP_columns + ['ERROR']
frc_columns = RTP_columns + ['f_x', 'f_y', 'f_z']

coarse_grain_columns_0 = ['C0', 'C1', 'C2', 'A']
coarse_grain_columns_1 = ['x', 'y', 'z']
coarse_grain_columns = pd.MultiIndex.from_product([ coarse_grain_columns_0, coarse_grain_columns_1 ])

def print_converged_zmat( gms_obj, gms_out_dict ):
    ## converged zmat
    if 'FINAL' in gms_out_dict.keys():
       zmat_dict = gms_out_dict['FINAL']['ZMAT']
    else:
       zmat_dict = gms_out_dict['ZMAT']
    with open( gms_obj.zmat_file, 'w+') as f:
       w = csv.DictWriter(f, zmat_dict.keys())
       w.writeheader()
       w.writerow(zmat_dict)

def print_coords_charges( gms_obj, gms_out_dict ):
    ## atomic coordinates and charges
    acac_df = pd.DataFrame( columns=['elem.', 'idx.', 'x', 'y', 'z', 'pop.', 'charge'] )
    cart_coords = gms_out_dict['CART.COORDS.']
    mull_charges = gms_out_dict['MULL.CHARGES']
    for cc, mc in zip( cart_coords.values(), mull_charges.values() ):
        acac_df = acac_df.append( pd.Series( { 'elem.' : cc['elem.'], 'idx.' : cc['idx.'],
                                               'x': cc['x'], 'y': cc['y'], 'z': cc['z'],
                                               'pop.' : mc['pop.'], 'charge': mc['charge'] } ), ignore_index = True )
    acac_df.to_csv( gms_obj.acac_file )
 
def print_internucl_dist( gms_obj, gms_out_dict ):
    ## internuclear distances
    natoms = len(gms_out_dict['CART.COORDS.'])
    ind_df = pd.DataFrame(columns = [ 'idx1', 'idx2', 'elem.1', 'elem.2', 'distance' ] )
    for at_idx_1 in range(1, natoms+1):
        for at_idx_2 in range(1, at_idx_1):
            [ (tmp_k, tmp_v) ]= [ (k,v) for k,v in gms_out_dict['INTERNUCL.DISTANCES'].items() if
                                int(v['at.1']['idx.1']) == at_idx_1 and int(v['at.2']['idx.2']) == at_idx_2 ]
            tmp_dict = {  'idx1'     : tmp_v['at.1']['idx.1'], #'idx{}'.format(v['at.1']['idx.1']),
                          'idx2'     : tmp_v['at.2']['idx.2'], #'idx{}'.format(v['at.2']['idx.2']),
                          'elem.1'   : tmp_v['at.1']['elem.1'],
                          'elem.2'   : tmp_v['at.2']['elem.2'],
                          'distance' : tmp_v['dist.'] }
            ind_df = ind_df.append( pd.Series( tmp_dict ), ignore_index=True )
 
    ind_df.idx2 = pd.to_numeric( ind_df.idx2 )
    ind_df.idx1 = pd.to_numeric( ind_df.idx1 )
    ind_df = ind_df.sort_values('idx2')
    ind_df = ind_df.sort_values('idx1')
    ind_df = ind_df.reset_index( drop=True )
    ind_df.to_csv( gms_obj.intd_file )

def print_beads_coords( gms_obj, gms_out_dict ):

    ## center of mass
    cart_coords = gms_out_dict['CART.COORDS.']
    com_df = pd.DataFrame( [cart_coords] )
    com_df.to_csv( gms_obj.com_file )  

    ## write position of two more extern Carbons 
    cat_coords = dict(( str(k), cart_coords[k]) for k in range(CAT_NAT) )
    ani_coords = dict(( str(k), cart_coords[k]) for k in range(CAT_NAT, CAT_NAT + ANI_NAT ) )
    cat_com = center_of_mass( cat_coords )
    ani_com = center_of_mass( ani_coords )
    CC_dist_dict = { k:v for (k,v) in gms_out_dict['INTERNUCL.DISTANCES'].items() if v['at.1']['elem.1'] == 'C' and v['at.2']['elem.2'] == 'C' }
    CC_max_dist = max(CC_dist_dict, key=lambda k: CC_dist_dict[k]['dist.']) 

    at1, at2, dist = CC_dist_dict[CC_max_dist].values()
    elem1, idx1 = at1.values()
    elem2, idx2 = at2.values()
    C1_dict = list( [ vv for (kk, vv) in cart_coords.items() if vv['elem.'] == elem1 and vv['idx.'] == idx1 ] )[0] 
    C2_dict = list( [ vv for (kk, vv) in cart_coords.items() if vv['elem.'] == elem2 and vv['idx.'] == idx2 ] )[0]
 
    C1_xyz = np.array( [C1_dict['x'], C1_dict['y'], C1_dict['z']] ).astype(np.float)
    C2_xyz = np.array( [C2_dict['x'], C2_dict['y'], C2_dict['z']] ).astype(np.float)

    ### angles between C1,C2 and Center of mass
    C1_xyz -= cat_com 
    C2_xyz -= cat_com 
    C1_C2_dist = np.linalg.norm( C1_xyz - C2_xyz )

    #print( angle_between(C1_xyz, C2_xyz) * 180/math.pi )
    cosine_angle = np.dot(C1_xyz, C2_xyz) / (np.linalg.norm(C1_xyz) * np.linalg.norm(C2_xyz))
    C1_C2_angle = np.arccos(cosine_angle)

    beads_dict = {
                  ('A','x')  : ani_com[0],   ('A','y')  : ani_com[1],   ('A','z')  : ani_com[2],
                  ('C0','x') : cat_com[0],   ('C0','y') : cat_com[1],   ('C0','z') : cat_com[2],
                  ('C1','x') : C1_dict['x'], ('C1','y') : C1_dict['y'], ('C1','z') : C1_dict['z'],
                  ('C2','x') : C2_dict['x'], ('C2','y') : C2_dict['y'], ('C2','z') : C2_dict['z'],
                }

    beads_df = pd.DataFrame( [beads_dict] )
    beads_df.to_csv( gms_obj.beads_file, index=False)
    

def calculate_com( gms_obj, gms_out_dict ):

    cart_coords = gms_out_dict['CART.COORDS.']
    mull_charges = gms_out_dict['MULL.CHARGES'] 

    com = center_of_mass( cart_coords ) 
    coc = center_of_charge( cart_coords, mull_charges ) 
    cat_coords = dict(( str(k), cart_coords[k]) for k in range(CAT_NAT) )
    ani_coords = dict(( str(k), cart_coords[k]) for k in range(CAT_NAT, CAT_NAT + ANI_NAT ) )
    cat_com = center_of_mass( cat_coords )
    ani_com = center_of_mass( ani_coords )
    dcom = LA.norm( np.array(cat_com)-np.array(ani_com ))

    return( com, coc, cat_com, ani_com, dcom )


def print_dft_results( gms_obj, gms_out_dict ):

    com, coc, cat_com, ani_com, dcom = calculate_com( gms_obj, gms_out_dict )

    dft_dict = {}
    dft_dict['TOT.EN.'] = gms_out_dict['TOT.EN.']
    dft_dict['INT.EN.'] = gms_out_dict['INT.EN.']
    dft_dict['CHARG.CAT.'] = gms_out_dict['CHARG.CAT.']
    dft_dict['CHARG.ANI.'] = gms_out_dict['CHARG.ANI.'] 
    dft_dict['DIST.COM'] = dcom
    dft_dict['COUL.EN.'] = Coulomb_Energy( dcom,float(dft_dict['CHARG.CAT.']), float(dft_dict['CHARG.ANI.']))
    dft_dict['DISP.EN.'] = dft_dict['INT.EN.'] - dft_dict['COUL.EN.']
    dft_dict['BASIS.DIM.'] = gms_out_dict['BASIS.DIM.']
    dft_dict['Run.Time'] = gms_out_dict['TIME']

    dft_dict = { **RTP_line, **dft_dict }

    tmp_dft_df = pd.DataFrame( [dft_dict] )
    tmp_dft_df.to_csv( gms_obj.dft_file )

    print( 'PRINTING' )
    for ii in tmp_dft_df.columns:
       print(tmp_dft_df[ii]  )
    return( tmp_dft_df )
    

def print_mp2_results( gms_obj, gms_out_dict ):

    com, coc, cat_com, ani_com, dcom = calculate_com( gms_obj, gms_out_dict )

    mp2_dict = {}
    mp2_dict['TOT.EN.'] = gms_out_dict['TOT.EN.']
    mp2_dict['MP2.EN.'] = gms_out_dict['MP2']['MP2.EN.'] - ZERO_MP2_ENER
    mp2_dict['CHARG.CAT.'] = gms_out_dict['CHARG.CAT.']
    mp2_dict['CHARG.ANI.'] = gms_out_dict['CHARG.ANI.'] 
    mp2_dict['DIST.COM'] = dcom
    mp2_dict['COUL.EN.'] = Coulomb_Energy( dcom, float(mp2_dict['CHARG.CAT.']), float(mp2_dict['CHARG.ANI.']))
    mp2_dict['DISP.EN.'] = gms_out_dict['MP2']['MP2.EN.'] - mp2_dict['COUL.EN.']
    mp2_dict['Run.Time'] = gms_out_dict['TIME']
 
    mp2_dict = { **RTP_line, **mp2_dict }
 
    tmp_mp2_df = pd.DataFrame( [mp2_dict] )
    tmp_mp2_df.to_csv( gms_obj.mp2_file )

    return( tmp_mp2_df )


#def write_pd_series( R, T, P, gms_out_dict, scan_inp_dict, post_proc=False, equil=False ):
#    t0 = time.time()
#    mull_charges = gms_out_dict['MULL.CHARGES'] 
#    cart_coords  = gms_out_dict['CART.COORDS.']
#    cart_dict = { 'Radius' : R, 'Theta' : T, 'Phi' : P, 'cart.coords.' : cart_coords, 'mull.charges' : mull_charges }
#    ##
#    com = center_of_mass( cart_coords ) 
#    coc = center_of_charge( cart_coords, mull_charges ) 
#    cat_coords = dict(( str(k), cart_coords[k]) for k in range(CAT_NAT) )
#    ani_coords = dict(( str(k), cart_coords[k]) for k in range(CAT_NAT, CAT_NAT + ANI_NAT ) )
#    cat_com = center_of_mass( cat_coords )
#    ani_com = center_of_mass( ani_coords )
#    dcom = LA.norm( np.array(cat_com)-np.array(ani_com ))
#
#    pd_dict = { 'Radius' : R, 'Theta' : T, 'Phi' : P }
#
#    ## common
#    pd_dict['TOT.EN.'] = gms_out_dict['TOT.EN.']
#    pd_dict['CHARG.CAT.'] = gms_out_dict['CHARG.CAT.']
#    pd_dict['CHARG.ANI.'] = gms_out_dict['CHARG.ANI.'] 
#    pd_dict['DIST.COM'] = dcom
#    pd_dict['COUL.EN.'] = Coulomb_Energy(float(pd_dict['DIST.COM']), float(pd_dict['CHARG.CAT.']), float(pd_dict['CHARG.ANI.']))
#
#    if not post_proc:
#       pd_dict['INT.EN.'] = gms_out_dict['INT.EN.']
#       pd_dict['DISP.EN.'] = gms_out_dict['INT.EN.'] - pd_dict['COUL.EN.']
#       pd_dict['BASIS.DIM.'] = gms_out_dict['BASIS.DIM.']
#       pd_dict['INERT.MOM']  = gms_out_dict['INERT.MOM.']
#       pd_dict['COM'] = com
#       pd_dict['COC'] = coc
#    
#    if post_proc == 'MP2':
#       pd_dict[ 'MP2.EN.'] = gms_out_dict['MP2']['MP2.EN.'] - ZERO_MP2_ENER
#       pd_dict['DISP.EN.'] = pd_dict['MP2.EN.'] - pd_dict['COUL.EN.']
#
#    elif post_proc in ['EDA', 'HES']:
#      print( 'write_out' )
#    else:
#      if equil:
#        pd_dict['Relax.Radius'] = gms_out_dict['FINAL']['ZMAT'][19]['STR']['val']
#
#    pd_dict['Run.Time'] = gms_out_dict['TIME']
#
#    pd_series = pd.Series( pd_dict ) #, dtype=object )
#    cart_series = pd.Series( cart_dict ) #, dtype=object )
#    t1 = time.time()
#    print( 'write pd', t1-t0)
#    return( pd_series, cart_series )

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

    cat_zmat, ani_zmat = {}, {}
    ZERO_DFT_ENER, ZERO_MP2_ENER = 0, 0

    for ion_label, ion_zmat in zip( [CAT_LABEL, ANI_LABEL], [cat_zmat, ani_zmat] ):

       ION_IL = IL.MONOMER( ion_label, basis, funct )
       ION_DFT_EN = ION_IL.mono_dict['OUT'][basis][funct]['DFT']['TOT.EN.'] 
       ION_MP2_EN = ION_IL.mono_dict['OUT'][basis][funct]['MP2']['MP2.EN.'] 
       ZERO_DFT_ENER += ION_DFT_EN
       ZERO_MP2_ENER += ION_MP2_EN

       # read ION_GMS.conv_zmat_file
       ION_GMS = GAMESS.GAMESS( inp_label = ion_label, run_dir = ION_IL.opt_dir )
       ion_reader = csv.DictReader(open(ION_GMS.zmat_file) )
       for ion_r in ion_reader:
          for k,v in ion_r.items():
            ion_zmat[k] = v

    return( DIMER, CAT_NAT, ANI_NAT, cat_zmat, ani_zmat, ZERO_DFT_ENER, ZERO_MP2_ENER )

def post_process( pp_label, pp_obj, dft_zmat, R, T, P, pp_jq = 'nodeshiq' ):
    pp_df = None
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
                if pp_label == 'MP2':
                   print_coords_charges( pp_obj, pp_out_dict )
                   pp_df = print_mp2_results( pp_obj, pp_out_dict ) 
                else:
                   print( update_pp_series )
                   ## need to output pp_df instead of pp_series 
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

    return pp_df

def read_object( read_obj, read_label='DFT', dft_zmat=None, jq = 'nodesloq' ):
    read_df = None
    read_exec_err = False 
    read_run = running_label( read_obj.inp_name )

    if read_run:
       print_tab( 4, '{}: Running'.format(read_label) )
    else:
       ## check all csv files
       if read_label == 'DFT' and os.path.exists( read_obj.dft_file ):
          pass
       elif read_label == 'MP2' and os.path.exists( read_obj.mp2_file ):
          pass
       else:
          print_tab( 4, 'Read output' )
          abort_dir = os.path.join( read_obj.run_dir, 'FAILED', 'ABORTED' )
          if os.path.exists( abort_dir ):
             print_tab( 4, 'ABORTED' )
             read_exec = 'ABORTED'
          else:
             if os.path.exists( read_obj.inp_file ):
                t0 = time.time()
                read_exec, read_exec_err = read_obj.get_job_exec()
                print( read_exec, read_exec_err ) 
                t1 = time.time()
                print( t1-t0 )
                exit()
                if [ read_exec, read_exec_err ] == ['MISSING.OUTPUT.FILE', False]:
                   read_obj.resubmit()
                elif [ read_exec, read_exec_err ] == ['TERMINATED.NORMALLY', False]:
                   read_inp_dict, read_out_dict, read_scf, read_geom = read_obj.get_job_results()
                   if read_scf == 'CONVERGED':
                      print_tab( 4, '{}: ok (new)'.format(read_label) )
                      print_coords_charges( read_obj, read_out_dict )
                      if read_label == 'DFT':
                         read_df = print_dft_results( read_obj, read_out_dict ) 
                         print_converged_zmat( read_obj, read_out_dict )
                         print_coords_charges( read_obj, read_out_dict )
                         print_internucl_dist( read_obj, read_out_dict )
                         print_beads_coords(   read_obj, read_out_dict )

                      elif read_label == 'MP2':
                         read_df = print_mp2_results( read_obj, read_out_dict ) 
                      else:
                         print( update_pp_series )
                         ## need to output pp_df instead of pp_series 
                         pp_series = write_pd_series( R, T, P, pp_out_dict, pp_inp_dict, post_proc = read_label )[0]
                   else:
                      print_tab( 4, '{}: not ok (new)'.format(read_label) )
                      read_obj.fix_error()
                else:
                   print_tab( 4, '{}: EXEC: {}, ERR: {}'.format(read_label, read_exec, read_exec_err ) )
                   read_obj.fix_error()
             else:
                print_tab( 4, 'Submitting {}'.format(read_label) )
                read_obj.run_new( zmat_dict = dft_zmat, job_queue = jq )

    return read_exec, read_exec_err, read_df



def main():

    global DIMER_LABEL
    global CAT_LABEL
    global ANI_LABEL
    global T_SCAN_LIST
    global P_SCAN_LIST
    global R_SCAN_LIST
    global CAT_NAT
    global ANI_NAT
    global ZERO_DFT_ENER
    global ZERO_MP2_ENER

    if len(sys.argv) == 1:
      DIMER_LABEL = 'EMIM_BF4'
    else:
      DIMER_LABEL = sys.argv[1]

    CAT_LABEL, ANI_LABEL = DIMER_LABEL.split('_')
 
    print_tab( 0, '>>>> {} <<<<'.format(DIMER_LABEL) )

    for tmp_basis in full_gbasis_list:
      for tmp_funct in full_functionals_list:
         print_tab( 1, '=== {} ==='.format(tmp_basis) )
         print_tab( 2, '=== {} ==='.format(tmp_funct) )

         if DIMER_LABEL in ['EMIM_BF4', 'EMIM_PF6']:
            if [ tmp_basis, tmp_funct ] == [ 'N311', 'B3LYP' ]:
               T_SCAN_LIST = [ '5', '45', '90', '135', '175' ]
               P_SCAN_LIST = [ '0', '45', '90', '135', '180', '225', '270', '315' ]
               R_SCAN_LIST += [ '2.6', '2.7', '2.8', '2.9',
                                '3.1', '3.2', '3.3', '3.4', 
                                '3.6', '3.7', '3.8', '3.9',
                                '4.1', '4.2', '4.3', '4.4', 
                                '4.6', '4.7', '4.8', '4.9' ]
               tmp_rad_list = [ float(r) for r in R_SCAN_LIST ]
               tmp_rad_list.sort()
               R_SCAN_LIST = [ str(r) for r in tmp_rad_list ]
            else:
              T_SCAN_LIST = [ '90' ]
              P_SCAN_LIST = [ '90' ]

         print_tab( 2, '  R_SCAN_LIST: {}'.format( R_SCAN_LIST ))
         print_tab( 2, '  T_SCAN_LIST: {}'.format( T_SCAN_LIST ))
         print_tab( 2, '  P_SCAN_LIST: {}'.format( P_SCAN_LIST ))

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
           crd_df = pd.DataFrame( columns = [ 'Radius', 'cart.coords.', 'mull.charges' ] ) 
           frc_df = pd.DataFrame( columns = frc_columns )
           crg_df = pd.DataFrame( columns = coarse_grain_columns )
              
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
           if os.path.exists( DIMER.scan_crd_csv ):
              crd_df = pd.read_csv( DIMER.scan_crd_csv, index_col=0, dtype=object )
           if os.path.exists( DIMER.scan_frc_csv ):
              frc_df = pd.read_csv( DIMER.scan_frc_csv, index_col=0, dtype=object )
           if os.path.exists( DIMER.scan_corase_grain_csv ):
              crg_df = pd.read_csv( DIMER.scan_corase_grain_csv, index_col=0, dtype=object )

#           ##############
#           if make_equil:
#              print_tab( 3, 'EQUILIBRIUM' )
#    
#              for T in equil_T_list:
#                for P in equil_P_list:
#                  for R in R_EQUIL_LIST:
#                    err_line = eq_err_df.loc[ eq_err_df['Radius']==float(R) ].loc[ eq_err_df[ 'Theta']==float(T) ].loc[ eq_err_df['Phi']==float(P) ]
#                    dft_line = eq_dft_df.loc[ eq_dft_df['Radius']==float(R) ].loc[ eq_dft_df[ 'Theta']==float(T) ].loc[ eq_dft_df['Phi']==float(P) ]
#                    mp2_line = eq_mp2_df.loc[ eq_mp2_df['Radius']==float(R) ].loc[ eq_mp2_df[ 'Theta']==float(T) ].loc[ eq_mp2_df['Phi']==float(P) ]
#                    if err_line.empty:
#                      if dft_line.empty or mp2_line.empty:
#                         run_mp2 = False
#                         run_hes = False
#                         print_tab( 3, 'T = {}, P = {}, R = {}'.format(T,P,R) )
#                         eq_label, eq_opt_obj = get_gms_object( tmp_basis, tmp_funct, T, P, R, equil = True )
#                         eq_label, eq_mp2_obj = get_gms_object( tmp_basis, tmp_funct, T, P, R, equil = True, post_scf = 'MP2', run_type = 'ENERGY')
#                         eq_label, eq_hes_obj = get_gms_object( tmp_basis, tmp_funct, T, P, R, equil = True, post_scf = 'DFTTYP', run_type = 'HESSIAN')
#                         running = running_label( eq_opt_obj.inp_name )
#                         if running:
#                            print_tab( 4, 'Running' )
#                            break ## skip till this is finished
#                         else:
#                            if os.path.exists( eq_opt_obj.inp_file ):
#                               eq_exec, eq_exec_err = eq_opt_obj.get_job_exec()
#                               if [ eq_exec, eq_exec_err ] == ['MISSING.OUTPUT.FILE', False]:
#                                  eq_obj.resubmit()
#                               if [eq_exec, eq_exec_err] == ['TERMINATED.NORMALLY', False]:
#                                  eq_inp_dict, eq_out_dict, eq_scf, eq_geom = eq_opt_obj.get_job_results()
#                                  if [ eq_scf, eq_geom ] == ['CONVERGED', 'LOCATED']:
#                                     print_tab( 4, 'OPT.EQ. ok' )
#                                     run_mp2 = True
#                                     if dft_line.empty:
#                                        eq_series = write_pd_series( R, T, P, eq_out_dict, eq_inp_dict, equil=True )[0]
#                                        eq_dft_df = eq_dft_df.append( eq_series , ignore_index=True )
#                                  else:
#                                     print_tab( 4, 'OPT.EQ. not ok' )
#                                     eq_opt_obj.fix_error()
#                                     eq_opt_err = eq_opt_obj.read_error()
#                               else:
#                                  print_tab( 4, 'OPT.EQ. FAILED' )
#                                  eq_opt_obj.fix_error()
#                            else:
#                               os.makedirs( eq_opt_obj.run_dir, exist_ok=True )
#                               comp_zmat = compose_zmatrices( cat_zmat, ani_zmat, radius=R , theta=T, phi=P )
#                               #eq_opt_obj.run_new( zmat_dict = comp_zmat, msg='equilibrium', job_queue='nodesloq' )
#                               eq_opt_obj.run_new( zmat_dict = comp_zmat, msg='equilibrium', job_queue='nodeshiq' )
#                               break ## skip till this is finished
#
#                            ## MP2 STARTS
#                            if run_mp2:
#                               eq_dft_zmat = eq_out_dict['FINAL']['ZMAT']
#                               eq_mp2_series = post_process( 'MP2', eq_mp2_obj, eq_dft_zmat, R, T, P )
#                               if isinstance(eq_mp2_series, pd.core.series.Series) :
#                                  eq_mp2_df = eq_mp2_df.append( eq_mp2_series, ignore_index=True )
#                            ## MP2 ENDS
#
#              eq_err_df.to_csv( DIMER.equil_err_csv )
#              eq_dft_df.to_csv( DIMER.equil_dft_csv )
#              eq_mp2_df.to_csv( DIMER.equil_mp2_csv )

           ##############
           if make_scan:
              print_tab(3, '==========================================')
              print_tab( 3, 'SCAN' )

              for T in T_SCAN_LIST:
                  print_tab(3, '==========================================')
                  for P in P_SCAN_LIST:
                      print_tab(3, '------------------------------------------')
                      for R in R_SCAN_LIST:
                        
                          #if tmp_basis == 'APCseg-1' and tmp_funct == 'PBE0' and R == '4.2':
                          #elif [T,P] == ['5', '0'] and float(R) in [2.0, 2.9]:
                          #elif [T,P] == ['5', '90'] and float(R) <= 3.8:
                          #elif [T,P] == ['5', '180'] and float(R) <= 3.1:
                          #elif [T,P] == ['45', '0'] and ( float(R) <= 3.8 or  4.3 <= float(R) <= 5.0 ) :
                          #elif [T,P] == ['45', '90'] and float(R) <= 2.6 :
                          #elif [T,P] == ['90', '0'] and  ( 2.9 <= float(R) <= 3.9 or  4.2 <= float(R) <= 4.6 ) :
                          #elif [T,P] == ['135', '0'] and float(R) <= 3.5 and float(R) != 3.2 :
                          #elif [T,P] == ['135', '90'] and float(R) in [6.5]:
                          #elif [T,P] == ['175', '0'] and float(R) <= 5.5:
                          #elif [T,P] == ['175', '90'] and float(R) <= 4.3:
                          #elif [T,P] == ['175', '180'] and float(R) <= 3.2:

                          RTP_time_start = time.time()
                          print_tab( 3, 'T = {}, P = {}, R = {}'.format(T,P,R) )
                          #print_tab( 3, 'T = {}, P = {}, R = {}'.format(type(T),type(P),type(R)) )
                          global RTP_line
                          RTP_line = { 'Radius' : R, 'Theta' : T, 'Phi' : P }
                          run_dft, run_mp2 = 2*[False]

                          err_line = err_df.loc[ err_df['Radius'] == R ].loc[ err_df['Theta'] == T ].loc[ err_df['Phi'] == P ]
                          dft_line = dft_df.loc[ dft_df['Radius'] == R ].loc[ dft_df['Theta'] == T ].loc[ dft_df['Phi'] == P ]
                          mp2_line = mp2_df.loc[ mp2_df['Radius'] == R ].loc[ mp2_df['Theta'] == T ].loc[ mp2_df['Phi'] == P ]
                         
                          if err_line.empty: 

                             if dft_line.empty or mp2_line.empty: 

                                scan_label, dft_obj = get_gms_object(tmp_basis,tmp_funct,T,P,R )
                                scan_label, mp2_obj = get_gms_object(tmp_basis,tmp_funct,T,P,R,post_scf='MP2',run_type='ENERGY')
     
                                dft_exec, dft_exec_err, dft_df = read_object( dft_obj ) 
                                if [ dft_exec, dft_exec_err ] == ['TERMINATED.NORMALLY', False]:
                                    print( dft_df )
                                    mp2_exec, mp2_exec_err, mp2_df = read_object( mp2_obj ) 
                                    exit()

                                abort_dir = os.path.join( mp2_obj.run_dir, 'FAILED', 'ABORTED' )
                                if os.path.exists( abort_dir ):
                                   print_tab( 4, 'ABORTED' )
                                else:
                                   if VERBOSE:
                                      print_tab( 3, '{}, {}'.format(mp2_obj.run_dir, mp2_obj.inp_name ) )
                                   #running = running_label( scan_label )
                                   run_check_scan = running_label( mp2_obj.inp_name )
                                   run_check_mp2  = running_label( mp2_obj.inp_name  )
                                   if run_check_scan:
                                      print_tab( 4, 'Running scan' )
                                   elif run_mp2: 
                                      print_check_tab( 4, 'Running mp2' )
                                   else:
                                      
                                      if os.path.exists( mp2_obj.inp_file ):

                                         running = False
                                         failed  = False

                                         if os.path.exists( mp2_obj.com_file   ) and \
                                            os.path.exists( mp2_obj.acac_file  ) and \
                                            os.path.exists( mp2_obj.intd_file  ) and \
                                            os.path.exists( mp2_obj.force_file ) and \
                                            os.path.exists( mp2_obj.zmat_file  ) and \
                                            os.path.exists( mp2_obj.beads_file ) and \
                                            os.path.exists( mp2_obj.dft_file   ):
                                            dft_df_line = pd.read_csv( mp2_obj.dft_file, index_col = 0 )
                                            print( 'READING' )
                                            for ii in dft_df_line.columns:
                                               print(dft_df_line[ii])

                                            exit()
                                         else:
                                            scan_exec, scan_exec_err = mp2_obj.get_job_exec()
                                            if [ scan_exec, scan_exec_err ] == ['MISSING.OUTPUT.FILE', False]:
                                               mp2_obj.resubmit()
                                               running = True
                                            elif [ scan_exec, scan_exec_err ] == ['TERMINATED.NORMALLY', False]:
                                               scan_inp_dict, scan_out_dict, scan_scf, scan_geom = mp2_obj.get_job_results()
                                               if [ scan_scf, scan_geom ] == ['CONVERGED', 'LOCATED']:
                                                  print_tab( 4, 'OPT: ok (new)' )
          
                                                  dft_df_line = print_dft_results( mp2_obj, scan_out_dict )

                                                  print_converged_zmat( mp2_obj, scan_out_dict )
                                                  print_coords_charges( mp2_obj, scan_out_dict )
                                                  print_internucl_dist( mp2_obj, scan_out_dict )
                                                  print_beads_coords( mp2_obj, scan_out_dict )

                                               else:
                                                  failed  = True
                                                  print_tab( 4, 'OPT: not ok' )
                                                  mp2_obj.fix_error()
                                                  mp2_obj_err = mp2_obj.read_error()
                                                  if mp2_obj_err in [ 'atoms.too.close', 'Stationary.Point.Location.failed', 'Gradient.out.of.range' ]:
                                                     print_tab( 4, 'OPT: EXEC: {}, ERR: {}'.format( scan_exec, mp2_obj_err ) )
                                                     failed_dict = { **RTP_line, **{ 'ERROR':mp2_obj_err } }
                                                     print( failed_dict )
                                                     failed_df   = pd.DataFrame( failed_dict )
                                                  else: 
                                                     print_tab( 4, 'Could not fix: {}'.format(mp2_obj_err) )

                                            else:
                                               failed  = True
                                               print_tab( 4, 'OPT: EXEC: {}, ERR: {}'.format( scan_exec, scan_exec_err ) )
                                               mp2_obj.fix_error()

                                         if not running: 
                                            if failed: 
                                               err_df = pd.concat( [ err_df, failed_df], ignore_index=True, sort=True ) 
                                            else:
                                               frc_new_line = pd.read_csv( mp2_obj.force_file, index_col = 0 )
                                               frc_new_line['Radius'], frc_new_line['Theta'], frc_new_line['Phi'] = R, T, P
                                               frc_df = pd.concat( [frc_df, frc_new_line], ignore_index=True, sort=True ) 

                                               dft_df = pd.concat( [dft_df, dft_df_line], ignore_index=True, sort=True ) 
                                               run_mp2 = True

                                      else:
                                        run_dft = True
                             else:
                                print_tab( 4, 'DFT, MP2 lines not empty ' )

                             ## RUN NEW DFT STARTS
                             if run_dft:
                                guess_zmat = compose_zmatrices( cat_zmat, ani_zmat, radius=R , theta=T, phi=P )
                                mp2_obj.run_new( zmat_dict = guess_zmat, msg = 'zmat.from.isolated.ions', job_queue='nodesloq' )
                             ## RUN NEW DFT ENDS

                             ## MP2 STARTS
                             if run_mp2 and mp2_line.empty and not run_check_mp2:
                                if os.path.exists( mp2_obj.mp2_file ) and os.path.exists( mp2_obj.acac_file ):
                                   mp2_df_line = pd.read_csv( mp2_obj.mp2_file, index_col = 0 )
                                else:
                                   scan_inp_dict, scan_out_dict, scan_scf, scan_geom = mp2_obj.get_job_results()
                                   dft_zmat = scan_out_dict['FINAL']['ZMAT']
                                   mp2_df_line = post_process( 'MP2', mp2_obj, dft_zmat, R, T, P )
                                mp2_df = pd.concat( [mp2_df, mp2_df_line], ignore_index=True, sort=True ) 
                                   #if isinstance(mp2_series, pd.core.series.Series) :
                                   #   mp2_df = mp2_df.append(mp2_series,ignore_index=True)
                             ## MP2 ENDS

                             ## CCSDT STARTS
                             run_ccs = False
                             if run_ccs and ccs_line.empty and not run_check_ccs and [T,P] == ['90', '90'] and R == '5.0':
                                try:
                                  scan_inp_dict, scan_out_dict, scan_scf, scan_geom = mp2_obj.get_job_results()
                                  ccs_zmat = scan_out_dict['FINAL']['ZMAT']
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
                                   scan_inp_dict, scan_out_dict, scan_scf, scan_geom = mp2_obj.get_job_results()
                                   dft_zmat = scan_out_dict['FINAL']['ZMAT']
                                   hes_series = post_process( 'HES', hes_obj, dft_zmat, R, T, P )
                                   if isinstance(hes_series, pd.core.series.Series) :
                                      hes_df = hes_df.append(hes_series,ignore_index=True)
                                   hes_obj.get_out_dict()
                             ## HES ENDS

                             ## EDA STARTS
                             run_eda = False
                             if run_eda and eda_line.empty and not run_check_eda:
                                scan_inp_dict, scan_out_dict, scan_scf, scan_geom = mp2_obj.get_job_results()
                                dft_zmat = scan_out_dict['FINAL']['ZMAT']
                                ### DFT-D METHODS ARE NOT SUPPORTED. ###
                                eda_series = post_process( 'EDA', eda_obj, dft_zmat, R, T, P )
                                if isinstance(eda_series, pd.core.series.Series) :
                                   eda_df = eda_df.append(eda_series,ignore_index=True)
                             ## EDA ENDS

#                            ## BSSE STARTS
#                            ## BSSE ENDS
                          else:
                            print_tab( 4, err_line['ERROR'].values[0] )

                          RTP_time_end = time.time()
                          if print_time:
                             print_tab( 3, 'Time = {}'.format( RTP_time_end - RTP_time_start ) )

           ## SORT CSV
           for df_obj in [dft_df, mp2_df]:
             
               print(df_obj)
               df_obj.sort_values(by=['Radius'], inplace=True )
               df_obj.reset_index() #drop=True)
 
           ##############
           if fix_errors:
              if not err_df.empty:
                 print_tab( 3, '==========================================' )
                 print_tab( 3, 'FIX ERRORS' )
                 for T in T_SCAN_LIST:
                     print_tab( 3, '==========================================' )
                     for P in P_SCAN_LIST:
                         print_tab( 3, '------------------------------------------' )

                         tmp_err_df = err_df.loc[ err_df['Theta'] == T ].loc[ err_df['Phi'] == P ]
                         tmp_dft_df = dft_df.loc[ dft_df['Theta'] == T ].loc[ dft_df['Phi'] == P ]

                         tmp_err_r_list = list(tmp_err_df['Radius'].astype(float).values)
                         tmp_dft_r_list = list(tmp_dft_df['Radius'].astype(float).values)
                         tmp_err_r_list.sort(reverse = True)
                         tmp_dft_r_list.sort(reverse = True)
                  
                         if len( tmp_err_r_list ) > 0:
                            tmp_err_r_list = list(set( [ tmp_err_r_list[0], tmp_err_r_list[-1] ] ) )
                            for tmp_err_r in tmp_err_r_list:
                                print_tab( 3, 'T = {}, P = {}, R = {} (failed)'.format(T,P,tmp_err_r) )
                                ##
                                print_tab( 4, 'Removing line from scan_err.csv' )
                                tmp_err_idx = err_df.loc[ err_df['Radius'] == float(tmp_err_r) ].loc[ err_df['Theta'] == float(T) ].loc[ 
                                                          err_df['Phi'] == float(P) ].index.values[0]
                                err_df.drop(tmp_err_idx, inplace=True)
                                ##
                                tmp_err_label, tmp_err_obj = get_gms_object( tmp_basis, tmp_funct, T, P, tmp_err_r )
                                run_tmp_err = running_label( tmp_err_obj.inp_name  )
                                if run_tmp_err:
                                   print_tab( 4, 'Running scan' )
                                else:  
                                   ## find closest dft R converged
                                   delta_list = []
                                   for tmp_r in tmp_dft_r_list:
                                       tmp_delta = abs( tmp_r - tmp_err_r )
                                       delta_list.append( (tmp_r, tmp_delta) )
                                   delta_list = sorted( delta_list , key=lambda x: x[1])
                                   for closest_r, closest_d in delta_list:
                                       print( closest_r )
                                       ## extract ZMAT and resubmit
                                       tmp_dft_label, tmp_dft_obj = get_gms_object( tmp_basis, tmp_funct, T, P, closest_r )
                                       tmp_dft_inp_dict, tmp_dft_out_dict, tmp_dft_scf, tmp_dft_geom = tmp_dft_obj.get_job_results()
                                       tmp_dft_zmat = tmp_dft_out_dict['FINAL']['ZMAT']
                                       ## modify value of R in converged ZMAT
                                       tmp_dft_zmat[CAT_NAT]['STR']['val'] = tmp_err_r 
                                       tmp_dft_zmat[CAT_NAT]['BEN']['val'] = T
                                       tmp_dft_zmat[CAT_NAT]['TOR']['val'] = P
                                       print_tab( 4, 'converged.zmat.from.R={}'.format(closest_r) )
                                       #print( skip_closest_test )
                                       #exit()
                                       #tmp_err_obj.fix_error( new_zmat_dict = tmp_dft_zmat, job_queue = 'nodeshiq',
                                       #                       msg = 'converged.zmat.from.R={}'.format(closest_r) )
                                       break
                                   break

           ## PRINT OUT TO CSV
           for df_obj, csv_obj in zip( [ dft_df,             mp2_df,             err_df,             crd_df,             frc_df], 
                                       [ DIMER.scan_dft_csv, DIMER.scan_mp2_csv, DIMER.scan_err_csv, DIMER.scan_crd_csv, DIMER.scan_frc_csv] ):
               df_obj.sort_values(by=['Radius'], inplace=True )
               df_obj.reset_index(drop=True)
               df_obj.to_csv(csv_obj)
               print( csv_obj, len(df_obj) )
 
if __name__ == '__main__':
   main()
