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
import datetime
import glob

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
import itertools

from collections import Counter

import GAMESS
import SLURM

import IONIC_LIQUID as IL

from Functions import print_tab, compose_zmatrices, running_label, center_of_charge, center_of_mass, Coulomb_Energy, Ang2Bohr, now_running
 
#from monomers import change_all_file_names 
from make_scan_list import *
import warnings

if user == 'mdi0316':
   work_dir = '/data/{}/WORK'.format(user)
else:
   work_dir = '/data/scratch-no-backup/{}/WORK'.format(user)

dimers_dir = os.path.join( work_dir, 'DIMERS' )
monomers_dir = os.path.join( work_dir, 'MONOMERS' )
os.makedirs( dimers_dir, exist_ok = True )
temp_dir = '/home/{}/Inputfiles/GAMESS/MONOMERS/AVOGADRO/'.format(user)
mono_json = os.path.join( work_dir, 'monomers_{}.json'.format(user) )

with open(mono_json,'r') as json_file:
  mono_dict = json.load(json_file)

RTP_columns  = ['Radius', 'Theta', 'Phi']
#'INERT.MOM', 'COM', 'COC', 'DIST.COM'

opt_columns = RTP_columns + [ 'exec.err','geom','gms.err','scf' ]
err_columns = RTP_columns + [ 'exec.err','geom','gms.err','scf' ]

crg_col0 = [ (c,c,c) for c in RTP_columns ]
beads_labels = [ 'ANION', 'RING', 'METHYL', 'HC_CHAIN' ]
beads_col    = [ 'pos', 'charge', 'mass' ]
crg_col1 = []
for bead in beads_labels:
  for col in beads_col:
    if col == 'pos':
      for direction in ['x','y','z']:
          crg_col1.append( (bead, col, direction) )
    else:
       crg_col1.append( (bead, col, col) )

crg_columns = pd.MultiIndex.from_tuples(crg_col0 + crg_col1) 

def print_converged_zmat( gms_obj, gms_out_dict ):
    ## converged zmat
    #print( gms_out_dict.keys() )
    if 'FINAL' in gms_out_dict.keys():
       zmat_dict = gms_out_dict['FINAL']['ZMAT']
    else:
       zmat_dict = gms_out_dict['ZMAT']
    zmat_dict = { str(k):v for k, v in zmat_dict.items() }
    zmat_df = dict_to_df( zmat_dict, gms_obj.zmat_file )
    return( zmat_dict, zmat_df )
    
    #fieldnames = zmat_dict.keys()
    #with open( gms_obj.zmat_file, 'w+', newline='') as f:
    #   w = csv.DictWriter(f, fieldnames ) #zmat_dict.keys())
    #   w.writeheader()
    #   w.writerow(zmat_dict)
    #### TO READ
    ##with open( gms_obj.zmat_file, 'r', newline='') as f:
    ##   reader = csv.DictReader(f)
    ##   for row in reader:
    ##       new_dict = dict(row)

def read_zmat( zmat_file ):
    zmat_df = pd.read_csv( zmat_file, index_col = 0 )
    zmat_dict = {} 
    for idx, row in zmat_df.iterrows():
      zmat_dict[str(idx)] = { 'idx.' : row['idx.'], 'elem.' : row['elem.']  }
      for k in [ 'STR', 'BEN', 'TOR' ]:
         if not str(row[k]) == 'nan':
           zmat_dict[str(idx)][k] = ast.literal_eval(row[k])
    return zmat_dict

#def read_zmat( zmat_file ):
#    with open( zmat_file, 'r', newline='') as f:
#       reader = csv.DictReader(f)
#       for row in reader:
#           tmp_dict = dict(row)
#    zmat_dict = {}
#    for i,k in tmp_dict.items():
#       zmat_dict[ i ] = ast.literal_eval(k)
#    return zmat_dict

def fill_df_line( df_line, radius, theta, phi ):
    df_line['Radius'] = float(radius)
    df_line['Theta'] = int(theta)
    df_line['Phi'] = int(phi)
    return df_line

def locate_df_line( df, radius, theta, phi ):
    if df.empty:
       return df
    else:
       return( df.loc[ (df['Radius']==float(radius)) & (df['Theta']==int(theta)) & (df['Phi']==int(phi)) ] )

def dict_to_df( read_dict, csv_file ):
    tmp_keys = [ list(v.keys()) for v in read_dict.values() ] 
    all_keys = list(set(list(itertools.chain(*tmp_keys)))) ## list of unique keys
    df = pd.DataFrame( columns = all_keys )
    for k,v in read_dict.items():
      df = df.append( pd.Series( v ), ignore_index = True )
    df.to_csv( csv_file )
    return df

def print_forces( gms_obj, gms_out_dict ):
    forces_dict = gms_out_dict['FORCES']
    forces_df = dict_to_df( forces_dict, gms_obj.forces_file )
    return forces_df

def print_charges( gms_obj, gms_out_dict ):
    ## atomic coordinates and charges
    print_tab( 4, 'Write charges_results.csv' )
    cart_coords = gms_out_dict['CART.COORDS.']
    charge_dict = gms_out_dict['CHARGE.ANALYSIS']
    charge_df = pd.DataFrame.from_dict( charge_dict, orient = 'index' )
    charge_df.to_csv( gms_obj.charges_file )
    return
 
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

def bead_pos_char_mass( bead_dict, gms_out_dict ):
    charges = gms_out_dict['CHARGE.ANALYSIS'] 
    bead_pos = [[],[],[]]
    bead_mass = 0
    bead_charge = 0
    for k,v in bead_dict.items():
      bead_charge += [ q for (p,q) in charges.items() if int(k)==int(p) and v['elem.'] == q['elem.'] ][0]['charge'] 
      v_elem = v['elem.']
      bead_mass += element(v_elem).mass
      bead_pos[0].append(v['x'])
      bead_pos[1].append(v['y'])
      bead_pos[2].append(v['z'])
    bead_pos[0] = np.array(bead_pos[0]).mean()
    bead_pos[1] = np.array(bead_pos[1]).mean()
    bead_pos[2] = np.array(bead_pos[2]).mean()

    return bead_pos, bead_mass, bead_charge
    
def print_beads_coords( gms_obj, gms_out_dict ):
    ## center of mass
    cart_coords = gms_out_dict['CART.COORDS.']

    C_cart_coords = { k:v for (k,v) in cart_coords.items() if v['elem.'] == 'C' } 
    H_cart_coords = { k:v for (k,v) in cart_coords.items() if v['elem.'] == 'H' } 

    ## define sphere corresponding to ring
    N_coords = [ [v['x'], v['y'], v['z']] for v in cart_coords.values() if v['elem.'] == 'N' ]
    ring_center = np.mean(N_coords, axis=0) # center of the ring
    C_dist_from_ring_center = []
    for v in C_cart_coords.values():
       tmp_coords = np.array( [v['x'], v['y'], v['z']] )
       tmp_dist = LA.norm( tmp_coords - ring_center ) 
       C_dist_from_ring_center.append( [v['idx.'],tmp_dist] )
    C_dist_from_ring_center.sort(key=lambda x : x[1] )

    C_ring_idxs = [ int(cc[0]) for cc in C_dist_from_ring_center[:3] ] 
    N_ring_idxs = [ int(v['idx.']) for v in cart_coords.values() if v['elem.'] == 'N' ]
    C_outer_idxs = [ int(cc[0]) for cc in C_dist_from_ring_center[3:] ] 

    ring_dict = { k:v for (k,v) in cart_coords.items() if int(v['idx.']) in C_ring_idxs + N_ring_idxs }

    ## write position of two more extern Carbons 
    CC_dist_dict = { k:v for (k,v) in gms_out_dict['INTERNUCL.DISTANCES'].items() if v['at.1']['elem.1'] == 'C' and v['at.2']['elem.2'] == 'C' }
    CC_max_dist = max(CC_dist_dict, key=lambda k: CC_dist_dict[k]['dist.']) 
    C_1, C_2, dist = CC_dist_dict[CC_max_dist].values()
    C_1_idx = int(C_1['idx.1'])
    C_2_idx = int(C_2['idx.2'])
    C_1_dict = [ v for v in C_cart_coords.values() if int(v['idx.']) == C_1_idx ][0]
    C_2_dict = [ v for v in C_cart_coords.values() if int(v['idx.']) == C_2_idx ][0]
    C_1_coords = np.array([ C_1_dict['x'], C_1_dict['y'], C_1_dict['z'] ])
    C_2_coords = np.array([ C_2_dict['x'], C_2_dict['y'], C_2_dict['z'] ])
    bead_1_dict = { k:v for (k,v) in C_cart_coords.items() if int(v['idx.']) == C_1_idx }
    bead_2_dict = { k:v for (k,v) in C_cart_coords.items() if int(v['idx.']) == C_2_idx }
    
    ## find ethyl/methyl groups 
    for C_k, C_v in C_cart_coords.items():
        C_idx = int(C_v['idx.'])
        if C_idx not in C_ring_idxs + [C_1_idx, C_2_idx]:
           C_dist = []
           for p,q in CC_dist_dict.items():
              idx1 = int(q['at.1']['idx.1'])
              idx2 = int(q['at.2']['idx.2'])
              if idx1 == C_idx and idx2 in [C_1_idx, C_2_idx]:
                 C_dist.append([ idx2, q['dist.'] ])
           C_dist.sort( key= lambda x:x[1] )
           closest_C_idx = C_dist[0][0]

           if closest_C_idx == C_1_idx:
              bead_1_dict[C_k] = C_v
           elif closest_C_idx == C_2_idx:
              bead_2_dict[C_k] = C_v

    ## attach H atoms
    for H_k, H_v in H_cart_coords.items():
        H_coords = np.array( [ H_v['x'], H_v['y'], H_v['z'] ]) 
        H_distances = [ LA.norm( H_coords-ring_center), LA.norm(H_coords-C_1_coords), LA.norm(H_coords-C_2_coords) ]
        if   H_distances.index(min(H_distances)) == 0: ## H belongs to ring
           ring_dict[H_k] = H_v
        elif H_distances.index(min(H_distances)) == 1: ## H belongs to bead_1
           bead_1_dict[H_k] = H_v
        elif H_distances.index(min(H_distances)) == 2: ## H belongs to bead_2
           bead_2_dict[H_k] = H_v

    if len(bead_1_dict) == 4 and len(bead_2_dict) >= len(bead_1_dict):
       methyl_dict = dict(bead_1_dict)
       hc_chain_dict = dict(bead_2_dict)
    elif len(bead_2_dict) == 4 and len(bead_1_dict) >= len(bead_2_dict):
       methyl_dict = dict(bead_2_dict)
       hc_chain_dict = dict(bead_1_dict)
    else:
       #print( bead_1_dict )
       #print( bead_2_dict )
       print( 'cannot_recognize_methyl_and_carbon_chain' )
       beads_df = pd.DataFrame( columns = crg_columns )
       beads_df.to_csv( gms_obj.beads_file ) 
       return beads_df

    ## Anion bead 
    anion_dict = { k:v for (k,v) in cart_coords.items() if int(k) >= CAT_NAT }

    ## write to out
    anion_pos , anion_mass , anion_charge  = bead_pos_char_mass( anion_dict,  gms_out_dict )
    ring_pos  , ring_mass  , ring_charge   = bead_pos_char_mass( ring_dict,   gms_out_dict )
    methyl_pos, methyl_mass, methyl_charge = bead_pos_char_mass( methyl_dict, gms_out_dict )
    hc_chain_pos, hc_chain_mass, hc_chain_charge = bead_pos_char_mass( hc_chain_dict, gms_out_dict )

    beads_dict = {
      ('ANION','pos','x') : anion_pos[0], ('ANION','pos','y') : anion_pos[1], ('ANION', 'pos', 'z') : anion_pos[2],    
      ('ANION','mass','mass') : anion_mass, ('ANION','charge','charge') : anion_charge,

      ('RING','pos','x') : ring_pos[0], ('RING','pos','y') : ring_pos[1], ('RING', 'pos', 'z') : ring_pos[2],
      ('RING','mass','mass') : ring_mass, ('RING','charge','charge') : ring_charge, 

      ('METHYL','pos','x') : methyl_pos[0], ('METHYL','pos','y') : methyl_pos[1], ('METHYL','pos','z') : methyl_pos[2],   
      ('METHYL', 'mass', 'mass') : methyl_mass, ('METHYL', 'charge', 'charge') : methyl_charge,   

      ('HC_CHAIN', 'pos', 'x') : hc_chain_pos[0], ('HC_CHAIN', 'pos', 'y') : hc_chain_pos[1], 
      ('HC_CHAIN', 'pos', 'z') : hc_chain_pos[2], ('HC_CHAIN', 'mass', 'mass') : hc_chain_mass, 
      ('HC_CHAIN', 'charge', 'charge') : hc_chain_charge
                  }
    beads_df = pd.DataFrame( [beads_dict], columns = crg_columns )
    beads_df.to_csv( gms_obj.beads_file, index=False ) 
     
    return( beads_df )

def calculate_dcom( gms_out_dict, cat_nat ):
    cart_coords  = gms_out_dict['CART.COORDS.']   #in Bohr
    charges = gms_out_dict['CHARGE.ANALYSIS'] 
    com = center_of_mass(   cart_coords, length_factor = 1/Ang2Bohr ) 
    coc = center_of_charge( cart_coords, charges, length_factor = 1/Ang2Bohr ) 
    cat_coords = { k:v for (k,v) in cart_coords.items() if int(k) <  cat_nat } 
    ani_coords = { k:v for (k,v) in cart_coords.items() if int(k) >= cat_nat }
    cat_com = center_of_mass( cat_coords, length_factor = 1/Ang2Bohr ) 
    ani_com = center_of_mass( ani_coords, length_factor = 1/Ang2Bohr ) 
    dcom = LA.norm( np.array(cat_com)-np.array(ani_com ))
    return( com, coc, cat_com, ani_com, dcom )

def print_ionic_charges( gms_out_dict, cat_nat ):
    mull_charge_cat = sum([float(v['Mull.charge']) for (k,v) in gms_out_dict['CHARGE.ANALYSIS'].items() if int(k)< cat_nat])
    mull_charge_ani = sum([float(v['Mull.charge']) for (k,v) in gms_out_dict['CHARGE.ANALYSIS'].items() if int(k)>=cat_nat])
    lowd_charge_cat = sum([float(v['Lowd.charge']) for (k,v) in gms_out_dict['CHARGE.ANALYSIS'].items() if int(k)< cat_nat])
    lowd_charge_ani = sum([float(v['Lowd.charge']) for (k,v) in gms_out_dict['CHARGE.ANALYSIS'].items() if int(k)>=cat_nat])
    return( mull_charge_cat, mull_charge_ani, lowd_charge_cat, lowd_charge_ani )

def print_dft_results( gms_obj, gms_out_dict, dimer, distance, zero_dft_ener ):
    if os.path.exists( gms_obj.dft_file ):
       print_tab( 4, 'Read dft_results.csv' )
       dft_df = pd.read_csv( gms_obj.dft_file, index_col = 0  )
    else:
       dft_dict = { k:v for (k,v) in gms_out_dict.items() if k in [ 'SCF', 'TOT.EN.', 'GEOM.', 'BANDGAP', 'BASIS.DIM.'] }
       if dimer:
          mull_charge_cat, mull_charge_ani = print_ionic_charges( gms_out_dict, gms_obj.nat_cat )[:2]
          com, coc, cat_com, ani_com, dcom = calculate_dcom( gms_out_dict, gms_obj.nat_cat )
          dft_dict['MULL.CHARG.CAT.'] = mull_charge_cat
          dft_dict['MULL.CHARG.ANI.'] = mull_charge_ani
          dft_dict['DIST.COM']     = dcom
          dft_dict['DFT.INT.EN.']  = dft_dict['TOT.EN.'] - zero_dft_ener
          dft_dict['COUL.EN.R.']   = Coulomb_Energy( distance, mull_charge_cat, mull_charge_ani )
          dft_dict['COUL.EN.COM.'] = Coulomb_Energy( dcom,  mull_charge_cat, mull_charge_ani )
          dft_dict['DISP.EN.R.']   = dft_dict['DFT.INT.EN.'] - dft_dict['COUL.EN.R.']
          dft_dict['DISP.EN.COM.'] = dft_dict['DFT.INT.EN.'] - dft_dict['COUL.EN.COM.']
       print_tab( 4, 'Write dft_results.csv' )
       dft_df = pd.DataFrame( [dft_dict] )
       dft_df.to_csv( gms_obj.dft_file )
    return( dft_df )

def print_mp2_results( gms_obj, gms_out_dict, dimer, distance, zero_mp2_ener ):
    if os.path.exists( gms_obj.mp2_file ):
       print_tab( 4, 'Read mp2_results.csv' )
       mp2_df = pd.read_csv( gms_obj.mp2_file, index_col = 0  )
    else:
       mp2_dict = gms_out_dict['MP2']
       if dimer:
          mull_charge_cat, mull_charge_ani = print_ionic_charges( gms_out_dict, gms_obj.nat_cat )[:2]
          com, coc, cat_com, ani_com, dcom = calculate_dcom( gms_out_dict, gms_obj.nat_cat )
          mp2_dict['MULL.CHARG.CAT.'] = mull_charge_cat
          mp2_dict['MULL.CHARG.ANI.'] = mull_charge_ani
          mp2_dict['DIST.COM']     = dcom
          mp2_dict['MP2.INT.EN.']  = mp2_dict['EN.MP2'] - zero_mp2_ener 
          mp2_dict['COUL.EN.R.']   = Coulomb_Energy( distance, mull_charge_cat, mull_charge_ani )
          mp2_dict['COUL.EN.COM.'] = Coulomb_Energy( dcom,  mull_charge_cat, mull_charge_ani )
          mp2_dict['DISP.EN.R.']   = mp2_dict['MP2.INT.EN.'] - mp2_dict['COUL.EN.R.']
          mp2_dict['DISP.EN.COM.'] = mp2_dict['MP2.INT.EN.'] - mp2_dict['COUL.EN.COM.']
       print_tab( 4, 'Write mp2_results.csv' )
       mp2_df = pd.DataFrame( [mp2_dict] )
       mp2_df.to_csv( gms_obj.mp2_file )
    return( mp2_df )

def print_ccsdt_results( gms_obj, gms_out_dict, dimer ):
    ccsdt_dict = {} 
    ccsdt_dict['REF.EN']  = gms_out_dict['CCSDT']['REF.EN.']
    ccsdt_dict['MBPT(2)'] = gms_out_dict['CCSDT']['MBPT(2)']['EN.']
    ccsdt_dict['CCSD']    = gms_out_dict['CCSDT']['CCSD']['EN.']
    ccsdt_dict['CCSD(T)'] = gms_out_dict['CCSDT']['CCSD(T)']['EN.']
    ccsdt_dict['CCSD[T]'] = gms_out_dict['CCSDT']['CCSD[T]']['EN.']
    ccsdt_df = pd.DataFrame( [ccsdt_dict] )
    ccsdt_df.to_csv( gms_obj.ccsdt_file )
    return ccsdt_df

def get_gms_object( dmr_label, basis, funct, T, P, R, post_scf = 'DFTTYP', run_type = 'OPTIMIZE', run_dir = False, opt_from = False ):
    
    tmp_ifreeze = None
    dmr_obj, cat_nat, ani_nat, zero_dft_en, zero_mp2_en = read_dimer( dmr_label, basis, funct )
    if run_type == 'OPTIMIZE' and post_scf == 'DFTTYP':
       cat_zmat_dim = 3*cat_nat-6
       tmp_ifreeze = '{},{},{}'.format( cat_zmat_dim + 1,  cat_zmat_dim + 2, cat_zmat_dim + 3 ) 
       run_dir = os.path.join( dmr_obj.opt_dir, 'DFT', basis, funct, 'T_%s'%T, 'P_%s'%P, 'R_%s'%R )
    elif run_type == 'ENERGY' and post_scf == 'DFTTYP':
       run_dir = os.path.join( dmr_obj.ene_dir, 'OPT_from_{}'.format(opt_from), 'DFT', basis, funct, 'T_%s'%T, 'P_%s'%P, 'R_%s'%R )
    elif run_type == 'ENERGY' and post_scf == 'MP2':
       run_dir = os.path.join( dmr_obj.ene_dir, 'OPT_from_{}'.format(opt_from), 'MP2', basis, 'T_%s'%T, 'P_%s'%P, 'R_%s'%R )
    elif run_type == 'EDA' and post_scf == 'DFTTYP':
       run_dir = os.path.join( dmr_obj.eda_dir, 'OPT_from_{}'.format(opt_from), 'DFT', basis, funct, 'T_%s'%T, 'P_%s'%P, 'R_%s'%R )
    elif run_type == 'EDA' and post_scf == 'MP2':
       run_dir = os.path.join( dmr_obj.eda_dir, 'OPT_from_{}'.format(opt_from), 'MP2', basis, 'T_%s'%T, 'P_%s'%P, 'R_%s'%R )

    if post_scf == 'DFTTYP':
       gms_inp = 'gms_SCAN_{}_{}_DFT_{}_{}_T_{}_P_{}_R_{}.inp'.format( DIM_LABEL.lower(), run_type[:3], basis, funct, T, P, R )
    elif post_scf == 'MP2':
       gms_inp = 'gms_SCAN_{}_{}_MP2_{}_T_{}_P_{}_R_{}.inp'.format(    DIM_LABEL.lower(), run_type[:3], basis, T, P, R )

    gms_obj = GAMESS.GAMESS( inp_name = gms_inp, run_dir = run_dir, 
                             natoms = cat_nat + ani_nat, nat_cat = cat_nat, nat_ani = ani_nat, 
                             icharge = 0, run_type = run_type, post_scf = post_scf,
                             basis = basis, functional = funct, ifreeze = tmp_ifreeze )

    return( gms_inp, gms_obj )

def read_dimer( dim_label, basis, funct, get_zmat = False ):

    dimer = IL.DIMER( dim_label, basis, funct )
    cat_label, ani_label = dim_label.split('_')
    cat_nat = dimer.cat_dict['nat']
    ani_nat = dimer.ani_dict['nat']

    cat_csv_file = os.path.join( monomers_dir, 'CSV', '{}.csv'.format( cat_label ) )
    cat_df = pd.read_csv( cat_csv_file, index_col = 0 )
    cat_line = cat_df.loc[ cat_df['BASIS'] == basis ].loc[ cat_df['FUNCT'] == funct ]

    cat_dft_en, cat_mp2_en = 0., 0.
    if cat_line.empty:
       print( 'WARNING: empty cat_line' )
    else:
       cat_dft_en = cat_line['DFT.TOT.EN.'].values[0] 
       cat_mp2_en = cat_line['MP2.TOT.EN.'].values[0]

    ani_csv_file = os.path.join( monomers_dir, 'CSV', '{}.csv'.format( ani_label ) )
    ani_df = pd.read_csv( ani_csv_file, index_col = 0 )
    ani_line = ani_df.loc[ ani_df['BASIS'] == basis ].loc[ ani_df['FUNCT'] == funct ]

    if ani_line.empty:
       print( 'WARNING: empty ani_line' )
    else:
       ani_dft_en = ani_line['DFT.TOT.EN.'].values[0]
       ani_mp2_en = ani_line['MP2.TOT.EN.'].values[0]

    zero_dft_en = cat_dft_en + ani_dft_en 
    zero_mp2_en = cat_mp2_en + ani_mp2_en 

    if get_zmat:
       cat_IL  = IL.MONOMER( cat_label, basis, funct )
       cat_opt = GAMESS.GAMESS( inp_label = '{}_{}_{}'.format( cat_label.lower(), basis, funct ), 
                                root_dir = cat_IL.root_dir, natoms = cat_nat, 
                                run_type = 'OPTIMIZE', post_scf = 'DFTTYP', basis = basis, functional = funct )
       cat_zmat = read_zmat( cat_opt.zmat_file )

       ani_IL  = IL.MONOMER( ani_label, basis, funct )
       ani_opt = GAMESS.GAMESS( inp_label = '{}_{}_{}'.format( ani_label.lower(), basis, funct ), 
                                root_dir = ani_IL.root_dir, natoms = ani_nat, 
                                run_type = 'OPTIMIZE', post_scf = 'DFTTYP', basis = basis, functional = funct )
       ani_zmat = read_zmat( ani_opt.zmat_file )
       return( dimer, cat_nat, ani_nat, zero_dft_en, zero_mp2_en, cat_zmat, ani_zmat )

    if np.isnan( zero_dft_en ) or np.isnan( zero_mp2_en ):
       warnings.warn("Monomer Energy is Nan")
       zero_dft_en, zero_mp2_en = 2*[False]

    return( dimer, cat_nat, ani_nat, zero_dft_en, zero_mp2_en )

def read_object( read_obj, dimer=True, read_dict=False, read_msg=False, read_template=False, print_msg=True ):
    print_tab( 3, '---- Read gms.obj starts ({}/{}) ----'.format(read_obj.run_type, read_obj.post_scf) )
    read_exec, read_exec_err, read_gms_err, read_scf, read_geom, read_time = 6*[False]
    t0 = time.time()
    ## check run
    read_run = running_label( read_obj.inp_name )
    if read_run:
       print_tab( 4, '{}: Running'.format(read_obj.run_type) )
       read_exec = 'Running'
    else:
       ## missing dir/inp_file/out_file
       if os.path.isdir( read_obj.run_dir ):
          if os.path.exists( read_obj.inp_file ):
             ## missing out
             if read_obj.out_file == 'MISSING':
                print_tab( 4, 'WARNING: missing output file' )
                read_obj.resubmit()
                read_exec = 'Running'
             ## check all csv files
             elif os.path.exists( read_obj.status_file ):
                print_tab( 4, 'read {}'.format(read_obj.status_file) )
                status_df     = pd.read_csv( read_obj.status_file, index_col=0 )
                read_exec     = status_df['EXEC.'    ].values[0]
                read_exec_err = status_df['EXEC.ERR.'].values[0]
                read_gms_err  = status_df['GMS.ERR.' ].values[0]
                read_scf      = status_df['SCF'      ].values[0]
                read_geom     = status_df['GEOM.'    ].values[0]
                read_time     = status_df['TIME'     ].values[0]
                
                if read_exec in ['Running', 'NEW']:
                   read_exec = None
          else:
             ## missing inp
             print_tab( 4, 'looking for {} in\n {}'.format( read_obj.inp_file, read_obj.run_dir ) )
             print_tab( 4, 'WARNING: missing input file' )
             read_exec = 'NEW'
       else:
          ## missing dir
          print_tab( 4, 'WARNING: missing folder' )
          read_exec = 'NEW'

    print( 'read_exec', 'read_exec_err', 'read_gms_err' )
    print(  read_exec,   read_exec_err,   read_gms_err  )
    ## read GAMESS and fix output 
    if   read_exec in [ 'Running', 'NEW' ]:
         pass #will write csv
    # OPT case not ok: will not write csv 
    elif read_obj.run_type == 'OPTIMIZE' and [ read_exec, read_exec_err, read_gms_err, read_scf, read_geom ] ==\
                                             [ 'TERMINATED.NORMALLY', False, False, 'UNCONVERGED', 'NOT.LOCATED']:
         pass #return read_exec, read_exec_err, read_gms_err, read_scf, read_geom, read_time
    # OPT case ok: will not write csv 
    elif read_obj.run_type == 'OPTIMIZE' and [ read_exec, read_exec_err, read_gms_err, read_scf, read_geom ] ==\
                                             [ 'TERMINATED.NORMALLY', False, False, 'CONVERGED', 'LOCATED']:
         pass #return read_exec, read_exec_err, read_gms_err, read_scf, read_geom, read_time
    # ENE case ok: will not write csv 
    elif read_obj.run_type != 'OPTIMIZE' and [ read_exec, read_exec_err, read_gms_err, read_scf, read_geom ] ==\
                                             [ 'TERMINATED.NORMALLY', False, False, 'CONVERGED', False]:
         pass #return read_exec, read_exec_err, read_gms_err, read_scf, read_geom, read_time
    elif read_exec in ['TERMINATED.ABNORMALLY', 'TERMINATED.ERROR']:
         read_gms_err = read_obj.read_error()
         if   read_gms_err == 'max. ang. momentum exceeded':
              read_obj.fix_ang_momentum_error()
              read_exec = 'Running'
         elif read_gms_err in [ 'insufficient.distributed.memory', 'insufficient.replicated.memory', 
                                'memory.request.exceeds.available.memory','not.enough.replicated.memory' ]:
              read_obj.fix_memory_error( read_exec_err )
              read_exec = 'Running'
         elif any( [ err in [read_gms_err, read_exec_err] for err in [ 
                               'Gradient.out.of.range', 'Stationary.Point.Location.failed', 
                               'Serious.Failure', 'Error.Einval', 'Error.Numerical.Gradient', 
                               'Error.Semop', 'Error.Shmget' ] ] ):
              pass
         elif 'dawrit' in str(read_gms_err) or 'dawrit' in str(read_exec_err):
              pass
         else: 
              print( 'dawrit' in [read_gms_err, read_exec_err] )
              print( 'read_exec', 'read_exec_err', 'read_gms_err' )
              print(  read_exec,   read_exec_err,   read_gms_err  )
              print( fix_error )
         return read_exec, read_exec_err, read_gms_err, read_scf, read_geom, read_time
    else:
         t3 = time.time()
         print_tab( 4, 'Read output file' )
         read_exec, read_exec_err = read_obj.get_job_exec()
         print( 'read_exec', 'read_exec_err', 'read_gms_err' )
         print(  read_exec,   read_exec_err,   read_gms_err  )
         read_time = read_obj.read_wall_clock()
         if [ read_exec, read_exec_err ] == ['MISSING.OUTPUT.FILE', False]:
              read_obj.resubmit()
              read_exec = 'RESUBMITTED'
         elif read_exec in ['TERMINATED.ABNORMALLY', 'TERMINATED.ERROR']:
              read_gms_err = read_obj.read_error()
         elif [ read_exec, read_exec_err ] == ['TERMINATED.NORMALLY', False]:
              read_inp_dict, read_out_dict, read_scf, read_geom = read_obj.get_job_results()
              if read_scf == 'CONVERGED':
                 if read_obj.run_type == 'OPTIMIZE':
                    print_converged_zmat( read_obj, read_out_dict )
              else:
                 print_tab(4, 'WARNING: read_scf: {}'.format(read_scf) )
         elif [ read_exec, read_exec_err ] == [False, False]:
              print_tab(4, [read_exec, read_exec_err] )
         else:
              print( read_exec, read_exec_err )
              raise RuntimeError("unknown read_exec {} in {}".format(read_exec, read_obj.run_dir))
         print_tab( 4, 'Read output file: {} sec.'.format(time.time()-t3) )

    ## MAKE NEW
    if read_exec == 'NEW':
       #go = input('proceed? (Y/N)')
       #go = 'Y' #input('proceed? (Y/N)')
       #if go == 'Y':
         running_now = now_running()
         if  running_now > 1000:
             print( 'too many jobs running' )
             return 6*[False]
         else:
             ## default, create zmat_dict from relaxed monomers
             if dimer and read_dict == False and read_template == False:
                cat_zmat, ani_zmat = read_dimer( DIM_LABEL, read_obj.basis, read_obj.functional, get_zmat = True )[5:7]
                read_dict = compose_zmatrices( cat_zmat, ani_zmat, radius=R, theta=T, phi=P )

             ## run monomers from template
             if bool(read_dict) == False and bool(read_template) == True:
                print_tab(4, 'printing new from _template_')
                read_obj.run_new( zmat_dat = read_template, msg = read_msg ) 
 
             ## run dimers/monomers from zmatrix
             elif bool(read_dict) == True  and bool(read_template) == False:
                print_tab(4, 'printing new from _dict_')
                read_obj.run_new( zmat_dict = read_dict, msg = read_msg ) 

             else:
               print( read_dict, read_template )
               raise NameError( 'read_dict and read_template both True' )
    
    t1 = time.time()
    if print_msg:
       print_tab( 4, 'Write status.csv' )
       print_tab( 4, 'exec. : {}, exec.err. : {}, gms.err. : {}'.format( read_exec, read_exec_err, read_gms_err ) )
       if read_obj.run_type == 'OPTIMIZE':
          print_tab( 4, 'scf : {}, geom : {}, time: {}'.format( read_scf, read_geom, read_geom ) )
       print_tab( 3, '---- time: {} sec.'.format(t1-t0) )

    status_df = pd.DataFrame( [ { 'TYPE' : read_obj.run_type, 'EXEC.' : read_exec, 'EXEC.ERR.' : read_exec_err, 
                                  'GMS.ERR.': read_gms_err, 'SCF' : read_scf,  'GEOM.' : read_geom, 'TIME' : read_time } ] )
    status_df.to_csv( read_obj.status_file )

    return read_exec, read_exec_err, read_gms_err, read_scf, read_geom, read_time

def single_point_calculations( dim_label, rlx_basis, rlx_funct, sp_basis, sp_funct, post_scf = 'DFTTYP' ):

    rlx_dimer, cat_nat, ani_nat, zero_dft_ener, zero_mp2_ener = read_dimer( DIM_LABEL, rlx_basis, rlx_funct )
    rlx_opt_df = pd.read_csv( rlx_dimer.scan_opt_csv, index_col=0 )
    
    sp_dimer = IL.DIMER( dim_label, sp_basis, sp_funct )
    r_sp_list, t_sp_list, p_sp_list = get_scan_list( DIM_LABEL, sp_basis, sp_funct )

   # t_sp_list = ['90'] 
   # p_sp_list = ['90'] 

    if post_scf == 'DFTTYP':
       ene_csv = os.path.join( sp_dimer.csv_dir, sp_basis, sp_funct, 'scan_dft_ene.csv' )
       eda_csv = os.path.join( sp_dimer.csv_dir, sp_basis, sp_funct, 'scan_dft_eda.csv' )
    elif post_scf == 'MP2':
       ene_csv = os.path.join( sp_dimer.csv_dir, sp_basis, 'scan_mp2_ene.csv' )
       eda_csv = os.path.join( sp_dimer.csv_dir, sp_basis, 'scan_mp2_eda.csv' )
      
    if os.path.exists( ene_csv ):
       ene_df = pd.read_csv( ene_csv, index_col=0 )
    else:
       os.makedirs( os.path.join( sp_dimer.csv_dir, sp_basis, sp_funct ), exist_ok = True )
       ene_df = pd.DataFrame() 

    if os.path.exists( eda_csv ):
       eda_df = pd.read_csv( eda_csv, index_col=0 )
    else:
       eda_df = pd.DataFrame() 
 
    print_tab(1, '=== SINGLE POINT {} ENERGY {}, {} (ON TOP OF: {}, {}) ==='.format( post_scf[:3], sp_basis, sp_funct, rlx_basis, rlx_funct) )
    print_tab(2, '{} has {} lines'.format(rlx_dimer.scan_opt_csv, len(rlx_opt_df)))
    print_tab(2, '{} has {} lines'.format(ene_csv, len(ene_df)))
    print_tab(2, '{} has {} lines'.format(eda_csv, len(eda_df)))
    
    for rlx_idx, rlx_row in rlx_opt_df.iterrows():
        rlx_R, rlx_T, rlx_P = str(rlx_row['Radius']), str(rlx_row['Theta']), str(rlx_row['Phi'])
        if rlx_R in r_sp_list and  rlx_T in t_sp_list and rlx_P in p_sp_list:
           ene_line = locate_df_line( ene_df, rlx_R, rlx_T, rlx_P )
           eda_line = locate_df_line( eda_df, rlx_R, rlx_T, rlx_P )
           print_tab( 3, '{}, {}, T = {}, P = {}, R = {}'.format( sp_basis, sp_funct, rlx_T, rlx_P, rlx_R) )
           if ene_line.empty or ( eda_line.empty and [ rlx_T , rlx_P ] == [ '90', '90'] ) :
              rlx_label, rlx_obj = get_gms_object( DIM_LABEL, rlx_basis, rlx_funct, rlx_T, rlx_P, rlx_R )
              rlx_out_zmat = read_zmat(rlx_obj.zmat_file) 
              print_tab( 3, '---- {}/ENE'.format(post_scf) )
              ene_obj = get_gms_object( DIM_LABEL, sp_basis, sp_funct, rlx_T, rlx_P, rlx_R, 
                                        post_scf = post_scf, run_type = 'ENERGY', 
                                        opt_from = 'DFT_{}_{}'.format(rlx_basis, rlx_funct))[1]
              ene_exec, ene_exec_err, ene_gms_err, ene_scf, ene_geom, ene_time = read_object( ene_obj, 
                                        read_dict = rlx_out_zmat, read_msg = 'OPT_from_{}_{}'.format(rlx_basis, rlx_funct) )
              if [ ene_exec, ene_exec_err, ene_gms_err, ene_scf ] == [ 'TERMINATED.NORMALLY', False, False, 'CONVERGED' ]:
                 if ene_line.empty:
                    ene_out_dict = ene_obj.get_job_results()[1]
                    if post_scf == 'DFTTYP':
                         ene_line = print_dft_results( ene_obj, ene_out_dict, dimer=True, distance=float(rlx_R), zero_dft_ener=zero_dft_ener ) 
                    elif post_scf == 'MP2':
                         ene_line = print_mp2_results( ene_obj, ene_out_dict, dimer=True, distance=float(rlx_R), zero_mp2_ener=zero_mp2_ener ) 
                    ene_line = fill_df_line( ene_line, rlx_R, rlx_T, rlx_P )
                    ene_df   = ene_df.append( ene_line, ignore_index=True, sort=True )
                    print_charges( ene_obj, ene_out_dict )

                 ## EDA starts
                 if rlx_T == '90' and rlx_P == '90' :
                 #if rlx_T == '90' and rlx_P in [ '0', '90', '180', '270' ]: 
                    print_tab( 3, '---- {}/EDA'.format(post_scf) )
                    eda_obj = get_gms_object( DIM_LABEL, sp_basis, sp_funct, rlx_T, rlx_P, rlx_R, 
                                              post_scf = post_scf, run_type = 'EDA', 
                                              opt_from = 'DFT_{}_{}'.format(rlx_basis, rlx_funct))[1]
                    eda_exec, eda_exec_err, eda_gms_err, eda_scf, eda_geom, eda_time = read_object( eda_obj, 
                                              read_dict = rlx_out_zmat, read_msg = 'OPT_from_{}_{}'.format(rlx_basis, rlx_funct) )
                    if [ eda_exec, eda_exec_err, eda_gms_err, eda_scf ] == [ 'TERMINATED.NORMALLY', False, False, 'CONVERGED' ]:
                       eda_line = pd.DataFrame( [eda_obj.get_job_results()[1]['EDA']] )
                       eda_line = fill_df_line( eda_line, rlx_R, rlx_T, rlx_P )
                       eda_df = eda_df.append( eda_line, ignore_index=True, sort=True )
                 ## EDA ends
           else :
              print_tab( 4, 'complete' )

    for df_obj, csv_file in zip( [ene_df, eda_df], [ene_csv, eda_csv] ):
        if not df_obj.empty:
           df_obj.sort_values(by=['Theta','Phi','Radius'], inplace=True )
           df_obj.drop_duplicates(inplace=True )
           df_obj.reset_index(drop=True, inplace=True)
           df_obj.to_csv(csv_file)
           print_tab( 3, 'printing to {}'.format(csv_file) )
        else:
           print_tab( 3, 'empty df, skip printing {}'.format(csv_file) )
    print_tab(2, '--- SINGLE POINT ENERGY --- ends' )
    return 

def read_df_line( df_line ):
    return df_line['Radius'], df_line['exec.err'], df_line['geom'], df_line['gms.err'], df_line['scf'], df_line['exec.']

def R_distances( radius_list, radius_val ):
    distance_list = [ (tmp_r, abs(float(tmp_r)-radius_val)) for tmp_r in list(set(radius_list) -set([str(radius_val)])) ]
    distance_list.sort(key=lambda x : x[1])
    return distance_list
    
def make_err_dir( fail_obj, fail_R, clst_r ):
    print_tab( 4, 'make_err_dir starts' )
    t5 = time.time()
    fail_exec, fail_exec_err, fail_gms_err, fail_scf, fail_geom, fail_time = read_object( fail_obj, print_msg=False ) 
    now = datetime.datetime.now()
    modify_inp_file = os.path.join( fail_obj.run_dir, 'modify_inp.dat' )

    if fail_exec == 'NEW':
       print( cannot_be_new )
    else:
      error_dir = os.path.join( fail_obj.run_dir, 'FAILED', 'zmat_from_{}'.format(clst_r) )
      if os.path.exists( error_dir ):
         print_tab( 5, 'warning error dir exists in {}'.format( fail_obj.run_dir ) )
         exec_exec = 'go_to_next_radius'

      else:
         if os.path.exists( modify_inp_file ):
            print_tab( 5, 'modify_inp.dat file exists in {}'.format( fail_obj.run_dir ) )
            open_mode = 'a'
            with open( modify_inp_file, 'r' ) as mod_inp:
                 modify_inp_lines = mod_inp.readlines()
            previous_zmat =  [ll for ll in modify_inp_lines if ll.startswith( 'zmat from') ][-1].strip().replace(' ','_')
         else:
            open_mode = 'w+'
            previous_zmat =  'zmat_from_monomers'

         error_dir = os.path.join( fail_obj.run_dir, 'FAILED', previous_zmat,
                                   now.strftime("%c").replace(' ', '_').replace('__', '_'), 
                                   str(fail_exec), str(fail_exec_err), 'gms_err_{}'.format(fail_gms_err) )
         print_tab( 3 , 'creating error dir = {}'.format( error_dir ) )
         os.makedirs( error_dir )
         print_tab( 3 , 'updating modify_inp.dat' )
         with open( modify_inp_file, open_mode ) as mod_inp:
              mod_inp.write( '{}\n'.format(now.strftime("%c")) )
              mod_inp.write( 'exec = {}, exec_err = {}, gms_err = {}, fail_scf = {}, fail_geom = {}\n'.format(
                              fail_exec, fail_exec_err,  fail_gms_err, fail_scf, fail_geom ))
              mod_inp.write( 'zmat from {}\n'.format(clst_r) )
         print_tab( 3 , 'moving files to error dir' )
         for ff in os.listdir( fail_obj.run_dir ):
             abs_ff = os.path.join( fail_obj.run_dir, ff )
             if os.path.isfile( abs_ff ):
                if ff.endswith('inp') or ff.startswith('log') or ff.startswith('status'):
                   shutil.move( abs_ff, error_dir )
                elif ff == 'modify_inp.dat':
                   pass
                else:
                   os.remove( abs_ff )
         print_tab( 3 , 'writing new file from clst_r = {}'.format(clst_r) )

         clst_obj = get_gms_object( DIM_LABEL, fail_obj.basis, fail_obj.functional, T, P, clst_r )[1]
         clst_zmat = read_zmat( clst_obj.zmat_file ) 
         clst_zmat[str( fail_obj.nat_cat )]['STR']['val'] = fail_R

         new_obj = get_gms_object(  DIM_LABEL, fail_obj.basis, fail_obj.functional, T, P, fail_R )[1]
         read_object( new_obj, read_dict = clst_zmat, read_msg = 'ZMAT_from_R={}'.format(clst_r) )
         exec_exec = 'Running'

    print_tab( 4, 'make_err_dir ends: {} sec.'.format( time.time() - t5 ) )
    return error_dir, exec_exec

def main():

    global DIM_LABEL

    global ZERO_DFT_ENER
    global ZERO_MP2_ENER

    global R
    global T
    global P

    if len(sys.argv) == 1:
      DIM_LABEL = 'EMIM_BF4'
    else:
      DIM_LABEL = sys.argv[1]

    CAT_LABEL, ANI_LABEL = DIM_LABEL.split('_')
 
    print_tab( 0, '>>>> {} <<<< begins'.format(DIM_LABEL) )
    T0 = time.time()

    opt_basis_list, opt_funct_list = minimal_basis_list, minimal_funct_list

    if DIM_LABEL == 'EMIM_BF4':
       ene_basis_list, ene_funct_list = limited_basis_list, limited_funct_list
    else:
       ene_basis_list, ene_funct_list = minimal_basis_list, minimal_funct_list

    #ene_basis_list, ene_funct_list = minimal_basis_list, minimal_funct_list

    make_opt = False
    make_opt = True

    fix_errors = True 
    fix_errors = False

    make_energy = False 
    make_energy = True

    if make_opt:
       for opt_basis in opt_basis_list:
         for opt_funct in opt_funct_list:

           print_tab(3, '=== OPTIMIZATION: {}, {} === begins'.format(opt_basis, opt_funct) )

           opt_dimer, cat_nat, ani_nat, zero_dft_en, zero_mp2_en = read_dimer( DIM_LABEL, opt_basis, opt_funct )

           os.makedirs( opt_dimer.runs_dir, exist_ok=True )
           os.makedirs( opt_dimer.csv_dir,  exist_ok=True )
           os.makedirs( os.path.join( opt_dimer.csv_dir, opt_basis, opt_funct ), exist_ok = True )

           err_df = pd.DataFrame( columns = err_columns )
           opt_df = pd.DataFrame( columns = opt_columns )

           if os.path.exists( opt_dimer.scan_err_csv ):
              err_df = pd.read_csv( opt_dimer.scan_err_csv, index_col=0 )
           if os.path.exists( opt_dimer.scan_opt_csv ):
              opt_df = pd.read_csv( opt_dimer.scan_opt_csv, index_col=0 )

           print_tab(3, 'len({}) = {}'.format( opt_dimer.scan_opt_csv, len(opt_df) ))
           print_tab(3, 'len({}) = {}'.format( opt_dimer.scan_err_csv, len(err_df) ))

           r_opt_list, t_opt_list, p_opt_list = get_scan_list( DIM_LABEL, opt_basis, opt_funct )
           #r_opt_list = [3.2]
           #t_opt_list = [175]
           #p_opt_list = [90]

           for T in t_opt_list:
             print_tab(3, '=========================================='.format(T))
             for P in p_opt_list:
               print_tab(3, '------------------------------------------'.format(P))
               for R in r_opt_list:
                   print_tab( 3, '{}, {}, T = {}, P = {}, R = {}'.format(opt_basis, opt_funct, T, P, R) )

                   err_line = locate_df_line( err_df, R, T, P )
                   opt_line = locate_df_line( opt_df, R, T, P )

                   if opt_line.empty and err_line.empty:
                      opt_label, opt_obj = get_gms_object( DIM_LABEL, opt_basis, opt_funct, T, P, R )
                      opt_exec, opt_exec_err, opt_gms_err, opt_scf, opt_geom, opt_time = read_object( opt_obj ) 
                      opt_Series = pd.Series( { 'Radius'   : float(R),     'Theta'   : int(T),  'Phi' : int(P), 
                                                'exec.'    : opt_exec,     'scf'     : opt_scf, 'geom' : opt_geom, 
                                                'exec.err' : opt_exec_err, 'gms.err' : opt_gms_err, 'time' : opt_time } )

                      if opt_exec == 'Running':
                         print( 'Running' )
                      elif opt_exec != 'TERMINATED.NORMALLY'  and [ opt_exec_err, opt_gms_err ] == [ False, False ]:
                         print( [ opt_exec, opt_exec_err, opt_gms_err, opt_scf, opt_geom ] )
                      elif [ opt_exec, opt_exec_err, opt_gms_err, opt_scf, opt_geom ] == [ 'TERMINATED.NORMALLY', False, False, 'CONVERGED', 'LOCATED' ]:
                         opt_df = opt_df.append( opt_Series, ignore_index=True )
                      else:
                         err_df = err_df.append( opt_Series, ignore_index=True )

               for [df_obj, csv_obj] in zip( [opt_df, err_df], [opt_dimer.scan_opt_csv, opt_dimer.scan_err_csv] ):
                    df_obj.sort_values(by=['Theta','Phi','Radius'], inplace=True )
                    df_obj.drop_duplicates(inplace=True)
                    df_obj.reset_index(drop=True, inplace=True)
                    df_obj.to_csv( csv_obj )
                    print_tab( 1, 'printing {} (len={})'.format( csv_obj, len(df_obj) ))

               ### CORRECT OPTIMIZATION
               if fix_errors:
                  print_tab(2, '=== FIX ERRORS === begins')

                  print_tab( 3, '{}, {}, T = {}, P = {}'.format(opt_basis, opt_funct, T, P) )
                  succ_block = opt_df.loc[ opt_df['Theta']==float(T) ].loc[ opt_df['Phi']==float(P) ]
                  fail_block = err_df.loc[ err_df['Theta']==float(T) ].loc[ err_df['Phi']==float(P) ]
                  succ_rad = [ float(x) for x in succ_block['Radius'].values ]
                  fail_rad = [ float(x) for x in fail_block['Radius'].values ]
                  succ_rad.sort()
                  fail_rad.sort()
                  print_tab( 3, '{} suceeded radius = {}\n'.format(len(succ_rad), succ_rad) )
                  print_tab( 3, '{} failed radius = {}\n'.format(len(fail_rad),fail_rad) )

                  # check for duplicates 
                  if [ r for r in succ_rad if r in fail_rad ] != []:
                       raise NameError( 'duplicate Radius in opt and fail dataframes\nplease check df below:\n{}\n{}'.format(succ_block, fail_block) )
 
                  ## make copy of extended_R_SCAN_LIST
                  all_R = [ float(r) for r in extended_R_SCAN_LIST ]
                  all_R.sort()

                  if not fail_block.empty:
                     fail_block.sort_values( by = 'Radius', ascending=False, inplace=True )
                     print( fail_block )
                     submitted_new = False
                     for fail_idx, fail_row in fail_block.iterrows():
                         fail_R, fail_exec_err, fail_geom, fail_gms_err, fail_scf, fail_exec = read_df_line( fail_row )
                         #if any( [ err in [ 'insufficient.distributed.memory', 'insufficient.replicated.memory', 
                         #                   'memory.request.exceeds.available.memory','not.enough.replicated.memory'] \
                         #          for err in [ fail_exec_err, fail_gms_err ] ] ) :
                         #   if fail_exec_err != fail_gms_err:
                         #      print( cannot_understand_error )
                         #   fail_obj = get_gms_object( DIM_LABEL, opt_basis, opt_funct, T, P, fail_R )[1]
                         #   fail_obj.fix_memory_error( fail_exec_err )
                         if submitted_new == False:
                            print_tab( 3, '>>>> trying to fix R = {}:'.format(fail_R) )
                            fail_obj = get_gms_object( DIM_LABEL, opt_basis, opt_funct, T, P, fail_R )[1]
                            fail_exec, fail_exec_err, fail_gms_err, fail_scf, fail_geom, fail_time = read_object( fail_obj ) 

                            ## make array of distances
                            all_R.remove(fail_R)

                            if fail_exec == 'Running':
                               print_tab( 3, '>>>> R = {} is Running, skipping'.format(fail_R) )
                               submitted_new = True
                            else:
                               distances_list = iter(R_distances( all_R, fail_R ))
                               while not submitted_new:
                                     clst_r, clst_d = next(distances_list)
                                     print_tab( 4, 'closest_R = {} (distance = {})'.format(clst_r, clst_d) )
                                     tmp_succ_line = succ_block.loc[ succ_block['Radius'] == float(clst_r) ] 
                                     tmp_fail_line = fail_block.loc[ fail_block['Radius'] == float(clst_r) ] 

                                     if tmp_succ_line.empty and tmp_fail_line.empty:
                                        print_tab( 5, 'closest_R = {} is neither finished nor un-finished'.format(clst_r) )
                                     elif tmp_succ_line.empty and not tmp_fail_line.empty:
                                        print_tab( 5, 'closest_R = {} is un-finished, will skip ... '.format(clst_r) )
                                     elif not tmp_succ_line.empty and tmp_fail_line.empty:
                                        print_tab( 5, 'closest_R = {} is finished, will copy zmat from there ... '.format(clst_r) )
                                        error_dir, exec_exec = make_err_dir( fail_obj, fail_R, clst_r ) 
                                        err_df.at[fail_idx, 'exec.'] = exec_exec 
                                        if exec_exec == 'Running':
                                           submitted_new = True
                                     else:
                                        raise NameError( 'unknown status in df' )
                            print_tab( 3, '>>>> fix finished.' )
                  print_tab(2, '--- FIX ERRORS --- ends')
           print_tab(3, '--- OPTIMIZATION: {}, {} --- ends'.format(opt_basis, opt_funct) )

    ##############
    ## MAKE SINGLE POINT ENERGY start
    if make_energy: 
       rlx_basis = 'N311'
       rlx_funct = 'B3LYP'

       for ene_basis in ene_basis_list:
           single_point_calculations( DIM_LABEL, rlx_basis, rlx_funct, ene_basis, 'B3LYP', post_scf = 'MP2' )
           for ene_funct in ene_funct_list:
               single_point_calculations( DIM_LABEL, rlx_basis, rlx_funct, ene_basis, ene_funct, post_scf = 'DFTTYP' )
    ## MAKE SINGLE POINT ENERGY ends

    T1 = time.time()
    print_tab( 0, '>>>> {} <<<< ends'.format(DIM_LABEL) )
    print_tab( 0, 'TOTAL TIME = {} seconds'.format(T1-T0) )

if __name__ == '__main__':
   main()

