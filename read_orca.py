#!/home/mdi0316/anaconda3/bin/python

import os, sys
import pandas as pd
import numpy as np
import warnings
import subprocess as sp
import shutil

class_dir = '/home/mdi0316/CLASSES'
sys.path.insert( 0, class_dir )

funct_dir = '/home/mdi0316/FUNCTIONS'
sys.path.insert( 0, funct_dir )

import ORCA
from Functions import print_tab

def read_obj( orca_obj ):

    mp2_ene, mp2_vt_ene, neg_freq, very_neg_freq = 4*[float('Nan')]

    tight_calculation = False
    if os.path.exists( orca_obj.inp_file ):
       if os.path.exists( orca_obj.out_file ):
          mp2_ene = orca_obj.read_output()
          if not np.isnan( mp2_ene ):
             tight_calculation = True
       else:
          print_tab( 2, 'tight calculation output MISSING' )

    if tight_calculation:
       if os.path.exists( orca_obj.vt_inp_file ):
          if os.path.exists( orca_obj.vt_out_file ):
             mp2_vt_ene, neg_freq, very_neg_freq = orca_obj.read_output( very_tight = True )
          else:
             print_tab( 2, 'very tight calculation output MISSING' )
       else:
          print( 'Make new vt calculation' )
          orca_obj.write_very_tight_input()

    if any( [ ee == 0. for ee in [ mp2_ene, mp2_vt_ene ] ] ):
       print( 'WARNING one energy is zero' )

    return( mp2_ene, mp2_vt_ene, neg_freq, very_neg_freq )

def main():
 
    ### options to run ###
    check_mono = True   ##
    check_mono = False  ##
                        ## 
    check_dime = False  ##
    check_dime = True   ##
                        ##
    run_bsse   = False  ##
    run_bsse   = True   ##
    ### options to run ###

    orca_dir =  '/home/rgi2972/runs_orca/fg/'
  
    fg_csv = '/home/rgi2972/INPUTFILES/ORCA/FG.csv' 
    fg_df = pd.read_csv( fg_csv, index_col=0 )

    #tmp_fg = '42_SO3H'
    #fg_df = fg_df.loc[ fg_df['FG'] == tmp_fg ]

    csv_dir = '/home/mdi0316/ORCA_CSV'
    csv_dir = '/home/rgi2972/ORCA_CSV'
    os.makedirs( csv_dir, exist_ok = True )
    all_mono_csv =  os.path.join( csv_dir, 'all_mono_energies.csv'  ) 
    best_mono_csv = os.path.join( csv_dir, 'best_mono_energies.csv' )

    all_dime_csv = os.path.join( csv_dir, 'all_dime_energies.csv' )
    best_dime_csv = os.path.join( csv_dir, 'best_dime_energies.csv' )
    incomplete_csv = os.path.join( csv_dir, 'incomplete.csv' )

    ### READ MONOMERS starts
    if check_mono:
       # Read output
       if os.path.exists( all_mono_csv ):
          all_mono_df = pd.read_csv( all_mono_csv, index_col = 0 )
       else:
          all_mono_df = pd.DataFrame()
          print( 'writing all_mono_df' )
          for fg_idx, fg_val in fg_df.iterrows():
             fg_k = fg_val['FG']
             fg_v = fg_val['mono.idx']
             fg_nat = int(fg_val['mono.nat'])
             funct = ORCA.ORCA( fg_k, fg_nat )
             for mono_idx in range(1, int(fg_v)+1):
                 mono_obj = ORCA.MONO( fg_k, fg_nat, mono_idx )
                 mono_mp2_ene, mono_mp2_vt_ene, neg_freq, very_neg_freq = read_obj( mono_obj )
                 all_mono_df = all_mono_df.append( { 'FG'         : fg_k, 
                                                     'mono.idx.'  : int(mono_idx), 
                                                     'MP2.EN.'    : float(mono_mp2_ene),
                                                     'MP2.VT.EN.' : float(mono_mp2_vt_ene),}, 
                                                     ignore_index = True )
          all_mono_df.to_csv( all_mono_csv )

       # Extract monomer conformation with lowest energy
       if os.path.exists( best_mono_csv ):
          best_mono_df = pd.read_csv( best_mono_csv, index_col = 0 )
       else:
          best_mono_df = pd.DataFrame()
          for fg_idx, fg_val in fg_df.iterrows():
              fg_k = fg_val['FG']
              fg_tmp_df = all_mono_df.loc[ all_mono_df['FG'] == fg_k ]
              min_ene = all_mono_df.loc[ all_mono_df['FG'] == fg_k ]['MP2.EN.'].min()
              min_vt_ene = all_mono_df.loc[ all_mono_df['FG'] == fg_k ]['MP2.VT.EN.'].min()
              min_vt_line = all_mono_df.loc[ (all_mono_df['FG'] == fg_k) & \
                                             #(all_mono_df['MP2.EN.'] == min_ene ) \
                                             (all_mono_df['MP2.VT.EN.'] == min_vt_ene ) \
                                           ].head(1).to_dict( orient='index')
              min_dict = list( min_vt_line.values() )[0] 
              best_mono_df = best_mono_df.append( [min_dict] , ignore_index = True )
          best_mono_df.to_csv( best_mono_csv )

    else:
       all_mono_df = pd.read_csv( all_mono_csv, index_col = 0 )
       best_mono_df = pd.read_csv( best_mono_csv, index_col = 0 )
    ### READ MONOMERS ends

    ### READ DIMERS starts
    incomplete_df = pd.DataFrame()
    # Read output
    if os.path.exists( all_dime_csv ):
       all_dime_df = pd.read_csv( all_dime_csv, index_col = 0 )
    else:
       all_dime_df = pd.DataFrame()

    if check_dime:
       for fg_idx, fg_val in fg_df.iterrows():
           if 1 < fg_idx < 46 and fg_idx != 33:
           #if fg_idx not in [ 0, 33] and fg_idx < 60:
              fg_k = fg_val['FG']
              fg_d = int(fg_val['dime.idx'])
              fg_nat = int(fg_val['mono.nat'])
              print_tab( 1, '===============' )
              print_tab( 1, 'Reading {} '.format( fg_k ) )
              print_tab( 1, '---------------' )
               
              func_gr_obj = ORCA.ORCA( fg_k, fg_nat )
              if len( [ ff for ff in os.listdir() if 'restart' in ff ] ):
                 print( 'warning restart file exists\nneed to check' )
              max_dime_idx, max_dime_dir = func_gr_obj.max_dimer_dir() 

              best_mono_dict = best_mono_df.loc[ best_mono_df['FG'] == fg_k ].to_dict('list')
              best_mono_idx = int( best_mono_dict['mono.idx.'][0] )
              best_mono_ene = float( best_mono_dict['MP2.EN.'][0] )
              best_mono_vt_ene = float( best_mono_dict['MP2.VT.EN.'][0] )

              ## check starts
              ## checking that best mono index (ie mono with lowest energy) 
              ## corresponds to highest dimer index found in fg_dir 
              if int( max_dime_idx ) != int( best_mono_idx ):
                 print( 'WARNING: max_dime_idx ({}) != best_mono_idx ({})'.format(max_dime_idx, best_mono_idx))
                 pd.options.display.float_format = '{:.10f}'.format
                 all_mono_tmp_dict = all_mono_df.loc[ all_mono_df['FG'] == fg_k ].sort_values( by = ['MP2.VT.EN.'] )
              ## check ends

              for dime_idx in range( 1, fg_d + 1 ):
                  print_tab( 2, '>> FG{}/M{}/D{}/{}'.format(fg_idx, best_mono_idx, best_mono_idx, dime_idx) )
                  if all_dime_df.empty:
                     dime_line = pd.DataFrame()
                  else:
                     dime_line = all_dime_df.loc[ (all_dime_df['FG'] == fg_k) & \
                                                  (all_dime_df['mono.idx.'] == best_mono_idx ) & \
                                                  (all_dime_df['dime.idx.'] == dime_idx ) ]
                  if dime_line.empty:
                     dime_obj = ORCA.DIME( fg_k, fg_nat, best_mono_idx, dime_idx )
                     dime_mp2_ene, dime_mp2_vt_ene, neg_freq, very_neg_freq = read_obj( dime_obj )
                     if dime_mp2_ene == False or np.isnan( dime_mp2_ene ) or np.isnan( dime_mp2_vt_ene ) :
                        print( 'WARNING >>>> come and fix me!!!')
                        incomplete_df = incomplete_df.append( 
                        { 'FG' : fg_k, 'mono.idx.'  : best_mono_idx, 
                                       'dime.idx.'  : dime_idx ,
                                       'MP2.EN.'    : dime_mp2_ene, 
                                       'MP2.VT.EN.' : dime_mp2_vt_ene },
                        ignore_index = True )
                     else:
                        dime_line =  { 'FG' : fg_k, 'mono.idx.': best_mono_idx, 'dime.idx.' : dime_idx, 
                                       'mono.nat.'      : int(fg_nat),
                                       'MP2.EN.'        : dime_mp2_ene, 
                                       'MP2.INT.EN.'    : dime_mp2_ene - best_mono_ene, 
                                       'MP2.VT.EN.'     : dime_mp2_vt_ene, 
                                       'MP2.VT.INT.EN.' : dime_mp2_vt_ene - best_mono_vt_ene, 
                                       'NEG.FREQ.'      : int(neg_freq), 
                                       'VERY.NEG.FREQ.' : int(very_neg_freq) } 
                        all_dime_df = all_dime_df.append( dime_line, ignore_index = True )
       all_dime_df.to_csv( all_dime_csv )
       incomplete_df.to_csv( incomplete_csv )
    ### READ DIMERS ends

    ### BSSE correction starts
    copy_dime_df = pd.DataFrame( all_dime_df )
    if run_bsse:
       print_tab( 1, 'BSSE' )
       for fg_r, fg_val in all_dime_df.iterrows():
           fg_kk = fg_val['FG']
           fg_nat = int(fg_val['mono.nat.'])
           fg_mono_idx = int(fg_val['mono.idx.'])
           fg_dime_idx = int(fg_val['dime.idx.'])
           fg_very_neg_freq = int(fg_val['VERY.NEG.FREQ.'])
           print_tab( 3, '>> FG{}/M{}/D{}/{}'.format(fg_kk, fg_mono_idx, fg_mono_idx, fg_dime_idx) )
           print(fg_val)
           if fg_very_neg_freq == 0:
              fg_obj = ORCA.DIME( fg_kk, fg_nat, fg_mono_idx, fg_dime_idx )
              print( fg_obj )
              cp_obj = ORCA.COUNTERPOISE( fg_kk, fg_nat, fg_mono_idx, fg_dime_idx )
              if os.path.exists( cp_obj.run_dir ):
                 copy_dime_df.at[fg_r, 'CP.CORR.'] = cp_obj.read_output( counterpoise = True )
              else:
                 xyz_coordinates = fg_obj.read_coordinates( very_tight = True )
                 cp_obj.write_counterpoise_file( xyz_coordinates )
           else:
             print( 'skipping CP since very negative frequency exists' )
    copy_dime_df.to_csv( all_dime_csv )
    ### BSSE correction ends

    ### rigid scan
        #all_coordinates = fg_obj.read_output()
        #print( all_coordinates ) 
        #traj_dict = fg_obj.read_scan_trajectories()
        #exit()
        #for kk, vv in all_coordinates.items():
        #    kk_rad = vv['RADIUS'] 
        #    vv['SP.MP2.EN.'] = traj_dict[kk_rad] 
        #    kk_dir = os.path.join( fg_obj.counterpoise_run_dir, 'R_{}'.format( kk_rad ) )
        #    if os.path.exists( kk_dir ):
        #       cp_corr = cp_obj.read_output( counterpoise = True )
        #       print( fg_kk, kk_rad, cp_corr )
        #    else:
        #       print(fasjlfsp)
        #       pass
    ### rigid scan

    ### get best DIMERS start
    # Extract monomer conformation with lowest energy
    all_dime_df = all_dime_df.loc[ all_dime_df['NEG.FREQ.'] == 0. ] 
    if os.path.exists( best_dime_csv ):
       best_dime_df = pd.read_csv( best_dime_csv, index_col = 0 )
    else:
       best_dime_df = pd.DataFrame()
       for fg_idx, fg_val in fg_df.iterrows():
           fg_k = fg_val['FG']
           fg_tmp_df = all_dime_df.loc[ all_dime_df['FG'] == fg_k ]
           if not fg_tmp_df.empty:
              min_vt_ene = all_dime_df.loc[ all_dime_df['FG'] == fg_k ]['MP2.VT.EN.'].min()
              min_vt_line = all_dime_df.loc[ (all_dime_df['FG'] == fg_k) & \
                                             (all_dime_df['MP2.VT.EN.'] == min_vt_ene ) \
                                           ].head(1).to_dict(orient='index')
              min_dict = list( min_vt_line.values() )[0] 
#              del min_dict['index']
              best_dime_df = best_dime_df.append( [min_dict] , ignore_index = True )
       best_dime_df.to_csv( best_dime_csv )
    ### get best DIMERS ends

if __name__ == '__main__':
  
  main()
