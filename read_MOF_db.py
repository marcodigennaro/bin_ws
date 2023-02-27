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

def write_best_mono_df( mono_results_df, best_mono_csv ):
    if os.path.exists( best_mono_csv ):
       best_mono_df = pd.read_csv( best_mono_csv, index_col = 0 )
    else:
       best_mono_df = pd.DataFrame()
       for fg, block in mono_results_df.groupby('FG'):
           block_min_ene = block['MP2.VT.EN.'].min()
           block_min_line = block.loc[ block['MP2.VT.EN.'] == block_min_ene ]
           best_mono_df = best_mono_df.append( block_min_line, ignore_index = True )
       best_mono_df.to_csv( best_mono_csv )
    print( 'Best mono csv file: ', best_mono_csv )
    return best_mono_df

def write_best_dime_df( dime_results_df, best_dime_csv ):
    if os.path.exists( best_dime_csv):
       best_dime_df = pd.read_csv( best_dime_csv, index_col = 0 )
    else:
       best_dime_df = pd.DataFrame()
       if not dime_results_df.empty:
          for fg, block in dime_results_df.groupby('FG'):
              block_min_ene = block['MP2.VT.INT.EN.'].min()
              block_min_line = block.loc[ block['MP2.VT.INT.EN.'] == block_min_ene ]
              best_dime_df = best_dime_df.append( block_min_line, ignore_index = True )
       best_dime_df.to_csv( best_dime_csv )
    print( 'Best dime csv file: ', best_dime_csv )
    return best_dime_df

def main():
 
    orca_dir =  '/home/rgi2972/runs_orca/fg/'
    fg_csv = '/home/rgi2972/INPUTFILES/ORCA/FG.csv' 
    fg_df = pd.read_csv( fg_csv, index_col=0 )
    fg_int_cols   = ['dime.idx','mono.idx','mono.nat']
    mono_int_cols = ['mono.idx','mono.nat']
    dime_int_cols = ['NEG.FREQ.','VERY.NEG.FREQ.','dime.idx','mono.idx','mono.nat']
    fg_df[fg_int_cols] = fg_df[fg_int_cols].astype(int)

    ##### options to run #####
    check_mono     = False 
    check_dime     = False 
    run_bsse       = False 
    run_rigid_scan = False 
    submit_new     = False 

    check_mono     = True
    check_dime     = True
    run_bsse       = True
    run_rigid_scan = True
    #submit_new     = True

    csv_dir = '/home/rgi2972/ORCA_CSV'
    csv_dir = '/data/mdi0316/WORK/MOFS/ORCA_CSV'
    os.makedirs( csv_dir, exist_ok = True )

    #tmp_fgs = [ '66_OSO2NH2', '67_SO2CN','68_CONHCN','69_CH2SO2H','70_COCHNN','71_CH2CN','72_SPO3H2','73_PSOH2','74_SO2NHOH','75_SO2NHNH2' ]
    #tmp_fgs += [ '33_OCONH2' ]
    #tmp_fgs = [ '66_OSO2NH2', '67_SO2CN' ]
    #tmp_fgs = [ '33_OCONH2', '75_SO2NHNH2' ]
    #fg_df = fg_df[ fg_df['FG'].isin( tmp_fgs ) ]

    ##### options to run #####

    mono_status_csv = os.path.join( csv_dir, 'mono_status.csv'  ) 
    mono_results_csv = os.path.join( csv_dir, 'mono_results.csv'  ) 
    best_mono_csv = os.path.join( csv_dir, 'best_mono.csv' )

    dime_status_csv = os.path.join( csv_dir, 'dime_status.csv'  ) 
    dime_results_csv = os.path.join( csv_dir, 'dime_results.csv' )
    best_dime_csv = os.path.join( csv_dir, 'best_dime.csv' )
 
    rigid_scan_csv = os.path.join( csv_dir, 'rigid_scan.csv' )

    ### READ MONOMERS starts
    if check_mono:
       print( 24*'#', '###  MONO  (starts)  ###',  24*'#', sep='\n' )

       print('reading ', mono_results_csv )
       if os.path.exists( mono_results_csv ):
          mono_status_df = pd.read_csv( mono_status_csv, index_col = 0 )
          mono_results_df = pd.read_csv( mono_results_csv, index_col = 0 )
          mono_results_df[mono_int_cols] = mono_results_df[mono_int_cols].astype(int)
       else:
          print('New dataframe')
          mono_results_df = pd.DataFrame()
          mono_status_df = pd.DataFrame()

       for fg_row, fg_line in fg_df.iterrows():
           mono_lab = fg_line['FG']
           mono_nat = fg_line['mono.nat']

           if mono_results_df.empty:
              mono_line = pd.DataFrame()
           else:
              mono_line = mono_results_df.loc[ ( mono_results_df['FG']       == mono_lab ) & \
                                               ( mono_results_df['mono.nat'] == mono_nat ) ]

           if mono_line.empty:
              mono_idx = fg_line['mono.idx']
              for tmp_idx in range(1, mono_idx+1):
                  mono_obj = ORCA.MONO( mono_lab, mono_nat, tmp_idx )
                  # tight
                  mono_obj.files_names( mono=True )
                  mono_obj.read_output_status()
                  if mono_obj.status == 'output.terminated.normally':
                     mono_mp2_ene = mono_obj.read_output_file()
                     # very tight
                     mono_obj.files_names( mono=True, very_tight=True )
                     mono_obj.read_output_status()
                     if mono_obj.status == 'output.terminated.normally':
                        mono_mp2_vt_ene, mono_neg_freq, mono_vneg_freq = mono_obj.read_output_file( very_tight=True )
                        mono_results_df = mono_results_df.append( { 'FG'             : mono_lab,
                                                                    'mono.idx'       : int(tmp_idx),
                                                                    'mono.nat'       : int(mono_nat),
                                                                    'MP2.EN.'        : float(mono_mp2_ene),
                                                                    'MP2.VT.EN.'     : float(mono_mp2_vt_ene), 
                                                                    'NEG.FREQ.'      : int(mono_neg_freq), 
                                                                    'VERY.NEG.FREQ.' : int(mono_vneg_freq) }, 
                                                                    ignore_index = True )
                     else:
                       print( mono_obj.status, mono_obj.out_file )
                  else:
                    print( mono_obj.status, mono_obj.out_file )

              mono_status_df = mono_status_df.append( { 'FG' : mono_lab, 'mono.idx' : int(tmp_idx), 
                                                        'status' : mono_obj.status, 'out.file': mono_obj.out_file }, ignore_index = True )
              
       mono_results_df.sort_values( by='FG', inplace=True )
       mono_results_df.reset_index( drop=True, inplace=True )
       print( 'Printing: ', mono_results_csv )
       mono_results_df.to_csv( mono_results_csv )
       mono_status_df.to_csv( mono_status_csv )

       ## Extract monomer conformation with lowest energy
       best_mono_df = write_best_mono_df( mono_results_df, best_mono_csv )
       print( 24*'-', '---   MONO  (ends)   ---',  24*'-', sep='\n' )

    else:
       mono_results_df = pd.read_csv( mono_results_csv, index_col = 0 )
    ### READ MONOMERS ends


    ### READ DIMERS starts
    if os.path.exists( dime_results_csv ):
       dime_status_df = pd.read_csv( dime_status_csv, index_col = 0 )
       dime_results_df = pd.read_csv( dime_results_csv, index_col = 0 )
       dime_results_df[dime_int_cols] = dime_results_df[dime_int_cols].astype(int)
    else:
       dime_status_df = pd.DataFrame()
       dime_results_df = pd.DataFrame()

    if check_dime:
       print( 24*'#', '###  DIME  (starts)  ###',  24*'#', sep='\n' )

       best_mono_df = write_best_mono_df( mono_results_df, best_mono_csv )
       for best_mono_row, best_mono_line in best_mono_df.iterrows():
           fg_label = best_mono_line['FG']
           best_mono_idx = int(best_mono_line['mono.idx'])
           best_mono_ene = best_mono_line['MP2.EN.']
           best_mono_vt_ene = best_mono_line['MP2.VT.EN.']
          
           if fg_label != '00_H2':
              print_tab( 1, '===============' )
              print_tab( 1, 'Reading {} '.format( fg_label ) )
              print_tab( 1, '---------------' )

              fg_line = fg_df.loc[ fg_df['FG'] == fg_label ]
              fg_idx = fg_label.split('_')[0]
              fg_nat = fg_line['mono.nat'].values[0]
              dime_idx = fg_line['dime.idx'].values[0]
        
              for dime_idx in range( 1, dime_idx + 1 ):
                  try:
                     print_tab( 2, '>> FG{}/M{}/D{}/{}'.format(fg_idx, best_mono_idx, best_mono_idx, dime_idx) )
                     if dime_results_df.empty:
                        dime_line = pd.DataFrame()
                     else:
                        dime_line = dime_results_df.loc[ ( dime_results_df['FG'] == fg_label) & \
                                                     ( dime_results_df['mono.idx'] == best_mono_idx ) & \
                                                     ( dime_results_df['dime.idx'] == dime_idx ) ]
                     if dime_line.empty:
                        dime_obj = ORCA.DIME( fg_label, fg_nat, best_mono_idx, dime_idx )
                        # tight
                        dime_obj.files_names( dime=True )
                        dime_obj.read_output_status()
                        if dime_obj.status == 'output.terminated.normally':
                           dime_mp2_ene = dime_obj.read_output_file()
                           # very tight
                           dime_obj.files_names( dime=True, very_tight=True )
                           dime_obj.read_output_status()
                           if dime_obj.status == 'output.terminated.normally':
                              dime_mp2_vt_ene, dime_neg_freq, dime_vneg_freq = dime_obj.read_output_file( very_tight=True )
                              # write energy + coordinates results
                              dime_dict =  { 'FG'             : fg_label, 
                                             'mono.idx'       : int(best_mono_idx), 
                                             'dime.idx'       : int(dime_idx),
                                             'mono.nat'       : int(fg_nat),
                                             'MP2.EN.'        : dime_mp2_ene, 
                                             'MP2.INT.EN.'    : dime_mp2_ene - best_mono_ene, 
                                             'MP2.VT.EN.'     : dime_mp2_vt_ene, 
                                             'MP2.VT.INT.EN.' : dime_mp2_vt_ene - best_mono_vt_ene, 
                                             'NEG.FREQ.'      : dime_neg_freq, 
                                             'VERY.NEG.FREQ.' : dime_vneg_freq } 
                              dime_vt_coords = dime_obj.read_coordinates()
                              for count, coord in enumerate(dime_vt_coords):
                                  k, elem, x, y, z = coord
                                  dime_dict[k] = { 'elem.' : elem, 'idx.' : count, 'x' : x, 'y' : y, 'z' : z }
                              dime_results_df = dime_results_df.append( [dime_dict], ignore_index = True )
                           else:
                             dime_dict =  { 'FG'             : fg_label, 
                                            'mono.idx'       : int(best_mono_idx), 
                                            'dime.idx'       : int(dime_idx),
                                            'mono.nat'       : int(fg_nat),
                                            'MP2.EN.'        : dime_mp2_ene, 
                                            'MP2.INT.EN.'    : dime_mp2_ene - best_mono_ene } 
                             dime_results_df = dime_results_df.append( [dime_dict], ignore_index = True )
                             if submit_new:
                                dime_obj.write_very_tight_input()
                        elif dime_obj.status in ORCA.error_dict.values():
                          print( dime_obj.status, dime_obj.out_file )
                        else:
                          if submit_new:
                             dime_obj.write_tight_input()
                  except(UnicodeDecodeError):
                    dime_obj.status = 'UnicodeDecodeError' 
                  dime_status_df = dime_status_df.append( { 'FG' : fg_label, 'mono.idx' : int(best_mono_idx), 'dime.idx' : int(dime_idx),
                                                            'status' : dime_obj.status, 'out.file' : dime_obj.out_file }, ignore_index = True )
                        
       print( 'Printing: ', dime_results_csv )
       print( dime_results_df )
       if not dime_results_df.empty:
          dime_results_df.sort_values( by='FG', inplace=True )
          dime_results_df.reset_index( drop=True, inplace=True )
       dime_results_df.to_csv( dime_results_csv )
       dime_status_df.to_csv( dime_status_csv )

       ### get best DIMERS start
       best_dime_df = write_best_dime_df( dime_results_df, best_dime_csv )
    ### get best DIMERS ends
    ### READ DIMERS ends

    ### BSSE correction starts
    copy_dime_df = pd.DataFrame( dime_results_df )
    if run_bsse:
       print( 24*'#', '###  BSSE  (starts)  ###',  24*'#', sep='\n' )
       for fg_r, fg_val in dime_results_df.iterrows():
           fg_label = fg_val['FG']
           fg_nat = int(fg_val['mono.nat'])
           fg_mono_idx = int(fg_val['mono.idx'])
           fg_dime_idx = int(fg_val['dime.idx'])
           fg_very_neg_freq = int(fg_val['VERY.NEG.FREQ.'])
           print_tab( 3, '>> FG{}/M{}/D{}/{}'.format(fg_label, fg_mono_idx, fg_mono_idx, fg_dime_idx) )
           if fg_very_neg_freq == 0:
              fg_obj = ORCA.DIME( fg_label, fg_nat, fg_mono_idx, fg_dime_idx )
              fg_obj.files_names( dime=True, very_tight=True )
              ## now select 
              cp_obj = ORCA.DIME( fg_label, fg_nat, fg_mono_idx, fg_dime_idx )
              cp_obj.files_names( dime=True, very_tight=True, counterpoise=True )
              if os.path.exists( cp_obj.counterpoise_run_dir ):
                 if os.path.exists( cp_obj.out_file ):
                    cpc = cp_obj.read_counterpoise()
                    if cpc:
                       copy_dime_df.at[fg_r, 'CP.CORR.'] = cp_obj.read_counterpoise() 
                    else:
                       print( 'WARNING: CPC failed in {}'.format(cp_obj.run_dir) )
                 else:
                    print( 'missing out file in ', cp_obj.counterpoise_run_dir )
              else:
                 if submit_new:
                    xyz_coordinates = fg_obj.read_coordinates()
                    cp_obj.write_counterpoise_file( xyz_coordinates )
           else:
             print( 'skipping CP since very negative frequency exists' )

       print( 24*'-', '---  BSSE (starts)  ---',  24*'-', sep='\n' )

    print( 'Printing: ', dime_results_csv )
    copy_dime_df.sort_values( by='FG', inplace=True )
    copy_dime_df.reset_index( drop=True, inplace=True )
    copy_dime_df.to_csv( dime_results_csv )
    ### BSSE correction ends

    ### rigid scan
    if run_rigid_scan:
       print( 27*'#', '##  RIGID SCAN (starts)  ##',  27*'#', sep='\n' )
       rigid_scan_idx_df = fg_df.loc[ fg_df['rigid.scan'] != 0 ]
       #best_dime_df = pd.read_csv( best_dime_csv, index_col = 0 )
       best_dime_df = write_best_dime_df( dime_results_df, best_dime_csv )
       rigid_scan_df = pd.DataFrame()
       for fg_idx, fg_read_line in rigid_scan_idx_df.iterrows():
           fg_label     = fg_read_line['FG']
           fg_rigid_idx = fg_read_line['rigid.scan']
           fg_nat       = fg_read_line['mono.nat']

           best_fg_dict = best_dime_df.loc[ best_dime_df['FG'] == fg_label ].to_dict('list')
           mono_fg_idx  = best_fg_dict['mono.idx'][0]
           dime_fg_idx  = best_fg_dict['dime.idx'][0]
           vneg_fg_freq = best_fg_dict['VERY.NEG.FREQ.'][0]

           for rigid_idx in range( 1, fg_rigid_idx+1 ): 
               print_tab( 3, '>> FG{}/M{}/D{}/{}'.format(fg_label, mono_fg_idx, mono_fg_idx, dime_fg_idx) )

               best_mono_dict = best_mono_df.loc[ best_mono_df['FG'] == fg_label ].to_dict('list')
               best_mono_vt_ene = float( best_mono_dict['MP2.VT.EN.'][0] )

               if vneg_fg_freq == 0:
     
                  ## rigid scan done by hand
                  rs_obj = ORCA.DIME( fg_label, fg_nat, mono_fg_idx, dime_fg_idx )
                  rs_obj.files_names( rigid_scan=rigid_idx, very_tight=True )
                  rs_trj = rs_obj.read_trajectories()
                  trj_dict = rs_obj.read_trj_file()
  
                  ## counterpoise object  
                  #cp_obj = ORCA.DIME( fg_label, fg_nat, fg_mono_idx, fg_dime_idx )
                  #cp_obj.files_names( rigid_scan=rigid_idx, counterpoise=True )

                  for kk, vv in rs_trj.items():
                      radius = vv['RADIUS'] 
                      radius_obj = ORCA.DIME( fg_label, fg_nat, mono_fg_idx, dime_fg_idx )
                      radius_obj.files_names( rigid_scan=rigid_idx, counterpoise=True, radius=radius )
                      if os.path.exists( radius_obj.counterpoise_run_dir ):
                        if os.path.exists( radius_obj.out_file ):
                           radius_cp_corr = radius_obj.read_counterpoise()
                           radius_dict = { 'FG' : fg_label, 'Radius' : radius, 'rigid.scan' : rigid_idx, 
                                           'MP2.EN.'     : trj_dict[radius], 
                                           'MP2.INT.EN.' : float(trj_dict[radius]) - float(best_mono_vt_ene),
                                           'CP.CORR.'    : radius_cp_corr }
                           rigid_scan_df = rigid_scan_df.append( [radius_dict], ignore_index=True )
                        else:
                           print_tab( 4, 'missing out file in ', radius_obj.counterpoise_run_dir )
                      else:
                         if submit_new:
                            xyz_coordinates = rs_obj.read_trajectories()
                            radius_coordinates =  [ cc_v for cc_v in xyz_coordinates.values() \
                                                 if cc_v['RADIUS'] == radius ][0]['CART.COORDS.']
                            radius_obj.write_counterpoise_file( radius_coordinates )
               else:
                 print_tab( 4, 'very negative frequencies found' )
 
       print( 27*'-', '--  RIGID SCAN (starts)  --',  27*'-', sep='\n' )
       print( 'Printing: ', rigid_scan_csv )
       rigid_scan_df.to_csv( rigid_scan_csv )
    ### rigid scan


if __name__ == '__main__':
  
  main()
