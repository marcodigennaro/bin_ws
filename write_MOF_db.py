#!/home/mdi0316/anaconda3/bin/python

import os, sys
import pandas as pd
import numpy as np
import math

class_dir = '/home/mdi0316/CLASSES'
sys.path.insert( 0, class_dir )

funct_dir = '/home/mdi0316/FUNCTIONS'
sys.path.insert( 0, funct_dir )

import ORCA
from Functions import print_tab, rotate_molecule, rotation_matrix
import mendeleev

def initiate_elem_dict( symbol, idx, coordinates, rotated_coordinates, group ): 
    elem = mendeleev.element( symbol )
    idx_dict =  { 'element'       : elem.symbol,
                  'elec_affinity' : elem.electron_affinity,
                  'ioniz_energy'  : elem.ionenergies[1], 
                  'group'         : elem.group_id,
                  'period'        : elem.period,
                  'vdw_radius'    : elem.vdw_radius, 
                  'atomic_weight' : elem.atomic_weight,
                  'xx'            : rotated_coordinates[idx][0],
                  'yy'            : rotated_coordinates[idx][1],
                  'zz'            : rotated_coordinates[idx][2],
                  'BENZ.RING'     : 0,
                  'FUN.GROUP'     : 0,
                  'H2.MOL'        : 0 }

    if   group == 'benz.ring': 
         idx_dict['BENZ.RING'] = 1
    elif group == 'fun.group': 
         idx_dict['FUN.GROUP'] = 1
    elif group == 'h2.mol': 
         idx_dict['H2.MOL']    = 1
    
    return idx_dict

def main():

    orca_dir =  '/home/rgi2972/runs_orca/fg/'

    fg_csv = '/home/rgi2972/INPUTFILES/ORCA/FG.csv'
    fg_df = pd.read_csv( fg_csv, index_col=0 )

    #tmp_fgs = ['04_CH3', '11_COOH', '42_SO3H']
    #fg_df = fg_df[ fg_df['FG'].isin( tmp_fgs ) ]

    csv_dir = '/home/rgi2972/ORCA_CSV'
    csv_dir = '/data/mdi0316/WORK/MOFS/ORCA_CSV'
    dime_results_csv = os.path.join( csv_dir, 'dime_results.csv' )
    dime_results_df = pd.read_csv( dime_results_csv, index_col = 0 )
    
    xyz_coords_dir = os.path.join( csv_dir, 'xyz_coords' )
    rot_coords_dir = os.path.join( csv_dir, 'rotated_coords' )
    benzene_csv_dir = os.path.join( csv_dir, 'cm_benzene' )
    hydrogen_csv_dir = os.path.join( csv_dir, 'cm_hydrogen' )
    fun_group_csv_dir = os.path.join( csv_dir, 'cm_fun_group' )
    average_coords_csv = os.path.join( csv_dir, 'average_benzene_H2_coords.csv' )
    extended_representation_csv = os.path.join( csv_dir, 'extended_representation.csv' )
    all_cms_csv = os.path.join( csv_dir, 'all_cms_and_aver.csv' )
    all_cms_df = pd.DataFrame()

    os.makedirs( xyz_coords_dir, exist_ok = True )
    os.makedirs( rot_coords_dir, exist_ok = True )
    os.makedirs( benzene_csv_dir, exist_ok = True )
    os.makedirs( hydrogen_csv_dir, exist_ok = True )
    os.makedirs( fun_group_csv_dir, exist_ok = True )

    C0_df, C1_df, C2_df, C3_df, C4_df, C5_df, H0_df, H1_df = 8*[pd.DataFrame()]
    #if os.path.exists( all_cms_csv ):
    #   all_cms_df = pd.read_csv( all_cms_csv )
    #else:
    #   all_cms_df = pd.DataFrame()

    for key, row in dime_results_df.iterrows():
        #fg = row['FG']
        #if fg in tmp_fgs:
        #
           fg = row['FG']
           dime_idx = row['dime.idx']
           mono_idx = row['mono.idx']
           mono_nat = row['mono.nat']
           mp2_en   = row['MP2.EN.']
           mp2_int_en    = row['MP2.INT.EN.']
           mp2_vt_en     = row['MP2.VT.EN.']
           mp2_vt_int_en = row['MP2.VT.INT.EN.']

           fg_label = '{}_{}_{}'.format( fg, mono_idx, dime_idx ) 

           benzene_csv = os.path.join( benzene_csv_dir, 'cm_{}.csv'.format(fg_label) ) 
           hydrogen_csv = os.path.join( hydrogen_csv_dir, 'cm_{}.csv'.format(fg_label) ) 
           fun_group_csv = os.path.join( fun_group_csv_dir, 'cm_{}.csv'.format(fg_label) ) 
           rot_coords_csv = os.path.join( rot_coords_dir, 'rotated_coords_{}.csv'.format(fg_label) ) 
           xyz_coords_dat = os.path.join( xyz_coords_dir, 'rotated_coords_{}.xyz'.format(fg_label) ) 

           ## READ/WRITE ALL ROTATED COORDINATES CSV starts
           if os.path.exists( rot_coords_csv ):
              print( 'File exists: ', rot_coords_csv )
              rotated_coords_df = pd.read_csv( rot_coords_csv, index_col = 0 )
           else:
              print( 'Write file: ', rot_coords_csv )
              coords = row.drop( [ 'FG', 'MP2.EN.', 'MP2.INT.EN.', 'MP2.VT.EN.', 'MP2.VT.INT.EN.',
                                   'NEG.FREQ.', 'VERY.NEG.FREQ.', 'dime.idx', 'mono.idx', 'mono.nat' ] )
              coords.dropna(inplace=True)
 
              dimer_obj = ORCA.DIME( fg, mono_nat, mono_idx, dime_idx )
              dimer_obj.files_names( dime=True, very_tight=True )

              coordinates = []
              coord_array = []
              print( 'Reading' , dimer_obj.xyz_file )

              with open(dimer_obj.xyz_file, 'r') as f:
                   coords_lines = f.readlines()
              for c_count, c_line in enumerate(coords_lines[2:]):
                  coordinates.append( [c_count] + c_line.split()) 
                  coord_array.append( [ float(cc) for cc in c_line.split()[1:] ] )
              coord_array = np.array( coord_array )
     
              ## recognize benzene, H2 and fun. group and rotate
              cc_idxs, hh_idxs, H2_idxs, fg_idxs, [C0_idx, C1_idx, C2_idx] = dimer_obj.find_benzene(coordinates)
              print( 'C/benzene idx = {}, H/benzene idx = {}, H2 idx = {}, FG idx = {}'.format( 
                                      cc_idxs, hh_idxs, H2_idxs, fg_idxs ))
              print( 'Rotating wrt: {}'.format( [C0_idx, C1_idx, C2_idx] ))
          
              rotated_coordinates = rotate_molecule( coord_array, C0_idx, C1_idx, C2_idx )
   
              ## sorting non fixed Cs in benzene according to y coordinate
              non_fixed_ccs_idxs = list( set(cc_idxs) - set( [C0_idx, C1_idx, C2_idx] ) ) 
              non_fixed_ccs = rotated_coordinates[ non_fixed_ccs_idxs ]
              non_fixed_ccs = non_fixed_ccs[ non_fixed_ccs[:,1].argsort() ]
              rotated_coordinates[ non_fixed_ccs_idxs ] = non_fixed_ccs 
         
              ## write rotated result to csv file 
              sorted_list_of_index  = [C0_idx, C1_idx, C2_idx] +  \
                 list( set(cc_idxs) - set( [C0_idx, C1_idx, C2_idx] ) ) +  \
                 hh_idxs + fg_idxs + H2_idxs 
              sorted_list_of_groups = len(cc_idxs)*['benz.ring'] + len(hh_idxs)*['benz.ring'] + \
                                      len(fg_idxs)*['fun.group'] + 2*['h2.mol']
              rotated_coords_df = pd.DataFrame()
              for idx, group in zip( sorted_list_of_index, sorted_list_of_groups ): 
                  symbol = coordinates[idx][1]
                  idx_dict = initiate_elem_dict( symbol, idx, coordinates, rotated_coordinates, group )
                  rotated_coords_df = rotated_coords_df.append( [idx_dict], ignore_index=True )
              rotated_coords_df.to_csv( rot_coords_csv )

              ## write xyz coordinates to file begins
              print( 'Writing file: ', xyz_coords_dat )
              with open( xyz_coords_dat, 'w+' ) as xyz:
                 xyz.write( '{}\n'.format(len(coordinates)) )
                 xyz.write( 'MOF - FG {}\n'.format(fg_label)) 
                 for row, line in rotated_coords_df.iterrows():
                     xyz.write( '{} {} {} {}\n'.format( line['element'], line['xx'], line['yy'], line['zz'] ))
              ## write xyz coordinates to file ends

           ## READ/WRITE ALL ROTATED COORDINATES CSV ends


           ## CALCULATE AVERAGE ELECTRONEGATIVITY/CENTER OF MASS FOR BENZENE/H2/FUNCTIONAL GROUP begins 
           for label, csv_file in zip( [ 'BENZ.RING', 'H2.MOL', 'FUN.GROUP' ],
                                       [ benzene_csv,  hydrogen_csv, fun_group_csv ] ):
               if os.path.exists( csv_file ):
                  print( 'File exists: ', csv_file )
               else:
                  print( 'Write file: ', csv_file )
                  df_block = rotated_coords_df.loc[ ( rotated_coords_df[label] == 1 ) ]

                  ion_ener = df_block['ioniz_energy'].mean()
                  elec_aff = df_block['elec_affinity'].mean()
                  at_weight = df_block['atomic_weight'] 
                  xx = df_block['xx'] 
                  yy = df_block['yy'] 
                  zz = df_block['zz'] 

                  cm_xx = (at_weight * xx).sum() / at_weight.sum()  
                  cm_yy = (at_weight * yy).sum() / at_weight.sum()  
                  cm_zz = (at_weight * zz).sum() / at_weight.sum()  

                  cm_df = pd.DataFrame( [ {'cm_xx' : cm_xx,  'cm_yy' : cm_yy,  'cm_zz' : cm_zz, 'aver_elec_affinity' : elec_aff, 'aver_ion_energy' : ion_ener } ] )
                  cm_df.to_csv( csv_file )
           ## CALCULATE AVERAGE ELECTRONEGATIVITY/CENTER OF MASS FOR BENZENE/H2/FUNCTIONAL GROUP ends 

           ## GATHER ALL AVERAGE FG/H2/BENZ INFO begins
           key_average_dict = { 'FG' : fg, 'mono.idx' : mono_idx, 'dime.idx' : dime_idx, 'MP2.VT.INT.EN.' : mp2_vt_int_en }
           for tmp_csv, tmp_label in zip( [ benzene_csv, hydrogen_csv, fun_group_csv ], ['BENZ', 'H2', 'FG'] ):
             tmp_dict = pd.read_csv( tmp_csv, index_col=0 ).to_dict('index')[0]
             for k, v in tmp_dict.items():
                 tmp_k = '{}.{}'.format(tmp_label, k)
                 key_average_dict[tmp_k] = v
           all_cms_df = all_cms_df.append( [ key_average_dict ], ignore_index = True )
           ## GATHER ALL AVERAGE FG/H2/BENZ INFO ends

           ## CALCULATE BENZENE/H2 AVERAGE POSITIONS begins 
           benzene_cc_df = rotated_coords_df.loc[ ( rotated_coords_df['BENZ.RING'] == 1 ) & ( rotated_coords_df['element'] == 'C' ) ] 
           hydrogen_cc_df = rotated_coords_df.loc[ rotated_coords_df['H2.MOL'] == 1 ] 

           C0_df = C0_df.append( benzene_cc_df.loc[0][['xx','yy','zz']], ignore_index = True ) 
           C1_df = C1_df.append( benzene_cc_df.loc[1][['xx','yy','zz']], ignore_index = True ) 
           C2_df = C2_df.append( benzene_cc_df.loc[2][['xx','yy','zz']], ignore_index = True ) 
           C3_df = C3_df.append( benzene_cc_df.loc[3][['xx','yy','zz']], ignore_index = True ) 
           C4_df = C4_df.append( benzene_cc_df.loc[4][['xx','yy','zz']], ignore_index = True ) 
           C5_df = C5_df.append( benzene_cc_df.loc[5][['xx','yy','zz']], ignore_index = True ) 
        
           H0_df = H0_df.append( hydrogen_cc_df.loc[hydrogen_cc_df.index[0]][['xx','yy','zz']], ignore_index = True ) 
           H1_df = H1_df.append( hydrogen_cc_df.loc[hydrogen_cc_df.index[1]][['xx','yy','zz']], ignore_index = True ) 
           ## CALCULATE BENZENE/H2 AVERAGE POSITIONS ends 


    average_coords_df = pd.DataFrame( [ { 
                                          'C0.xx.mean' : C0_df['xx'].mean(), 'C0.xx.std' : C0_df['xx'].std() ,  
                                          'C0.yy.mean' : C0_df['yy'].mean(), 'C0.yy.std' : C0_df['yy'].std() ,  
                                          'C0.zz.mean' : C0_df['zz'].mean(), 'C0.zz.std' : C0_df['zz'].std() ,  
                                          'C1.xx.mean' : C1_df['xx'].mean(), 'C1.xx.std' : C1_df['xx'].std() ,  
                                          'C1.yy.mean' : C1_df['yy'].mean(), 'C1.yy.std' : C1_df['yy'].std() ,  
                                          'C1.zz.mean' : C1_df['zz'].mean(), 'C1.zz.std' : C1_df['zz'].std() ,  
                                          'C2.xx.mean' : C2_df['xx'].mean(), 'C2.xx.std' : C2_df['xx'].std() ,  
                                          'C2.yy.mean' : C2_df['yy'].mean(), 'C2.yy.std' : C2_df['yy'].std() ,  
                                          'C2.zz.mean' : C2_df['zz'].mean(), 'C2.zz.std' : C2_df['zz'].std() ,  
                                          'C3.xx.mean' : C3_df['xx'].mean(), 'C3.xx.std' : C3_df['xx'].std() ,  
                                          'C3.yy.mean' : C3_df['yy'].mean(), 'C3.yy.std' : C3_df['yy'].std() ,  
                                          'C3.zz.mean' : C3_df['zz'].mean(), 'C3.zz.std' : C3_df['zz'].std() ,  
                                          'C4.xx.mean' : C4_df['xx'].mean(), 'C4.xx.std' : C4_df['xx'].std() ,  
                                          'C4.yy.mean' : C4_df['yy'].mean(), 'C4.yy.std' : C4_df['yy'].std() ,  
                                          'C4.zz.mean' : C4_df['zz'].mean(), 'C4.zz.std' : C4_df['zz'].std() ,  
                                          'C5.xx.mean' : C5_df['xx'].mean(), 'C5.xx.std' : C5_df['xx'].std() ,  
                                          'C5.yy.mean' : C5_df['yy'].mean(), 'C5.yy.std' : C5_df['yy'].std() ,  
                                          'C5.zz.mean' : C5_df['zz'].mean(), 'C5.zz.std' : C5_df['zz'].std() ,  
                                          'H0.xx.mean' : H0_df['xx'].mean(), 'H0.xx.std' : H0_df['xx'].std() ,  
                                          'H0.yy.mean' : H0_df['yy'].mean(), 'H0.yy.std' : H0_df['yy'].std() ,  
                                          'H0.zz.mean' : H0_df['zz'].mean(), 'H0.zz.std' : H0_df['zz'].std() ,  
                                          'H1.xx.mean' : H1_df['xx'].mean(), 'H1.xx.std' : H1_df['xx'].std() ,  
                                          'H1.yy.mean' : H1_df['yy'].mean(), 'H1.yy.std' : H1_df['yy'].std() ,  
                                          'H1.zz.mean' : H1_df['zz'].mean(), 'H1.zz.std' : H1_df['zz'].std() 
                                         } ] )  

    average_coords_df.to_csv( average_coords_csv )
    all_cms_df.to_csv( all_cms_csv )


    if os.path.exists( extended_representation_csv ):
       pass
    else:
      ext_repr_df = pd.DataFrame()
      
      for row, line in dime_results_df.iterrows():
          fg = line['FG']
          if fg not in ['01_H']:
          #if fg in sample_fgs:
              
              dime = line['dime.idx']
              mono = line['mono.idx']
      
              row_dict = line.to_dict() 
              
              label = '{}_{}_{}.csv'.format(fg, int(mono), int(dime))
          
              cm_hydrogen_csv_file = os.path.join( csv_dir, 'cm_hydrogen', 'cm_{}'.format(label) )
              cm_hydrogen_df = pd.read_csv( cm_hydrogen_csv_file, index_col = 0 )
              
              cm_fun_group_csv_file = os.path.join( csv_dir, 'cm_fun_group', 'cm_{}'.format(label) )
              cm_fun_group_df = pd.read_csv( cm_fun_group_csv_file, index_col = 0 )
              
              rotated_csv_file = os.path.join( csv_dir, 'rotated_coords', 'rotated_coords_{}'.format(label) )
              rotated_df = pd.read_csv( rotated_csv_file, index_col = 0 )
      
              for column in ['cm_xx', 'cm_yy', 'cm_zz']:
                  row_dict['H2_{}'.format(column)] = cm_hydrogen_df[column].values[0]
                  row_dict['FG_{}'.format(column)] = cm_fun_group_df[column].values[0]
      
              row_dict['FG_aver_elec_affinity'] = cm_fun_group_df['aver_elec_affinity'].values[0]
              row_dict['FG_aver_ion_energy'] = cm_fun_group_df['aver_ion_energy'].values[0]
              
              fun_group_coords = rotated_df.loc[ rotated_df['FUN.GROUP'] == 1 ][
                  ['element', 'atomic_weight', 'elec_affinity', 'ioniz_energy', 
                  'group', 'period', 'vdw_radius', 'xx', 'yy', 'zz' ]]
              fun_group_coords.sort_values( by = 'elec_affinity', ascending=False, inplace = True )
              fun_group_coords.reset_index( drop = True, inplace = True )
      
              for row, line in fun_group_coords.iterrows():
                  elem_dict = line.to_dict()
                  del elem_dict['element']
                  for k,v in elem_dict.items():
                      row_dict['fg.{}.{}'.format(row, k)] = v
              
              ext_repr_df = ext_repr_df.append( [row_dict], ignore_index=True )
      
      ext_repr_df = ext_repr_df.fillna(0)
      
      ext_repr_df.to_csv(extended_representation_csv)

if __name__ == '__main__':
  main()
