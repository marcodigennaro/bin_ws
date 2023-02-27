#!/home/mdi0316/anaconda3/bin/python

### common input start
import time
import os, sys, re
import shutil
import numpy as np
from numpy import linalg as LA
import pandas as pd
import subprocess as sp
import csv 
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

from Functions import print_tab, running_jobs, compose_zmatrices, running_label, center_of_charge, center_of_mass, Coulomb_Energy, angle_between
 
from make_scan_list import *
from dimers import read_dimer

def find_nearest(array, value):
    array = np.asarray(array)
    value = float(value)
    print( np.where(array == array.min()) )
    if value in array:
       array = np.array(list(filter(lambda x : x !=value, array)))
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def neighbors(arr,x,y,n=3):
    ''' Given a 2D-array, returns an nxn array whose "center" element is arr[x,y]'''
    arr=np.roll(np.roll(arr,shift=-x+1,axis=0),shift=-y+1,axis=1)
    return arr[:n,:n]

def extract_line( read_df, radius, theta, phi ):
    return read_df.loc[ read_df['Radius'] == radius ].loc[ read_df['Theta'] == int(theta) ].loc[ read_df['Phi'] == int(phi) ]

def main():

    for read_scf, read_col in zip( [ 'dft', 'mp2'], [ 'DFT.INT.EN.', 'MP2.INT.EN.' ] ):

        all_csvs = sp.check_output( 'ls DIMERS/*/CSV/N311/B3LYP/SCAN/scan_{}.csv'.format(read_scf), shell=True ).decode("utf-8").split('\n')
   
        all_csvs.remove('')
        for read_csv in all_csvs:

            print_tab(1, 'read file: {}'.format(read_csv))
            read_df = pd.read_csv( read_csv, index_col = 0 )
            min_ener_path_df = pd.DataFrame()

            print_tab(1, 'read file has {} lines'.format(len(read_df)))
            print_tab(0, 'read file columns: ')
            print_tab(0, read_df.dtypes)

            all_radius = list(set(read_df['Radius']))
            all_thetas = list(set(read_df['Theta']))
            all_phis   = list(set(read_df['Phi']))

            all_radius.sort(reverse=True)
            all_thetas.sort()
            all_phis.sort()

            print_tab( 1, 'all R = {}'.format(all_radius) )
            print_tab( 1, 'all T = {}'.format(all_thetas) )
            print_tab( 1, 'all P = {}'.format(all_phis  ) )

            start_r = max(all_radius)
            
            start_directions = [ (  '5',  '0'),
                                 ('175',  '0'), 
                                 ( '90',  '0'), 
                                 ( '90', '90'), 
                                 ( '90','180'), 
                                 ( '90','270')  ]

            ## create T/P matrix  
            full_mat = np.empty( [len(all_thetas), len(all_phis), 2 ] )
            for tdx, t in enumerate(all_thetas):
              for pdx, p in enumerate(all_phis):
                full_mat[tdx, pdx] = [int(t),int(p)]

            for init_t, init_p in start_directions:
            
                print_tab(2, 'new initial direction: T={}, P={}'.format(init_t, init_p))

                center_t, center_p = init_t, init_p
  
                for center_r in all_radius:
   
                    center_line = extract_line( read_df, center_r, int(center_t), int(center_p) ) 
                    if not  center_line.empty:
                       center_ener = center_line[read_col].values[0]
                       print_tab( 3, 'scan around T={}, P={}, R={}, E={}'.format( center_t, center_p, center_r, center_ener ))
  
                       ### get reduced matrix from full T/P matrix 3x3 around center
                       #for tdx, tmp_t in enumerate(all_thetas):
                       #  for pdx, tmp_p in enumerate(all_phis):
                       #    if int(tmp_t) == int(center_t) and int(tmp_p) == int(center_p):
                       #      reduced_mat = neighbors( full_mat, tdx, pdx ).reshape(9,2) 
                       #      break

                       ## get reduced matrix from full T/P matrix only 4 (up/down, left/right)
                       tdx =  all_thetas.index(int(center_t))
                       pdx =  all_phis.index(  int(center_p))
                       prev_t = np.roll( all_thetas, -tdx )[-1]
                       next_t = np.roll( all_thetas, -tdx )[1]
                       prev_p = np.roll( all_phis, -pdx )[-1]
                       next_p = np.roll( all_phis, -pdx )[1]
                       reduced_mat = [ (t,p) for t in [prev_t, next_t] for p in [prev_p, next_p] ]


                       print_tab( 3, 'sourrounding matrix: {}'.format(reduced_mat) ) #.flatten()) )
  
                       min_E = 1e6
                       min_T, min_P = None, None
                       for (close_t, close_p) in reduced_mat:
                           close_line = extract_line( read_df, center_r, int(close_t), int(close_p) ) 
                           if close_line.empty:
                              close_ener = 'result not finished' 
                           else:  
                              close_ener = close_line[read_col].values[0]
                              if close_ener < -10:
                                 print( close_line )
                                 exit() 
                              if close_ener < min_E:
                                 #print_tab( 3, 'new minimum found')
                                 min_E = close_ener
                                 min_T, min_P = close_t, close_p
                                 
                           #print_tab( 2, 'T={}, P={}, R={}, E={}'.format( close_t, close_p, center_r, close_ener ))
  
                       print_tab( 3, 'min energy point T={}, P={}, R={}, (E={})'.format( center_t, center_p, center_r, center_ener ))
                       min_ener_serie = pd.Series( { 'init.T' : init_t, 'init.P' : init_p, 'R' : center_r, 
                                                     'min.T'  : min_T,  'min.P'  : min_P,  'min.{}'.format(read_col) : min_E } )
                       min_ener_path_df = min_ener_path_df.append( min_ener_serie, ignore_index = True )
  
                       center_t, center_p = min_T, min_P 
  
            min_ener_csv = read_csv.replace( 'scan', 'min_path' )
            min_ener_path_df.to_csv( min_ener_csv )

if __name__ == '__main__':
   main()
