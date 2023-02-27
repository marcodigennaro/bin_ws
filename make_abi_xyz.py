#!/home/mdi0316/anaconda3/bin/python
import os, sys, re
import shutil
import numpy as np
from numpy import linalg as LA
import pandas as pd
import csv 
import time

import matplotlib as mlp
import matplotlib.pyplot as plt

main_dir = '/home/mdi0316/WORK/SEP/ABINIT/EMIM_BF4'
scripts_dir = '/home/mdi0316/FUNCTIONS'
sys.path.insert(0, scripts_dir)

from Functions import print_tab, rotate_molecule, center_of_mass

xyz_file = 'emim_bf4.xyz'
xyz_lines = open(xyz_file, 'r').readlines()


#xyz_np = np.array( [ np.array( [ float(xx) for xx in line.split()[1:]] ) for line in xyz_lines[2:] ] )
cxyz_np = np.array( [ np.array( line.split() ) for line in xyz_lines[2:] ] )
elem_list = cxyz_np[:,0]

emim_np = cxyz_np[:19, 1:].astype(float)
bf4_np  = cxyz_np[19:, 1:].astype(float)

emim_np=rotate_molecule(emim_np)
bf4_np=rotate_molecule(bf4_np)

emim_dict = {}
for k,v in enumerate(emim_np):
   emim_dict[k] = {}
   emim_dict[k]['elem.'] = 'C'
   emim_dict[k]['x'] = v[0]
   emim_dict[k]['y'] = v[1]
   emim_dict[k]['z'] = v[2]

emim_com = center_of_mass( emim_dict )
bf4_np += emim_com

def make_fig(np_array):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #ax.scatter(xyz_np[:,0], xyz_np[:,1], xyz_np[:,2], label='orig')
    ax.scatter(np_array[:,0], np_array[:,1], np_array[:,2], label='emim')
    plt.legend()
    plt.show()


for Z in range( 10, 20 ):

   bf4_np[:,2] += Z
   emim_bf4_np = np.vstack( [emim_np, bf4_np] )
   #make_fig( emim_bf4_np )

   new_dir = os.path.join( main_dir, 'SCAN', 'Z_{}'.format(Z) )
   os.makedirs( new_dir, exist_ok=True )
   new_xyz_file = os.path.join( new_dir, 'emim_bf4.xyz' )
   with open( new_xyz_file, 'w+' ) as f:
      f.write( '24\nemim_bf4\n' )
      for elem, coord in zip(elem_list, emim_bf4_np):
          f.write( '{}  {:8.8f}  {:8.8f}  {:8.8f}\n'.format(elem, coord[0], coord[1], coord[2]) )

   
