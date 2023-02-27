#!/home/mdi0316/anaconda3/bin/python

### common input start
import os, sys, re
import numpy as np
import pandas as pd
import shutil
import subprocess as sp
import datetime

import getpass
user = getpass.getuser()

scripts_dir = '/home/mdi0316/FUNCTIONS'
classes_dir = '/home/mdi0316/CLASSES'
sys.path.insert(0, scripts_dir)
sys.path.insert(0, classes_dir)

import IONIC_LIQUID as IL
import GAMESS
from Regression import pred_GPR_3D, pred_KRR, get_dataset
from Functions  import *

import json

from qml.fchl import generate_representation

csv_file = 'cart_coords.csv'
df = pd.read_csv(csv_file, index_col = 0, dtype=object)

for idx, row in df.iterrows():
    cart_coords = ast.literal_eval(row['cart.coords.'])
    mull_charges = ast.literal_eval(row['mull.charges'])
    
    at_charg = np.array([float(v['charge']) for v in mull_charges.values() ])
    at_coord = np.empty([0,3])
    for kk, tmp_dict in cart_coords.items():
        x = float(tmp_dict['x']) 
        y = float(tmp_dict['y']) 
        z = float(tmp_dict['z'])
        print( kk, tmp_dict, [x,y,z] )
        at_coord = np.vstack( ( at_coord, [x,y,z] ) )
        rep = generate_representation( at_coord, at_charg[0:int(kk)+1] )
        print( kk, at_coord, len(at_charg[0:int(kk)+1]) )
    
    print( at_coord.shape ) 
    print( at_charg.shape ) 
    rep = generate_representation( at_coord, at_charg )
    print(rep)
