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
import glob

import getpass
user = getpass.getuser()

scripts_dir = '/home/{}/FUNCTIONS'.format(user)
classes_dir = '/home/{}/CLASSES'.format(user)
zmat_converter_dir = '/home/{}/CLASSES/zmatrix-master'.format(user)

sys.path.insert(0, scripts_dir)
sys.path.insert(0, classes_dir)
sys.path.insert(0, zmat_converter_dir)


from Functions import print_tab, running_jobs, compose_zmatrices, running_label, center_of_charge, center_of_mass, Coulomb_Energy, angle_between

dime_path='/data/mdi0316/WORK/DIMERS'
for r, d, f in os.walk(dime_path):
    if r.endswith('DFT' ):
       if 'FAILED' in os.listdir(r):
         try:
            inp_file = [ f for f in os.listdir(r) if f.endswith('inp') ][0]
            running = running_label( inp_file )
            located = False
            not_loc = False
            aborted = False
            abnorma = False
            if not running:
               log_file = [ f for f in os.listdir(r) if f.startswith('log') ][0]
               log_lines = open( os.path.join( r, log_file ), 'r' ).readlines()
               log_lines.reverse()
               for l in log_lines:
                   if 'LOCATED' in l:
                      located = True
                      break
                   if 'FAIL' in l:
                      not_loc = True
                      print( 'FAIL',  r )
                      break
                   if 'ABNORMALLY' in l:
                      aborted = True
                      print( 'ABNOR',  r )
                      for ll in log_lines:
                           if 'ERROR OPENING' in ll:
                              print( ll, r )
                              break
                      break
                   if 'ABORTED' in l:
                      aborted = True
                      print( 'ABORT',  r )
                      break
         except(IndexError):
            pass

 
