#!/home/mdi0316/anaconda3/bin/python

import os, sys, re
import shutil
import numpy as np
import pandas as pd

user = 'mdi0316'
classes_dir = '/home/{}/CLASSES'.format(user)

sys.path.insert(0, classes_dir)

import GAMESS

minimal_T_SCAN_LIST = [ '90' ]
minimal_P_SCAN_LIST = [ '90' ]
minimal_P_SCAN_LIST = [ '7.0' ]

limited_T_SCAN_LIST = [ '5', '90', '175' ]
limited_P_SCAN_LIST = [ '0', '90', '180', '270' ]
limited_R_SCAN_LIST = [ '2.0', '2.5', '3.0', '3.5',
                        '4.0', '4.5', '5.0', '5.5', 
                        '6.0', '6.5', '7.0', '7.5', 
                        '8.0', '8.5', '9.0', '9.5', 
                        '10.0', '11.0', '12.0', '13.0', '15.0', '17.0', '20.0', '25.0', '30.0'  ]

extended_T_SCAN_LIST = limited_T_SCAN_LIST + [ '45', '135' ]
extended_P_SCAN_LIST = limited_P_SCAN_LIST + [ '45', '135', '225', '315' ]
extended_R_SCAN_LIST = limited_R_SCAN_LIST + [ '2.3',
                                               '2.6', '2.7', '2.8', '2.9', 
                                               '3.1', '3.2', '3.3', '3.4', '3.6', '3.7', '3.8', '3.9',
                                               '4.1', '4.2', '4.3', '4.4', '4.6', '4.7', '4.8', '4.9' ]

#extended_T_SCAN_LIST = [ '90' ]
#extended_P_SCAN_LIST = [ '90' ]
#extended_P_SCAN_LIST = [ '7.0' ]

minimal_basis_list = ['N311'] 
minimal_funct_list = ['B3LYP']

limited_basis_list = minimal_basis_list + [ 'APCseg-1', 'STO', 'CCQ' ]
limited_funct_list = minimal_funct_list + [ 'PBE0', 'M11' , 'wB97x-D' ]

#limited_basis_list = minimal_basis_list 
#limited_funct_list = minimal_funct_list 

extended_basis_list = GAMESS.full_gbasis_list
extended_funct_list = GAMESS.full_functionals_list

def get_scan_list( label, basis, funct ):
    t_list = limited_T_SCAN_LIST
    p_list = limited_P_SCAN_LIST
    r_list = limited_R_SCAN_LIST
    if [basis, funct] == ['N311', 'B3LYP']:
       t_list = extended_T_SCAN_LIST
       p_list = extended_P_SCAN_LIST
       r_list = extended_R_SCAN_LIST
    else:
       t_list = ['90']
       p_list = ['90'] #limited_P_SCAN_LIST
       r_list = limited_R_SCAN_LIST

    for tmp_list in [t_list, r_list, p_list]:
      try:
        tmp_list = [ int(x) for x in tmp_list ]
      except(ValueError):
        tmp_list = [ float(x) for x in tmp_list ]
      tmp_list.sort()
      tmp_list = [ str(x) for x in tmp_list ]
    return( r_list, t_list, p_list )
