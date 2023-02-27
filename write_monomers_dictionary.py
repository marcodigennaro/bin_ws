#!/home/mdi0316/anaconda3/bin/python

### common input start
import os, sys, re
import numpy as np
import pandas as pd
import shutil
import subprocess as sp
import datetime

scripts_dir = '/home/mdi0316/FUNCTIONS'
classes_dir = '/home/mdi0316/CLASSES'
sys.path.insert(0, scripts_dir)
sys.path.insert(0, classes_dir)

import IONIC_LIQUID
import GAMESS
#from Regression import pred_GPR_3D, pred_KRR, get_dataset
#from Functions  import *
#from Plot       import *

from Functions import running_jobs

from IONIC_LIQUID import mono_dir

from GAMESS import functionals_list, gbasis_list

run_ids, run_job_labels = running_jobs()

mono_dict = {}

mono_dict['EMIM'] =  { 'nat':19, 'type':'cation', 'input':'emim.inp',
                       'composition' : {'C':6, 'N':2, 'H':11, 'p-atoms' : 8, 'd-atoms'  : 0} }
mono_dict['BMIM'] =  { 'nat':25, 'type':'cation', 'input':'bmim.inp' }
mono_dict['BF4']  =  { 'nat': 5, 'type':'anion',  'input':'bf4.inp', 
                       'composition' : {'B':1, 'F':4, 'p-atoms' : 5, 'd-atoms': 0} }
mono_dict['SCN']  =  { 'nat': 3, 'type':'anion',  'input':'scn.inp' }
mono_dict['DEP']  =  { 'nat':19, 'type':'anion',  'input':'dep.inp' }
mono_dict['PF6']  =  { 'nat': 7, 'type':'anion',  'input':'pf6.inp', 
                       'composition' : {'P':1, 'F':5, 'p-atoms': 7, 'd-atoms': 0} }

from bs_conv_monomers import post_dict

for mono in mono_dict.keys():
    mono_dict[mono]['OUT'] = {}
    natoms = mono_dict[mono]['nat']
    root_dir = os.path.join( mono_dir, mono ) 
    #for gb in ['APCseg-1']: #gbasis_list:
    for gb in gbasis_list:
        mono_dict[mono]['OUT'][gb] = {}
        for fun in functionals_list:
            mono_dict[mono]['OUT'][gb][fun] = {}
            
            for k,v in post_dict.items():
                tmp_run  = v['runtyp']
                run_lab = tmp_run[:3]

                tmp_post = v['post_typ']
                if tmp_post == 'DFTTYP':
                   post_lab = 'DFT'
                else:
                   post_lab = tmp_post

                OUT_dir = os.path.join( root_dir, post_lab, gb, fun, run_lab, 'ZMT' ) 
                mono_dict[mono]['OUT'][gb][fun][post_lab] = {}
                mono_dict[mono]['OUT'][gb][fun][post_lab][run_lab] = {}

                if os.path.exists( OUT_dir ):
                  inp_file_name = [ ff for ff in  os.listdir( OUT_dir ) if ff.endswith( 'inp' ) ][0]
                  inp_file_lab = inp_file_name.replace('.inp', '')
                  mono_dict[mono]['OUT'][gb][fun][post_lab][run_lab] = { 'inp' : inp_file_name }
                  if inp_file_name in run_job_labels:
                     mono_dict[mono]['OUT'][gb][fun][post_lab][run_lab]['err'] = 'Running'
                  else:
                     logs_list = [ ff for ff in  os.listdir( OUT_dir ) if ff.startswith( 'log' ) ]
                     if len(logs_list) == 1:
                        mono_obj = GAMESS.GAMESS_calculation( inp_file_lab,
                                                              root_dir, zero_en = 0, icharge = 0,
                                                              runtyp=tmp_run, post_scf=post_lab, basis=gb, functional=fun,
                                                              natoms=natoms )
                        mono_exec = mono_obj.get_execution(run_job_labels)
                        if mono_exec == 'NORMALLY':
                           
                           if run_lab == 'OPT':
                              out_dict = mono_obj.get_out_dict()[post_lab][run_lab]['FINAL']
                           else:
                              out_dict = mono_obj.get_out_dict()[post_lab][run_lab]
                           for kk, vv in out_dict.items():
                               if kk != 'NSERCH':
                                 mono_dict[mono]['OUT'][gb][fun][post_lab][run_lab][kk] = vv
                        elif mono_exec == '-ABNORMALLY-':
                           mono_dict[mono]['OUT'][gb][fun][post_lab][run_lab]['ERR'] = mono_obj.read_error()
                           print( mono_obj.run_dir )
                           print( mono_obj.read_error() )
                        else:
                           out_dict = mono_obj.get_out_dict()[post_lab][run_lab]
                           print( mono_obj.run_dir, mono_obj.read_error(), out_dict.keys() )
                           

                if mono_dict[mono]['OUT'][gb][fun][post_lab] == {}:
                   del mono_dict[mono]['OUT'][gb][fun][post_lab]
            if mono_dict[mono]['OUT'][gb][fun] == {}:
               del mono_dict[mono]['OUT'][gb][fun]
        if mono_dict[mono]['OUT'][gb] == {}:
           del mono_dict[mono]['OUT'][gb]
        
            
    
import json
with open('monomers.json', 'w') as fp:
    json.dump(mono_dict, fp)

shutil.copy2( 'monomers.json', '/home/mdi0316/Inputfiles/GAMESS' )


